# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 15:34:47 2023

@author: YuxiaoLan
"""
from HTMACat.Extract_info import *
import os
import numpy as np
import math
import ase
from ase.build import bulk
from ase.visualize import view
from ase.io.vasp import write_vasp
from catkit.gen.surface import SlabGenerator
from catkit.build import surface
from catkit.gen.adsorption import Builder
from catkit.gratoms import *
from catkit.gen.adsorption import AdsorptionSites


class Bulk(object):
    def __init__(self, main_element='Pt', lattice_type='fcc', lattice_constant=None,
                 ele_dop='Cu', natom_dop='0', super_cell=None):
        if super_cell is None:
            super_cell = [3, 3, 1]
        if lattice_constant is None:
            lattice_constant = {'a': 3.96}
        if isinstance(natom_dop, int):
            natom_dop = str(natom_dop)
        if not isinstance(lattice_constant, dict):
            self.lattice_constant = {}
            self.set_lattice_constant(lattice_constant)
        else:
            self.lattice_constant = lattice_constant  # 不同的lattice_type有不同的constant，hcp有个两个constant

        self.main_element = main_element
        self.lattice_type = lattice_type
        self.ele_dop = ele_dop
        self.natom_dop = natom_dop
        if self.natom_dop[0] == 'b':
            self.super_cell = [2, 2, 1]
        else:
            self.super_cell = super_cell

    def set_lattice_constant(self, latcont):
        if not isinstance(latcont, list):
            self.lattice_constant['a'] = latcont
        else:
            self.lattice_constant['a'] = latcont[0]
            self.lattice_constant['c'] = latcont[1]

    def get_ele_dop(self):
        return self.ele_dop

    def get_super_cell(self):
        return self.super_cell

    def get_natom_dop(self):
        return self.natom_dop

    def get_main_element(self):
        return self.main_element

    def construct(self):
        mname = self.main_element
        mlat = self.lattice_type
        if self.lattice_type == 'hcp':
            mlatcon = self.lattice_constant['a']
            mc = self.lattice_constant['c']
            mbulk = bulk(mname, mlat, a=float(mlatcon), covera=float(mc) / float(mlatcon), cubic=False)
        else:
            mlatcon = self.lattice_constant['a']
            mbulk = bulk(mname, mlat, a=mlatcon, cubic=True)
        mbulk = self.dop_bulk(mbulk)
        return mbulk

    def dop_bulk(self, mbulk):
        N_dop_bulk = self.get_ele_dop()
        natom = self.get_natom_dop()
        if N_dop_bulk == [] or natom[0] != 'b':
            return mbulk

        dop_atoms_number = int(natom[1])
        assert dop_atoms_number < len(mbulk), \
            'bulk dop element should smaller than atoms number in one cell:%d,' \
            ' but %d was given' % (len(mbulk), dop_atoms_number)
        for i in range(dop_atoms_number):
            mbulk[i].symbol = N_dop_bulk
        return mbulk

    def get_dop_element(self):
        return self.ele_dop

    @classmethod
    def from_dict(cls, init_dict):
        main_element = init_dict['element']
        lattice_type = init_dict['lattype']
        lattice_const = init_dict['latcont']
        ele_dop = init_dict['element_dop']
        natom_dop = init_dict['dop_type']
        return cls(main_element=main_element, lattice_type=lattice_type, lattice_constant=lattice_const,
                   ele_dop=ele_dop, natom_dop=natom_dop)


class Slab(object):
    def __init__(self, in_bulk=Bulk(), facet='100'):
        self.bulk = in_bulk
        self.facet = facet
        self.property = {}

    def get_miller_index(self):
        miller_index = tuple(list(map(int, list(self.facet))))
        return miller_index

    def is_dope(self) -> bool:
        natom = self.bulk.get_natom_dop()
        if not natom.isdigit():
            return False
        if natom == '0':
            return False
        return True

    def get_facet(self):
        return self.facet

    def out_file_name(self):
        mname = self.bulk.get_main_element()
        if self.bulk.get_natom_dop() == '0':
            return '_'.join([mname, self.facet])
        else:
            ele_dop = self.bulk.get_dop_element()
            natom_dop = self.bulk.get_natom_dop()
            return '_'.join([mname, ele_dop, self.facet, natom_dop])

    def get_dis_inter(self):
        latcon = self.bulk.lattice_constant
        if self.facet == '111':
            dis_min = float(latcon['a']) * math.sqrt(2) / 2 * math.sqrt(3) / 2
            dis_max = float(latcon['a']) * math.sqrt(2) / 2 * math.sqrt(3)
            dis_inter = [dis_min, dis_max]
        elif (self.facet == '100') or (self.facet == '0001'):
            dis_min = float(latcon['a']) * 1
            dis_max = float(latcon['a']) * 2
            dis_inter = [dis_min, dis_max]
        else:
            raise ValueError('Do not support facet: %s for generate inter distance' % self.facet)
        return dis_inter

    def construct(self):
        natom = self.bulk.get_natom_dop()
        if natom == '0' or natom[0] == 'b':
            slabs = self.Construct_slab()
        elif natom == '1L':
            slabs = self.Construct_1stLayer_slab()
        elif natom.isdigit():
            slabs = self.Construct_doped_slab()
        else:
            raise ValueError('Do not support %s dop type to construct slab' % natom)
        return slabs

    def Construct_slab(self):
        slab = []
        miller_index = self.get_miller_index()
        mbulk = self.bulk.construct()
        super_cell = self.bulk.get_super_cell()
        ##### generate the surfaces #####   
        gen = SlabGenerator(
            mbulk,
            miller_index=miller_index,
            layers=4,
            fixed=2,
            layer_type='trim',
            vacuum=8,
            standardize_bulk=True)
        terminations = gen.get_unique_terminations()
        for i, t in enumerate(terminations):
            slab += [gen.get_slab(iterm=i) * super_cell]
        return slab

    def Construct_1stLayer_slab(self):
        slabs = self.Construct_slab()
        Ele_dop = self.bulk.get_dop_element()
        slabs_dop = []
        for i, slab in enumerate(slabs):
            surface_atoms = slab.get_surface_atoms()
            slb = slab.copy()
            for j, surf_atom in enumerate(surface_atoms):
                slb[surf_atom].symbol = Ele_dop
            slabs_dop += [slb]
            # view(slb*(2,2,1))
        return slabs_dop

    def Construct_doped_slab(self):
        # generate surface adsorption configuration
        slabs = self.Construct_slab()
        slabs_dop = []
        p1 = []
        p1_symb = []
        Natom = int(self.bulk.get_natom_dop())
        Ele_dop = self.bulk.get_dop_element()
        for i, slab in enumerate(slabs):
            site = AdsorptionSites(slab)
            site_typ = site.get_connectivity()
            topo = site.get_topology()
            # atom_number_dop =[]
            for j in range(len(site_typ)):
                if site_typ[j] == Natom:
                    atom_number_dop = topo[j]
                    p1 += [slab.get_positions()[k] for k in atom_number_dop]

                    slb = slab.copy()
                    symbol = slb.get_chemical_symbols()
                    for k, item in enumerate(atom_number_dop):
                        symbol[item] = Ele_dop
                    slb.set_chemical_symbols(symbol)
                    p1_symb += [slb.get_chemical_symbols()[k] for k in atom_number_dop]
                    slabs_dop += [slb]
                    # view(slb*(2,2,1))
        self.property['p1'] = p1
        self.property['p1_symb'] = p1_symb
        return slabs_dop

    @classmethod
    def from_dict(cls, init_dict):
        in_bulk = Bulk.from_dict(init_dict)
        facet = init_dict['facet']
        return cls(in_bulk, facet)


def substrate_from_input(init_dict):
    substrate = Slab.from_dict(init_dict)
    return substrate
