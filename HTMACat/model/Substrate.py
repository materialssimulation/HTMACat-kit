"""Created on Fri Mar 17 15:34:47 2023.

@author: YuxiaoLan
"""

import math
from ase.build import bulk
from ase.io import read
from ase.constraints import FixAtoms
from HTMACat.catkit.gen.surface import SlabGenerator
from HTMACat.catkit.gen.adsorption import AdsorptionSites
from HTMACat.catkit.gratoms import Gratoms
import numpy as np
from HTMACat.model.Structure import Structure
from HTMACat.Base_tools import *


class Bulk:
    def __init__(
        self,
        element="Pt",
        lattice_type="fcc",
        lattice_constant=None,
        ele_dop="Cu",
        natom_dop="0",
        supercell=None,
    ):
        if supercell is None:
            supercell = [3, 3]
        if lattice_constant is None:
            lattice_constant = {"a": 3.96}
        if isinstance(natom_dop, int):
            natom_dop = str(natom_dop)
        if not isinstance(lattice_constant, dict):
            self.lattice_constant = {}
            self.set_lattice_constant(lattice_constant)
        else:
            self.lattice_constant = lattice_constant  # 不同的lattice_type有不同的constant，hcp有个两个constant

        self.main_element = element
        self.lattice_type = lattice_type
        self.ele_dop = ele_dop
        self.natom_dop = natom_dop
        if self.natom_dop[0] == "b":
            self.supercell = [2, 2, 1]
        else:
            self.supercell = supercell + [1]

    def set_lattice_constant(self, latcont):
        if not isinstance(latcont, list):
            self.lattice_constant["a"] = latcont
        else:
            self.lattice_constant["a"] = latcont[0]
            self.lattice_constant["c"] = latcont[1]

    def get_ele_dop(self):
        return self.ele_dop

    def get_supercell(self):
        return self.supercell

    def get_natom_dop(self):
        return self.natom_dop

    def get_main_element(self):
        return self.main_element

    def construct(self):
        mname = self.main_element
        mlat = self.lattice_type
        if self.lattice_type == "hcp":
            mlatcon = self.lattice_constant["a"]
            mc = self.lattice_constant["c"]
            mbulk = bulk(
                mname, mlat, a=float(mlatcon), covera=float(mc) / float(mlatcon), cubic=False
            )
        else:
            mlatcon = self.lattice_constant["a"]
            mbulk = bulk(mname, mlat, a=mlatcon, cubic=True)
        mbulk = self.dop_bulk(mbulk)
        return mbulk

    def dop_bulk(self, mbulk):
        N_dop_bulk = self.get_ele_dop()
        natom = self.get_natom_dop()
        if N_dop_bulk == [] or natom[0] != "b":
            return mbulk

        dop_atoms_number = int(natom[1])
        assert dop_atoms_number < len(mbulk), (
            "bulk dop element should smaller than atoms number in one cell:%d,"
            " but %d was given" % (len(mbulk), dop_atoms_number)
        )
        for i in range(dop_atoms_number):
            mbulk[i].symbol = N_dop_bulk
        return mbulk

    def get_dop_element(self):
        return self.ele_dop




class Slab(Structure):
    def __init__(self, in_bulk=Bulk(), facet="100", layers=4, layers_relax=2):
        self.bulk = in_bulk
        self.file = None
        self.facet = facet
        self.property = {} #掺杂元素序号以及元素符号p1,p1_symb
        self.layers = layers
        self.layers_relax = layers_relax
        if "p1" not in self.property or "p1_symb" not in self.property:
            self.property["p1"] = []
            self.property["p1_symb"] = []

    def get_miller_index(self):
        miller_index = tuple(list(map(int, list(self.facet))))
        return miller_index

    def is_dope(self) -> bool:
        natom = self.bulk.get_natom_dop()
        if not natom.isdigit():
            return False
        if natom == "0":
            return False
        return True

    def get_facet(self):
        return self.facet

    def get_layers(self):
        return self.layers
    
    def get_layers_relax(self):
        return self.layers_relax

    def out_file_name(self):
        mname = self.bulk.get_main_element()
        if self.bulk.get_natom_dop() == "0":
            return "_".join([mname, self.facet])
        else:
            ele_dop = self.bulk.get_dop_element()
            natom_dop = self.bulk.get_natom_dop()
            return "_".join([mname, ele_dop, self.facet, natom_dop])

    def out_print(self):
        mname = self.bulk.get_main_element()
        if self.bulk.get_natom_dop() == "0":
            return f"{mname} ({self.facet}) substrate"
        else:
            ele_dop = self.bulk.get_dop_element()
            natom_dop = self.bulk.get_natom_dop()
            return f"{ele_dop} {natom_dop} doped {mname} ({self.facet}) substrate"

    def get_dis_inter(self):
        latcon = self.bulk.lattice_constant
        if self.facet == "111":
            dis_min = float(latcon["a"]) * math.sqrt(2) / 2 * math.sqrt(3) / 2
            dis_max = float(latcon["a"]) * math.sqrt(2) / 2 * math.sqrt(3)
            dis_inter = [dis_min, dis_max]
        elif (self.facet == "100") or (self.facet == "0001"):
            dis_min = float(latcon["a"]) * 1
            dis_max = float(latcon["a"]) * 2
            dis_inter = [dis_min, dis_max]
        else:
            raise ValueError("Do not support facet: %s for generate inter distance" % self.facet)
        return dis_inter

    def construct(self):
        natom = self.bulk.get_natom_dop()
        if natom == "0" or natom[0] == "b":
            slabs = self.Construct_slab()
        elif natom == "1L":
            slabs = self.Construct_1stLayer_slab()
        elif natom.isdigit():
            slabs = self.Construct_doped_slab()
        else:
            raise ValueError("Do not support %s dop type to construct slab" % natom)
        return slabs

    def Construct_slab(self):
        slab = []
        miller_index = self.get_miller_index()
        mbulk = self.bulk.construct()
        supercell = self.bulk.get_supercell()
        layers = self.get_layers()
        layers_relax = self.get_layers_relax()
        ##### generate the surfaces #####
        gen = SlabGenerator(
            mbulk,
            miller_index=miller_index,
            layers=layers,
            fixed=layers-layers_relax,
            layer_type="trim",
            vacuum=8,
            standardize_bulk=True,
        )
        terminations = gen.get_unique_terminations()
        for i, t in enumerate(terminations):
            slab += [gen.get_slab(iterm=i) * supercell]
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
        Natom = int(self.bulk.get_natom_dop())
        for i, slab in enumerate(slabs):  # 修改了fcc的(100)面无法生成dope type 3的bug
            site = AdsorptionSites(slab)
            site_typ = site.get_connectivity()
            topo = site.get_topology()
            # atom_number_dop =[]
            if Natom in site_typ:
                j = np.argwhere(site_typ == Natom)[0][0]
                atom_number_dop = topo[j]
            elif Natom < len(topo[-1]):
                atom_number_dop = topo[-1]
            else:
                raise ValueError("the max dope number for this system is %d" % len(topo[-1]))
            slb = self.dope_slab(slab, atom_number_dop[0:Natom])
            slabs_dop += [slb]
        return slabs_dop

    def dope_slab(self, slab, atom_number_dop):
        self.property["p1"] += [slab.get_positions()[k] for k in atom_number_dop]
        Ele_dop = self.bulk.get_dop_element()
        slb = slab.copy()
        symbol = slb.get_chemical_symbols()
        for k, item in enumerate(atom_number_dop):
            symbol[item] = Ele_dop
        slb.set_chemical_symbols(symbol)
        self.property["p1_symb"] += [slb.get_chemical_symbols()[k] for k in atom_number_dop]
        return slb

    @classmethod
    def init_one_slab(cls, init_dict):
        # initialize bulk
        bulk_list = ['element','lattice_type','lattice_constant','ele_dop','natom_dop','supercell']
        bulk_init_dict = get_new_dict(bulk_list, init_dict)
        in_bulk = Bulk(**bulk_init_dict)
        # get the parameters for slab initialization
        slab_list = ['facet','layers','layers_relax']
        slab_init_dict = get_new_dict(slab_list, init_dict)

        return cls(in_bulk, **slab_init_dict)

    @classmethod
    def init_all_slab(cls, struct_Info: dict):
        assert isinstance(
            struct_Info, dict
        ), "Substrates init by self defined struct should be dict!"
        if struct_Info == {}:
            return []
        substrates = []
        # get the parameters for struct initialization
        struct_list = ['element','lattice_type','lattice_constant','supercell','layers','layers_relax']
        struct_init_dict = get_new_dict(struct_list,struct_Info)
        # get the parameters for dope and surface initialization
        dope_system = struct_Info["dope"]
        dope_init_list = []
        surface_init_list = []
        for key, value in dope_system.items():
            for i in value:
                dope_init_list.append({"ele_dop": key, "natom_dop": i})
        for i in struct_Info["facet"]:
            surface_init_list.append({"facet": i})

        # substrates initialization
        for i in range(len(dope_init_list)):
            for j in range(len(surface_init_list)):
                init_dict = {**struct_init_dict, **dope_init_list[i], **surface_init_list[j]}
                substrates.append(cls.init_one_slab(init_dict))

        return substrates


class FileSlab(Structure):
    def __init__(self, filename):
        self.filename = filename
        try:
            file_slab = read(self.filename, format="vasp")
        except:
            file_slab = read(self.filename, format="cif")
        self.gratoms = Gratoms(
            positions=file_slab.positions,
            numbers=file_slab.get_atomic_numbers(),
            magmoms=file_slab.get_initial_magnetic_moments(),
            cell=file_slab.cell,
            pbc=[True, True, False],
        )

    def get_dis_inter(self):
        return [2.5, 5]

    def is_dope(self):
        return False

    def find_surface_atoms(
        self, atomstype, tol_zdiff=0.7, tol_zangle_min=0
    ):  # (Last modified: 20230416, zjwang)
        """
        Generate the list of surface atoms (top surface).
        Parameters
        ----------
        atomstype: str
            'top' means find the top surface of slab
            'bottom' means find the bottom surface of slab
        tol_zdiff: number
            If the z_coord of an atom is higher than zmax-tol_zdiff, this atom is recognized as a "surface atom".
        Returns
        ----------
        index_topsurf: list
            List of the indices of the surface atoms.
        """
        if atomstype == "top":
            coords = self.gratoms.positions
        elif atomstype == "bottom":
            coords = -self.gratoms.positions
        else:
            raise ValueError("Only support top and bottom for finding surface atoms")
        index_topsurf = []
        zmax = np.max(coords[:, 2])
        for i, icoord in enumerate(coords):
            if icoord[2] > zmax - tol_zdiff:
                index_topsurf.append(i)
        return index_topsurf

    def construct(self):
        topsurf_atoms = self.find_surface_atoms("top")
        bottomsurf_atoms = self.find_surface_atoms("bottom")
        self.gratoms.set_surface_atoms(top=topsurf_atoms, bottom=bottomsurf_atoms)
        c_fix = FixAtoms(
            indices=[
                atom.index
                for atom in self.gratoms
                if (not atom.index in topsurf_atoms) and (not atom.index in bottomsurf_atoms)
            ]
        )
        self.gratoms.set_constraint(c_fix)
        return [self.gratoms]

    def out_print(self) -> str:
        return "%s substrate" % self.filename

    def out_file_name(self) -> str:
        return self.filename


    @classmethod
    def init_all_slab(cls, input_list: list):
        assert isinstance(input_list, list), "Substrates reading from file should be list!"
        substrates = []
        for filename in input_list:
            assert isinstance(filename, str), "Filename should be str!"
            substrates.append(cls(filename))
        return substrates


def substrate_from_input(init_dict):
    substrates = []
    substrates = substrates + FileSlab.init_all_slab(init_dict["file"])
    substrates = substrates + Slab.init_all_slab(init_dict["struct"])
    return substrates
