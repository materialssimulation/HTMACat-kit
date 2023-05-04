from HTMACat.model.Substrate import Bulk,Slab,File_Substrate
from HTMACat.catkit.gen.surface import SlabGenerator
from HTMACat.catkit.gen.adsorption import AdsorptionSites
from HTMACat.catkit.gratoms import Gratoms
import pytest
from ase.build import bulk
import numpy as np
import math

@pytest.fixture
def slab1():
     slab1 = Slab()
     return slab1

@pytest.fixture
def slab2():
     bulk = Bulk(natom_dop='1L')
     slab2 = Slab(in_bulk=bulk)
     return slab2

@pytest.fixture
def slab3():
     bulk = Bulk(natom_dop='2')
     slab3 = Slab(in_bulk=bulk,facet='111')
     return slab3

def test_get_miller_index(slab1):
     assert slab1.get_miller_index() == (1, 0, 0)

def test_is_dope(slab1, slab2, slab3):
     assert slab1.is_dope() == False
     assert slab2.is_dope() == False
     assert slab3.is_dope() == True

def test_out_file_name(slab1, slab2):
     assert slab1.out_file_name() =='Pt_100'
     assert slab2.out_file_name() =='Pt_Cu_100_1L'

def test_out_print(slab1, slab2):
     assert slab1.out_print() == 'Pt (100) substrate'
     assert slab2.out_print() == 'Cu 1L doped Pt (100) substrate'

def test_get_dis_inter(slab1,slab3):
     dis_min1 = 3.96
     dis_max1 = 3.96 * 2
     assert np.allclose(slab1.get_dis_inter(), [dis_min1,dis_max1])
     dis_min3 = 3.96 * math.sqrt(2) / 2 * math.sqrt(3) / 2
     dis_max3 = 3.96 * math.sqrt(2) / 2 * math.sqrt(3)
     assert np.allclose(slab3.get_dis_inter(), [dis_min3,dis_max3])


@pytest.fixture
def construct_slab1(slab1):
     slab = []
     miller_index = (1, 0, 0)
     mbulk = slab1.bulk.construct()
     super_cell = [3, 3, 1]
     gen = SlabGenerator(
          mbulk,
          miller_index=miller_index,
          layers=4,
          fixed=2,
          layer_type='trim',
          vacuum=8,
          standardize_bulk=True
     )
     terminations = gen.get_unique_terminations()
     for i,t in enumerate(terminations):
          slab += [gen.get_slab(iterm=i) * super_cell]
     return slab

def test_Construct_slab(construct_slab1, slab1):
     slab1 = slab1.construct()
     for i in range(len(slab1)):
          assert np.allclose(slab1[i].symbols.numbers, construct_slab1[i].symbols.numbers)
          assert np.allclose(slab1[i].positions, construct_slab1[i].positions)


@pytest.fixture
def construct_slab2(slab2):
     slabs = slab2.Construct_slab()
     Ele_dop = 'Cu'
     slabs_dop = []
     for i, slab in enumerate(slabs):
          surface_atoms = slab.get_surface_atoms()
          slb = slab.copy()
          for j, surf_atom in enumerate(surface_atoms):
               slb[surf_atom].symbol = Ele_dop
          slabs_dop += [slb]
     return slabs_dop

def test_Construct_1stLayer_slab(construct_slab2, slab2):
     slab2 = slab2.construct()
     for i in range(len(slab2)):
          assert np.allclose(slab2[i].symbols.numbers, construct_slab2[i].symbols.numbers)
          assert np.allclose(slab2[i].positions, construct_slab2[i].positions)

@pytest.fixture
def construct_slab3(slab3):
     slabs = slab3.Construct_slab()
     slabs_dop = []
     Natom = 2
     for i, slab in enumerate(slabs):
          site = AdsorptionSites(slab)
          site_typ = site.get_connectivity()
          topo = site.get_topology()
          if Natom in site_typ:
               j = np.argwhere(site_typ == Natom)[0][0]
               atom_number_dop = topo[j]
          elif Natom < len(topo[-1]):
               atom_number_dop = topo[-1]
          slb = slab3.dope_slab(slab, atom_number_dop[0:Natom])
          slabs_dop += [slb]
     return slabs_dop

def test_construct_slab3(construct_slab3, slab3):
     slab3 = slab3.construct()
     for i in range(len(slab3)):
          assert np.allclose(slab3[i].symbols.numbers, construct_slab3[i].symbols.numbers)
          assert np.allclose(slab3[i].positions, construct_slab3[i].positions)

def test_file_substrate():
    symbols_numbers = 0
    positions = 0
    file_slab = File_Substrate('./Pt.vasp')
    slab = file_slab.construct()
    assert np.allclose(slab.symbols.numbers, symbols_numbers)
    assert np.allclose(slab.positions, positions)