from HTMACat.catkit.gen.adsorption import Builder
from HTMACat.catkit.gratoms import *
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import math
from HTMACat.Extract_info import *
from HTMACat.model.Substrate import Slab
from HTMACat.model.Ads import Species
from HTMACat.catkit.gen.adsorption import AdsorptionSites
from HTMACat.model.Structure import Structure
import networkx.algorithms.isomorphism as iso
from ase import Atoms
import pytest

def test_out_file_name():
    species = Species(form='[NH2]',sml=True)
    assert species.out_file_name() == 'H2N'


def test_MolToNXGrath():
    pass

@pytest.fixture
def species1():
    species = Species(form='N')
    return species

@pytest.fixture
def species2():
    species = Species(form='[NH2]',sml=True)
    return species

def test_get_molecule(species1, species2):
    ads_molecule1 = species1.get_molecule()
    assert ads_molecule1.number == 7
    assert np.allclose(ads_molecule1.position, [0,0,0])
    ads_molecule2 = species2.get_molecule()
    assert np.allclose(ads_molecule2.numbers, [7, 1, 1])
    print(ads_molecule2.positions)
    assert np.allclose(ads_molecule2.positions,[[-8.70817145e-04,  3.80919534e-01, -0.00000000e+00],
                                                [-8.70678129e-01, -1.90758492e-01, -0.00000000e+00],
                                                [ 8.71548946e-01, -1.90161042e-01,  0.00000000e+00]])