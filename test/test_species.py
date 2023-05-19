from HTMACat.catkit.gen.adsorption import Builder
from HTMACat.catkit.gratoms import *
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import math
from HTMACat.Extract_info import *
from HTMACat.model.Substrate import Slab
from HTMACat.model.Species import (
    Sim_Species,
    File_Species,
    Sml_Species,
    species_from_input,
    init_from_ads,
)
from HTMACat.catkit.gen.adsorption import AdsorptionSites
from HTMACat.model.Structure import Structure
import networkx.algorithms.isomorphism as iso
from ase import Atoms
import pytest


@pytest.fixture
def sim_species():
    species = Sim_Species("NH3")
    return species


@pytest.fixture
def file_species():
    species = File_Species("./test/NH3+.xyz")
    return species


@pytest.fixture
def sml_species():
    species = Sml_Species("[NH2]")
    return species


def test_get_molecule(sim_species, file_species, sml_species):
    # sim_species
    sim_molecule = sim_species.get_molecule()
    assert np.allclose(sim_molecule.numbers, [7, 1, 1, 1])
    assert np.allclose(
        sim_molecule.positions,
        [
            [0.0, 0.0, 0.116489],
            [0.0, 0.939731, -0.271808],
            [0.813831, -0.469865, -0.271808],
            [-0.813831, -0.469865, -0.271808],
        ],
    )
    # file_species
    file_molecule = file_species.get_molecule()
    assert np.allclose(file_molecule.numbers, [7, 1, 1, 1])
    assert np.allclose(
        file_molecule.positions,
        [
            [0.0, 0.0, 0.116489],
            [0.0, 0.939731, 0.408],
            [0.813831, -0.469865, 0.40808],
            [-0.813831, -0.469865, 0.40808],
        ],
    )
    # sml_species
    sml_molecule = sml_species.get_molecule()
    print(sml_species.form)
    # assert np.allclose(sml_molecule.numbers, [7, 1, 1])
    # assert np.allclose(sml_molecule.positions, [[-0.00569876,  0.40600421, -0.        ],
    #                                             [-0.84031275, -0.20509521, -0.        ],
    #                                             [ 0.84601151, -0.200909  ,  0.        ]])
    assert sml_molecule.get_chemical_formula() == "H2N"


def test_species_from_input():
    init_dict = {"sml": {"O2": "O=O"}, "file": {"NH3+": "NH3+.xyz"}}
    species_dict = species_from_input(init_dict)
    assert species_dict["O2"].form == "O=O"
    assert species_dict["NH3+"].form == "NH3+.xyz"


def test_init_from_ads():
    init_dict = {"sml": {"O2": "O=O"}, "file": {"NH3+": "NH3+.xyz"}}
    species_dict = species_from_input(init_dict)

    init_str = "O"
    species = init_from_ads(init_str, species_dict)
    assert species.form == "O"

    init_dict = {"s": "N"}
    species = init_from_ads(init_dict, species_dict)
    assert species.form == "N"

    init_str_dict = {"f": "NH3+.xyz"}
    species = init_from_ads(init_str_dict, species_dict)
    assert species.form == "NH3+.xyz"


species = Sml_Species("[NH2]")
ads_molecule = species.get_molecule()
print(ads_molecule)
# init_dict = {'sml':{'O2': 'O=O', 'CO': '[C]=O'}, 'file':{'NH3+': 'NH3+.xyz'}}
# species_dict = species_from_input(init_dict)
# print(species_dict)
# print(species_dict['O2'].form)
# print(species_dict['NH3+'].form)
# species = File_Species('./test/NH3+.xyz')
# file_molecule =  species.get_molecule()
# print()
