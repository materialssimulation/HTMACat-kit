import shutil
from ruamel.yaml import YAML
from HTMACat.model.Substrate import substrate_from_input
from HTMACat.model.Ads import ads_from_input
from HTMACat.model.Species import species_from_input
from ase.io.vasp import write_vasp
from HTMACat.model.Structure import Structure
from HTMACat.IO import substrate_part, species_part, adsorption_part, Input
from pathlib import Path
from rich import print
import pytest


def test_input():
    f = "./test/test_IO.yaml"
    substrates, species_dict, ads = Input(f)
    assert len(ads) == 5


@pytest.fixture()
def result():
    result = {
        "StrucInfo": {
            "struct": {
                "element": "Pt",
                "lattype": "fcc",
                "latcont": 3.92,
                "facet": ["111"],
                "dope": {"Cu": [0]},
            }
        },
        "Species": {"file": {"NH3+": "NH3+.xyz"}, "sml": {"O2": "O=O"}},
        "Model": {
            "ads": [["NH3+", 1], [{"f": "NH3-.xyz"}, 1], ["NO", 2], ["O2", 2], [{"s": "O=O"}, 1]],
            "coads": [],
        },
    }
    return result


def test_substrate_part(result):
    substrates = substrate_part(result)
    assert substrates[0].facet == "111"
    assert substrates[0].bulk.lattice_constant["a"] == 3.96
    assert substrates[0].bulk.main_element == "Pt"
    assert substrates[0].bulk.lattice_type == "fcc"
    assert substrates[0].bulk.ele_dop == "Cu"
    assert substrates[0].bulk.natom_dop == "0"


def test_species_part(result):
    species_dict = species_part(result)
    assert species_dict["NH3+"].form == "NH3+.xyz"
    # assert species_dict['NH3+']._ABS_Speices__type =='file'
    assert species_dict["O2"].form == "O=O"
    # assert species_dict['O2']._ABS_Speices__type =='sml'


def test_adsorption_part(result):
    ads = adsorption_part(result)
    assert len(ads) == 5


# def
# result = {
#     'StrucInfo': {
#         'struct': {
#             'element': 'Pt',
#             'lattype': 'fcc',
#             'latcont': 3.92,
#             'facet': ['111'],
#             'dope': {'Cu': [0]}
#         }
#     },
#     'Species': {'file': {'NH3+': 'NH3+.xyz'}, 'sml': {'O2': 'O=O'}},
#     'Model': {
#         'ads': [
#             ['NH3+', 1],
#             [{'f': 'NH3-.xyz'}, 1],
#             ['NO', 2],
#             ['O2', 2],
#             [{'s': 'O=O'}, 1]
#         ],
#         'coads': []
#     }
# }
# test = adsorption_part(result)

# f = './test/test_IO.yaml'
# substrates, species_dict, ads = Input(f)
# print()
