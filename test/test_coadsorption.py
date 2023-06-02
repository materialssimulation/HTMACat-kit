from HTMACat.catkit.gen.adsorption import Builder
from HTMACat.catkit.gratoms import *
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import math
from HTMACat.Extract_info import *
from HTMACat.model.Substrate import Slab
from HTMACat.model.Ads import Adsorption, Coadsorption
from HTMACat.model.Species import Sml_Species
from HTMACat.catkit.gen.adsorption import AdsorptionSites
from HTMACat.model.Structure import Structure
import networkx.algorithms.isomorphism as iso
from ase import Atoms
import pytest
import numpy as np


@pytest.fixture
def species():
    species1 = Sml_Species(form="N")
    species2 = Sml_Species(form="[O]")
    return [species1, species2]


@pytest.fixture
def coads11(species):
    sites = ["1", "1"]
    ads = Coadsorption(species=species, sites=sites)
    return ads


def test_out_file_name(coads11):
    assert coads11.out_file_name() == "Pt_100_H3N_O"


def test_construct_coadsorption_11(coads11):
    slab_ad_final = coads11.Construct_coadsorption_11()
    assert len(slab_ad_final) == 40
    assert np.allclose(
        slab_ad_final[0].numbers,
        [
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            78,
            7,
            1,
            1,
            1,
            8,
        ],
    )
    print(slab_ad_final[0].positions)
    # assert np.allclose(slab_ad_final[0].positions, [[ 1.40007143e+00,  1.40007143e+00,  8.00000000e+00],
    #                                                 [-1.03582051e-16, -1.47387144e-16,  9.98000000e+00],
    #                                                 [ 1.40007143e+00,  1.40007143e+00,  1.19600000e+01],
    #                                                 [ 0.00000000e+00,  0.00000000e+00,  1.39400000e+01],
    #                                                 [ 1.40007143e+00,  4.20021428e+00,  8.00000000e+00],
    #                                                 [-1.03582051e-16,  2.80014285e+00,  9.98000000e+00],
    #                                                 [ 1.40007143e+00,  4.20021428e+00,  1.19600000e+01],
    #                                                 [ 0.00000000e+00,  2.80014285e+00,  1.39400000e+01],
    #                                                 [ 1.40007143e+00,  7.00035713e+00,  8.00000000e+00],
    #                                                 [-1.03582051e-16,  5.60028571e+00,  9.98000000e+00],
    #                                                 [ 1.40007143e+00,  7.00035713e+00,  1.19600000e+01],
    #                                                 [ 0.00000000e+00,  5.60028571e+00,  1.39400000e+01],
    #                                                 [ 4.20021428e+00,  1.40007143e+00,  8.00000000e+00],
    #                                                 [ 2.80014285e+00, -1.47387144e-16,  9.98000000e+00],
    #                                                 [ 4.20021428e+00,  1.40007143e+00,  1.19600000e+01],
    #                                                 [ 2.80014285e+00,  0.00000000e+00,  1.39400000e+01],
    #                                                 [ 4.20021428e+00,  4.20021428e+00,  8.00000000e+00],
    #                                                 [ 2.80014285e+00,  2.80014285e+00,  9.98000000e+00],
    #                                                 [ 4.20021428e+00,  4.20021428e+00,  1.19600000e+01],
    #                                                 [ 2.80014285e+00,  2.80014285e+00,  1.39400000e+01],
    #                                                 [ 4.20021428e+00,  7.00035713e+00,  8.00000000e+00],
    #                                                 [ 2.80014285e+00,  5.60028571e+00,  9.98000000e+00],
    #                                                 [ 4.20021428e+00,  7.00035713e+00,  1.19600000e+01],
    #                                                 [ 2.80014285e+00,  5.60028571e+00,  1.39400000e+01],
    #                                                 [ 7.00035713e+00,  1.40007143e+00,  8.00000000e+00],
    #                                                 [ 5.60028571e+00, -1.47387144e-16,  9.98000000e+00],
    #                                                 [ 7.00035713e+00,  1.40007143e+00,  1.19600000e+01],
    #                                                 [ 5.60028571e+00,  0.00000000e+00,  1.39400000e+01],
    #                                                 [ 7.00035713e+00,  4.20021428e+00,  8.00000000e+00],
    #                                                 [ 5.60028571e+00,  2.80014285e+00,  9.98000000e+00],
    #                                                 [ 7.00035713e+00,  4.20021428e+00,  1.19600000e+01],
    #                                                 [ 5.60028571e+00,  2.80014285e+00,  1.39400000e+01],
    #                                                 [ 7.00035713e+00,  7.00035713e+00,  8.00000000e+00],
    #                                                 [ 5.60028571e+00,  5.60028571e+00,  9.98000000e+00],
    #                                                 [ 7.00035713e+00,  7.00035713e+00,  1.19600000e+01],
    #                                                 [ 5.60028571e+00,  5.60028571e+00,  1.39400000e+01],
    #                                                 [-1.57500777e-15, -3.67705866e-15,  1.60100000e+01],
    #                                                 [ 9.79326751e-01, -2.15252258e-01,  1.57214005e+01],
    #                                                 [-6.78657353e-01, -7.24867079e-01,  1.57201624e+01],
    #                                                 [-2.91398819e-01,  9.60519651e-01,  1.57219977e+01],
    #                                                 [ 7.00035713e+00,  2.58653660e-15,  1.53960907e+01]])


# @pytest.fixture
# def coads12(species):
#     sites=['1','2']
#     coads = Coadsorption(species=species, sites=sites)
#     return coads

# def test_construct_coadsorption_12(coads12):
#     slab_ad_final = coads12.Construct_coadsorption_12()

# @pytest.fixture
# def coads22(species):
#     sites=['2','2']
#     coads = Coadsorption(species=species, sites=sites)
#     return coads

# def test_construct_coadsorption_22(coads22):
#     slab_ad_final = coads22.Construct_coadsorption_22()
