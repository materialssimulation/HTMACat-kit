from HTMACat.model.Substrate import Bulk
import pytest
from ase.build import bulk
import numpy as np


def test_bulk():
    a = Bulk(natom_dop=1)
    assert a.lattice_constant["a"] == 3.96
    assert a.super_cell == [3, 3, 1]
    assert a.natom_dop == "1"
    b = Bulk(natom_dop=["b"], lattice_constant=3)
    assert b.super_cell == [2, 2, 1]
    assert b.lattice_constant["a"] == 3


def test_bulk_set_lattice_constant():
    a = Bulk(lattice_constant=3)
    assert a.lattice_constant["a"] == 3
    b = Bulk(lattice_constant=[3.96, 4.96])
    assert b.lattice_constant["a"] == 3.96
    assert b.lattice_constant["c"] == 4.96


@pytest.fixture
def bulk_a():
    mbulk = bulk("Pt", "hcp", a=3.96, covera=4.96 / 3.96, cubic=False)
    return mbulk


@pytest.fixture
def bulk_b():
    mbulk = bulk("Pt", "fcc", a=3.96, cubic=True)
    return mbulk


def test_construct(bulk_a, bulk_b):
    a = Bulk(lattice_type="hcp", lattice_constant=[3.96, 4.96], main_element="Pt")
    a = a.construct()
    assert np.allclose(a.arrays["numbers"], bulk_a.arrays["numbers"])
    assert np.allclose(a.arrays["positions"], bulk_a.arrays["positions"])
    b = Bulk(lattice_type="fcc", lattice_constant=3.96, main_element="Pt")
    b = b.construct()
    assert np.allclose(b.arrays["numbers"], bulk_b.arrays["numbers"])
    assert np.allclose(b.arrays["positions"], bulk_b.arrays["positions"])


@pytest.fixture
def dop_bulk(bulk_b):
    N_dop_bulk = "Cu"
    dop_atoms_number = 1
    for i in range(dop_atoms_number):
        bulk_b[i].symbol = N_dop_bulk
    return bulk_b


def test_dop_bulk(dop_bulk):
    a = Bulk(
        lattice_type="fcc", lattice_constant=3.96, ele_dop="Cu", natom_dop="b1", main_element="Pt"
    )
    a = a.construct()
    assert np.allclose(a.arrays["numbers"], dop_bulk.arrays["numbers"])
    assert np.allclose(a.arrays["positions"], dop_bulk.arrays["positions"])
