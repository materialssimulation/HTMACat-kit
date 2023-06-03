**Installation Guide**
======================

.. _environment:

**Environment**
---------------

| Python >= 3.6, <=3.9
| ase >= 3.19
| numpy >= 1.20.0
| networkx >=2.1
| scipy >= 0.1
| spglib >= 0.1
| ruamel.yaml > 0.15
| rdkit >= 2022.3.3
| typer[all] >= 0.6


If you haven't installed **Python3.x** yet, [download](https://www.python.org) the specified version of package and install it.

Recommended for use in virtual environments created by anaconda.

.. _installation:

**Installation**
----------------

Installation from releases

.. code-block:: console

    pip install -U HTMACat

or install from source code

.. code-block:: console

    # master branch
    pip install -U git+https://github.com/stanfordbshan/HTMACat-kit.git@master
    # or dev branch
    pip install -U git+https://github.com/stanfordbshan/HTMACat-kit.git@dev
