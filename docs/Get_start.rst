**Get_start**
============

.. _environment:

**Environment**
---------------

Python 3.6-3.9; ASE 3.22.1; CatKit 0.5.4; 1.20.0 <= Numpy <= 1.23.5 

If you haven‘t installed **Python3.x** yet, [download](https://www.python.org) the specified version of package and install it.

Recommended for use in virtual environments created by anaconda.

.. _installation:

**Installation**
----------------

Automatic installation
~~~~~~~~~~~~~~~~~~~~~~~~

download the install script:[install_wheel.bat](https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/install_wheel.bat) or [install_wheel.sh](https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/install_wheel.sh)


For windows

.. code-block:: console
    
    ./install_wheel.bat

For linux

.. code-block:: console
    
    ./install_wheel.sh 

Manual installation
~~~~~~~~~~~~~~~~~~~~~~~~~

install requires package

.. code-block:: console

    pip install ase numpy
    pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/FireWorks-2.0.3-py3-none-any.whl
    pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/CatKit-0.5.4-py3-none-any.whl

from setup.py

.. code-block:: console

    python setup.py install

from wheel

.. code-block:: console

    pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/HTMACat-1.0.0-py3-none-any.whl

