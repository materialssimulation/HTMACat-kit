**Getting started**
===================

**1.High-throughput adsorption modeling**
------------------------------------------------

(1) Prepare the input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``config.yaml`` should be prepared in your working folder.

``Config.yaml`` file Contains two parts:

* ``StuctInfo`` parts
* ``Species`` parts (optional)
* ``Model`` parts

The format of the whole `config.yaml` is as follows:

.. code-block:: console

    StrucInfo:
      file: POSCARfile
      struct:
        element: Au
        lattice_type: fcc
        latttice_constant: 4.16
        facet: ['111','100'] #This is comment
        supercell: [3,3]     #Default: [3,3], can be undefined
        layers: 4            #Default: 4, can be undefined
        dope:
          Cu: [3]

    # This part alias for Species
    Species:
      file:
        'NH3+': NH3+.xyz #struct_file
        'NH3-': NH3-.xyz #struct_file
      sml:
        'O2': O=O #smiles

    Model:
      ads:
        - ['NH3+',1]
        - ['NO', 2]
        - [s: 'C(=O)O',1]
        - [f: 'NH3+',1]
        - [s: 'C[Al](C)C',1Al]
      coads:
        - ['NH3','O',1,1]

``StuctInfo`` parts contains the information about substrate,
the keywords should be chosen as following two types: ``file`` or ``struct``

* ``file``: Path , this keywords will read substrate from POSCAR file
* ``struct``: this keywords will generate substrate according to following parameter:
    * ``element``: bulk phase element
    * ``lattice_type``: lattice type
    * ``latticd_constant``: lattice parameter
    * ``facet``: crystal plane, should be a *list*, start with ``[``, separate with ``,`` and end with ``]``,
      facet index should be a str start with ``'`` end with ``'``, like ``'100'``, ``'111'``
    * ``supercell``: supercell of the substrate in the xy-plane
    * ``layers``: layers of the substrate in z axis
    * ``dope``: dope of the substrate, the formate is ``dope element : [dope type1, dope type2]``
      before ``:`` is the doped element, after the ``:`` is the dope type, the dope type can be
      chosen as follows: ``0`` corresponds to no doping, ``1``, ``2`` and ``3`` correspond to surface
      layers doped with ``1``, ``2`` and ``3`` atoms, respectively. ``1L`` and ``b1`` represent surface layer
      substitution and bulk equivalent proportional substitution.

``Species`` part is an optional part, which used to alias the species name, including:

* ``sml``: using *Smiles* toe represent the adsorption species
* ``file``: adsorption species reading from structure file, such as *xyz* format and *vasp* format

``Model`` part contains the adsorption modeling parameters, including:

* ``ads``: using ``- [ adsorbate formular , adsorption sites type]`` to represent one adsorption status, where
  adsorption formular should be a str start with ``'`` end with ``'``. Different adsorption status should start with
  new line, adsorption sites type can be chosen from ``1`` or ``2``.
* ``coads``: using ``- [ads1, ads2, ads1 sites, ads2 sites]``, ``ads1`` and ``ads2`` is two adsorbate species formular,
  ``ads1 sites`` and ``ads2 sites`` is adsorption sites type, can be chosen from ``1`` and ``2``.
* Each *adsorbate species formular* has three representation:
    * *str* format such as ``'ads1'``, this format will automatically generate species according to
      chemical formular
    * *smile* format such as ``s: 'ads1'`` or ``sml: 'ads1'``, this format will generate spcecies according to
      the *smile* formular
    * *file* format such as ``f: 'ads1'`` or ``file: 'ads1'``, this format will read structure from
      file 'ads1', now only 'xyz' and 'vasp' structure file format is supported
* New in version 1.0.5: Adding the symbol of the element contained in the adsorbate followed by the type of adsorption site 
  generates a configuration in which all atoms corresponding to that element are adsorbed atoms, and automatically adjusts the 
  orientation so that the adsorbed atom is at the bottom. Such as ``- [s: 'C[Al](C)C',1Al]`` will generatte configurations with 
  Al atoms as adsorbed stoms

**To avoid ambiguity, it is recommended that SMILES be used when declaring complex species.** Users can modify the
corresponding parameters to achieve customized modeling according to their research needs.

(2) Run script
~~~~~~~~~~~~~~

Running command ``htmat ads`` can automate the enumeration and construction of all possible configurations, and ultimately
output structure files in the VASP format, like Au_Cu_111_1_NH_0.vasp,
Furthermore, for the adsorption models on doped surfaces, we use the ``get_binding_adatom`` function to determine the
types of surface atoms that the adsorbate binds to in the preliminary configurations. We only select configurations
where the adsorbate is bound to a surface atom that contains a dopant atom, in order to reduce the number of final
output structures.

(3) Completion of the input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the command ``htmat complete``, it will add the default parameters according to the current config.yaml and output 
the complete_config.json file.


**2.Automated construction of reaction transition state calculation process**
-----------------------------------------------------------------------------

**3.Automated extraction of calculation results**
-------------------------------------------------

**4.Automated extraction of descriptors**
-----------------------------------------
