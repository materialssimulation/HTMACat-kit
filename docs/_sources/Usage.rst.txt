**Usage**
===========

**1.High-throughput single adsorption modeling**
-----------

(1) Prepare the input files: *StrucInfo* and *Model*
~~~~~~~~~~~~~~

"StrucInfo" file Contains the carrier of lattice type and the crystal parameters, the format is as follows:

.. code-block:: console

    Pt fcc 3.97 111 # Element Lattice_type Lattice_parameter Crystal_plane

The four parameters above are the bulk phase element, lattice type, lattice parameter and crystal plane. In addition, by adding more rows, more systems can be constructed continuously.  

"Model" file contains the surface doping and adsorption modeling parameters, the format is as follows:

.. code-block:: console

    Cu # Doping elements
    b1 # Doping type
    SML C(=O)O # Adsorption species (SMILES of formic acid)
    1 2 # Adsorption type

Where the first line of the file is the doping element and the second line is the doping type, 0 corresponds to no doping, 1, 2 and 3 correspond to surface layers doped with 1, 2 and 3 atoms, respectively, and 1L and b1 represent surface layer substitution and bulk equivalent proportional substitution. The third and fourth lines respectively represent the adsorbate species and the adsorption types. 
Please note that the adsorbates should be declared by their SMILES expressions if the third line starts with 'SML'. 
If no 'SML' at the beginning of this line, the species are represented in their chemical formula. 
To avoid ambiguity, it is recommended that SMILES be used when declaring complex species.
Users can modify the corresponding parameters to achieve customized modeling according to their research needs.

(2) Run script
~~~~~~~~~~~~

Running command "ads" can automate the enumeration and construction of all possible configurations, and ultimately output structure files in the VASP format, like Pt_Cu_111_1_NH_0.vasp.
Furthermore, for the adsorption models on doped surfaces, we use the get_binding_adatom function to determine the types of surface atoms that the adsorbate binds to in the preliminary configurations. We only select configurations where the adsorbate is  bound to a surface atom that contains a dopant atom, in order to reduce the number of final output structures.

**2.High-throughput co-adsorption modeling**
-------------

(1) Prepare the input files:  *StrucInfo* and *Model*
~~~~~~~~~~~~~~

The two input files for the co-adsorption model construction module are still StrucInfo and Model, and StrucInfo has the same file format as the single adsorption model construction module. However, the Model files for the two differ. 
The format of the Model file for the co-adsorption modeling module is as follows:

.. code-block:: console

    Cu # Doping elements
    3 # Doping type
    SML N [O];[NH2] [OH];[N] [N];[N] [O];[N] [N]=O # Co-adsorption species
    1 1;1 1;1 1;1 1;1 1 # Adsorption type

The first line represents the dopant element, the second line represents the dopant type, the third and fourth lines respectively represent the co-adsorbate species (also in SMILES format) and the co-adsorption type. Different co-adsorption combinations are separated by ";".

(2) Run script
~~~~~~~~~~~~~~~~
Running command "coads" can automate the enumeration and construction of all possible configurations, and ultimately output structure files in the VASP format, like Pt_Cu_111_2_N_H_0.vasp.

**3.Automated construction of reaction transition state calculation process**
-------------------

**4.Automated extraction of calculation results**
-------------------

**5.Automated extraction of descriptors**
-------------------