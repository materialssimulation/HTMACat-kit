# HTMACat

# Introduction

A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes, it provides key tools for high-throughput design and screening of catalytic materials. The software mainly includes functional modules such as surface structure analysis and information extraction, catalytic surface and various adsorption model construction, automatic construction of primitive reaction processes, automatic extraction of computational data and automatic extraction and construction of descriptors.The software can perform the following computational workflows: adsorption energy calculation and analysis workflow, primitive reaction calculation and analysis workflow, high throughput calculation and automated analysis of adsorption energy and reaction potential of catalytic primitive reaction processes, etc.

# Installation Guide


## 1.  Environment
Python 3.6-3.9; ASE 3.22.1; CatKit 0.5.4; 1.20.0 <= Numpy <= 1.23.5 

If you haven‘t installed **Python3.x** yet, [download](https://www.python.org) the specified version of package and install it.

## 2.  Installation

```shell
pip install ase numpy==1.23.5 scikit-learn
pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/FireWorks-2.0.3-py3-none-any.whl
pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/CatKit-0.5.4-py3-none-any.whl
pip install -U https://github.com/stanfordbshan/HTMACat-kit/releases/download/v1.0.2/HTMACat-1.0.2-py3-none-any.whl

```

# Getting started

1. **High-throughput single adsorption modeling**

    (1) Prepare the input files: ​*StrucInfo ​*and ​*Model*

    "StrucInfo" file Contains the carrier of lattice type and the crystal parameters, the format is as follows:

    ```shell
    Pt fcc 3.97 111 # Element Lattice_type Lattice_parameter Crystal_plane
    ```
    The four parameters above are the bulk phase element, lattice type, lattice parameter and crystal plane. In addition, by adding more rows, more systems can be constructed continuously.  
    "Model" file contains the surface doping and adsorption modeling parameters, the format is as follows:

    ```shell
    Cu # Doping elements
    b1 # Doping type
    SML C(=O)O # Adsorption species (SMILES of formic acid)
    1 2 # Adsorption type
    ```
    Where the first line of the file is the doping element and the second line is the doping type, 0 corresponds to no doping, 1, 2 and 3 correspond to surface layers doped with 1, 2 and 3 atoms, respectively, and 1L and b1 represent surface layer substitution and bulk equivalent proportional substitution. The third and fourth lines respectively represent the adsorbate species and the adsorption types. <b>Please note that the adsorbates should be declared by their SMILES expressions if the third line starts with 'SML'!</b> If no 'SML' at the beginning of this line, the species are represented in their chemical formula. <b>To avoid ambiguity, it is recommended that SMILES be used when declaring complex species.</b> Users can modify the corresponding parameters to achieve customized modeling according to their research needs.

    (2) run script

    Running command "ads" can automate the enumeration and construction of all possible configurations, and ultimately output structure files in the VASP format, like Pt_Cu_111_1_NH_0.vasp  
    Furthermore, for the adsorption models on doped surfaces, we use the get_binding_adatom function to determine the types of surface atoms that the adsorbate binds to in the preliminary configurations. We only select configurations where the adsorbate is bound to a surface atom that contains a dopant atom, in order to reduce the number of final output structures.
2. **High-throughput co-adsorption modeling**

    (1) Prepare the input files: ​*StrucInfo ​*and ​*Model*

    The two input files for the co-adsorption model construction module are still StrucInfo and Model, and StrucInfo has the same file format as the single adsorption model construction module. However, the Model files for the two differ. The format of the Model file for the co-adsorption modeling module is as follows:

    ```shell
    Cu # Doping elements
    3 # Doping type
    SML N [O];[NH2] [OH];[N] [N];[N] [O];[N] [N]=O # Co-adsorption species
    1 1;1 1;1 1;1 1;1 1 # Adsorption type
    ```
    The first line represents the dopant element, the second line represents the dopant type, the third and fourth lines respectively represent the co-adsorbate species (also in SMILES format) and the co-adsorption type. Different co-adsorption combinations are separated by ";".

    (2) run script

    Running command "coads" can automate the enumeration and construction of all possible configurations, and ultimately output structure files in the VASP format, like Pt_Cu_111_2_N_H_0.vasp.
3. **Automated construction of reaction transition state calculation process**
4. **Automated extraction of calculation results**
5. **Automated extraction of descriptors**

‍

# Author

* Jiaqiang Yang email:[jqyang_hust@hust.edu.cn](mailto:jqyang_hust@hust.edu.cn)
* Feifeng Wu email:[wufeifeng_hust@163.com](wufeifeng_hust@163.com)
* Bin Shan email:[bshan@mail.hust.edu.cn](bshan@mail.hust.edu.cn)

‍

# Links

* Materials Design and Nano-Manufacturing Center@HUST:http://www.materialssimulation.com/
