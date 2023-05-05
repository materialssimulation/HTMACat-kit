# HTMACat

# 1. Introduction

A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes, it provides key tools
for high-throughput design and screening of catalytic materials. The software mainly includes functional modules such as
surface structure analysis and information extraction, catalytic surface and various adsorption model construction,
automatic construction of primitive reaction processes, automatic extraction of computational data and automatic
extraction and construction of descriptors.The software can perform the following computational workflows: adsorption
energy calculation and analysis workflow, primitive reaction calculation and analysis workflow, high throughput
calculation and automated analysis of adsorption energy and reaction potential of catalytic primitive reaction
processes, etc.

# 2. Installation Guide

## 2.1 Environment

```requirements.txt
Python >= 3.6, <=3.9
ase >= 3.19
numpy >= 1.20.0
networkx >=2.1
scipy >= 0.1
spglib >= 0.1
ruamel.yaml > 0.15
rdkit >= 2022.3.3
typer[all] >= 0.6
```
## 2.2  Installation
install wheel file from releases

```shell
pip install -U HTMACat

```
or install from source code
```shell
# master branch
pip install -U git+https://github.com/stanfordbshan/HTMACat-kit.git@master
# or dev branch
pip install -U git+https://github.com/stanfordbshan/HTMACat-kit.git@dev
```

# 3. Getting started

## 3.1 High-throughput adsorption modeling

### 1) Prepare the input files

`config.yaml` should be prepared in your working folder.
`Config.yaml` file Contains three parts: 
- `StuctInfo` part 
- `Species` part (optional)
- `Model` part

The format of the whole `config.yaml` is as follows:
   ```yaml
   StrucInfo:
     file: POSCARfile
     struct:
        element: Au
        lattype: fcc
        latcont: 4.16
        facet: ['111','100'] #This is comment
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
    coads: 
      - ['NH3','O',1,1]
   ```

`StuctInfo` parts contains the information about substrate, 
the keywords should be chosen as following two types: `file` or `struct`
- `file`: Path , this keywords will read substrate from POSCAR file 
- `struct`: this keywords will generate substrate according to following parameter:
  - `element`: bulk phase element
  - `lattype`: lattice type
  - `latcont`: lattice parameter
  - `facet`: crystal plane, hould be a *list*, start with `[`, separate with `,` and end with `]`, facet index
    should
    be a str start with `'` end with `'`, like `'100'`,`'111'`
  - `dope`: dope of the substrate, the formate is `dope element : [dope type1, dope type2]`
    before `:` is the doped element, after the `:` is the dope type, the dope type can be chosen as follows：
    `0` corresponds to no doping, `1`, `2` and `3`correspond to surface layers doped with `1`, `2` and `3` atoms,
    respectively.
    `1L` and `b1` represent surface layer substitution and bulk equivalent proportional substitution.

`Species` part is an optional part, which used to alias the species name, including:
- `sml`: using *Smiles* toe represent the adsorption species
- `file`: adsorption species reading from structure file, such as *xyz* format and 
          *vasp* format

`Model` part contains the adsorption modeling parameters, including:
- `ads`: using `- [ adsorbate formular , adsorption sites type]` to represent one adsorption status, where
  adsorption formular should be a str start with `'` end with `'`. Different adsorption status should start with
  new line, adsorption sites type can be chosen from `1` or `2`.
- `coads`: using `- [ads1, ads2, ads1 sites, ads2 sites]`, *ads1* and *ads2* is two adsorbate species formular,
  *ads1 sites* and *ads2 sites* is adsorption sites type, can be chosen from `1` and `2`.
- Each *adsorbate species formular* has three representation:
  - *str* format such as `'ads1'`, this format will automatically generate species according to
    chemical formular
  - *smile* format such as `s: 'ads1'` or `sml: 'ads1'`, this format will generate spcecies according to
    the *smile* formular
  - *file* format such as `f: 'ads1'` or `file: 'ads1'`, this format will read structure from
    file 'ads1', now only 'xyz' and 'vasp' structure file format is supported  

**To avoid ambiguity, it is recommended that SMILES be used when declaring complex species.** Users can modify the
corresponding parameters to achieve customized modeling according to their research needs.

### 2) run script

Running command `htmat ads` can automate the enumeration and construction of all possible configurations, and ultimately
output structure files in the VASP format, like Au_Cu_111_1_NH_0.vasp, 
Furthermore, for the adsorption models on doped surfaces, we use the `get_binding_adatom` function to determine the
types of surface atoms that the adsorbate binds to in the preliminary configurations. We only select configurations
where the adsorbate is bound to a surface atom that contains a dopant atom, in order to reduce the number of final
output structures.


## 3.2 Automated construction of reaction transition state calculation process
## 3.3 Automated extraction of calculation results
## 3.4 Automated extraction of descriptors

‍

# Author

* Jiaqiang Yang email:[jqyang_hust@hust.edu.cn](mailto:jqyang_hust@hust.edu.cn)
* Feifeng Wu email:[wufeifeng_hust@163.com](wufeifeng_hust@163.com)
* Bin Shan email:[bshan@mail.hust.edu.cn](bshan@mail.hust.edu.cn)

‍

# Links

* Materials Design and Nano-Manufacturing Center@HUST:http://www.materialssimulation.com/
* Pypi homepage for HTMACat project:https://pypi.org/project/HTMACat/
* HTMACat-kit’s documentation:https://stanfordbshan.github.io/HTMACat-kit/