# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 11:47:30 2023

@author: YuxiaoLan
"""
import shutil
from ruamel.yaml import YAML
from HTMACat.model.Substrate import substrate_from_input
from HTMACat.model.Ads import ads_from_input
from HTMACat.model.Species import species_from_input
from ase.io.vasp import write_vasp
from HTMACat.model.Structure import Structure
from pathlib import Path
from rich import print


def Input(filename):
    yaml = YAML(typ='safe')
    with open(filename, 'r', encoding='utf-8') as f:
        result: dict = yaml.load(f.read())

    substrates = substrate_part(result)
    species_dict = species_part(result)
    ads = adsorption_part(result)
    return substrates, species_dict, ads


def substrate_part(result):
    # A substrate is one facet with one dop element with on dop_type
    struct_Info = result['StrucInfo']
    file_default = []
    struct_default = {'element': 'Pt', 'lattype': 'fcc', 'latcont': 3.92, 'facet': ['111'], 'dope': {}, 'supercell':[3,3,1]}
    if 'file' in struct_Info:
        if isinstance(struct_Info['file'], list):
            file_default = struct_Info['file']
        else:
            file_default = [struct_Info['file']]
    if 'struct' in struct_Info:
        struct_default.update(struct_Info['struct'])
    else:
        struct_default = {}
    substrate_init_dict = {'file': file_default, 'struct': struct_default}
    substrates = substrate_from_input(substrate_init_dict)
    return substrates


def species_part(result: dict):
    if 'Species' not in result.keys():
        return {}
    species_info = result['Species']
    species_dict = species_from_input(species_info)
    return species_dict


def adsorption_part(result):
    # Adsorption Part
    substrates = substrate_part(result)
    species_dict = species_part(result)
    ads_model = result['Model']
    ads_default_dict = {'ads': [], 'coads': []}
    for key in ['ads', 'coads']:
        if key not in ads_model:
            ads_model[key] = ads_default_dict[key]
    ads = ads_from_input(ads_model, substrates, species_dict)
    return ads


def get_templator():
    templator = \
        """#Templator[YAML-format] for HTMACat-kit input
StrucInfo:
# The followting two keywords should be chosen only one
    file: POSCARfile #read substrate from POSCAR file
    struct:
        element: Au #bulk phase element
        lattype: fcc #lattice type
        latcont: 4.16 #lattice parameter
        facet: ['111','100'] #crystal plane
        dope: #dopeing part
            Cu: [3] #dopeing element and dope type
            Ag: [1,2,3,'b1','1L']

Model:
    SML: False #whether the species is in Smile style
    ads: #single adsorption for one species
        - ['NH3',1] #species name and adsorption sites
        - ['NO', 2]
    coads: #co-adsorption for one species
        - ['NH3','O',1,1] #species name and adsorption sites"""
    return templator


def print_templator():
    templator = get_templator()
    print(templator)


def out_templator_file():
    templator = get_templator()
    workdir = Path('./config.yaml')
    if workdir.is_file():
        while True:
            choose = input('Warning: Replace config.yaml file in your folder? [y/n(default)]:')
            if choose == 'y':
                break
            elif choose == 'n' or choose == '':
                return
            else:
                continue

    with open('config.yaml', 'w') as f:
        f.write(templator)


def out_vasp(struct_class: Structure):
    file_name = struct_class.out_file_name()
    print_str = struct_class.out_print()
    structures = struct_class.construct()
    if structures:
        for i, structure in enumerate(structures):
            write_vasp("%s_%s.vasp" % (file_name, i), structure, direct=True, sort=[''], vasp5=True)
        print(f"[i][u]{print_str}[/u][/i] are construct with [i][u]{i + 1}[/u][/i] configurations!")
    else:
        print(f"[i][u]{print_str}[/u][/i] are construct with !No! configurations!")
