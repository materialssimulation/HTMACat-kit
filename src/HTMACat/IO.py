# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 11:47:30 2023

@author: YuxiaoLan
"""
import shutil

from ruamel.yaml import YAML
from HTMACat.model.Substrate import substrate_from_input, substrate_from_file
from HTMACat.model.Ads import ads_from_input
from ase.io.vasp import write_vasp
from HTMACat.model.Structure import Structure
from pathlib import Path


def Input(filename):
    yaml = YAML(typ='safe')
    with open(filename, 'r', encoding='utf-8') as f:
        result = yaml.load(f.read())
    struct_Info = result['StrucInfo']
    substrates = substrate_part(struct_Info)
    ads_model = result['Model']
    ads = adsorption_part(ads_model, substrates)
    return substrates, ads


def substrate_part(struct_Info):
    # A substrate is one facet with one dop element with on dop_type
    substrates = []
    if 'file' in struct_Info:
        struct_filename = struct_Info['file']
        substrates.append(substrate_from_file(struct_filename))
        return substrates

    struct_filename = None
    struct_Info = struct_Info['struct']
    struct_init_dict = {'element': struct_Info['element'],
                        'lattype': struct_Info['lattype'],
                        'latcont': struct_Info['latcont']}
    dope_init_dict = {'element_dop': [], 'dop_type': []}
    surface_init_dict = {'facet': []}
    dope_system = struct_Info['dope']
    dope_init_list = []
    surface_init_list = []
    for key, value in dope_system.items():
        for i in value:
            dope_init_list.append({'element_dop': key, 'dop_type': i})
    for i in struct_Info['facet']:
        surface_init_list.append({'facet': i})

    ## substrates initialization
    if not struct_filename:
        for i in range(len(dope_init_list)):
            for j in range(len(surface_init_list)):
                init_dict = {**struct_init_dict, **dope_init_list[i], **surface_init_list[j]}
                substrates.append(substrate_from_input(init_dict))
    return substrates


def adsorption_part(ads_model, substrates):
    ## Adsorption Part
    ads_default_dict = {'SML': False, 'ads': [], 'coads': []}
    ads = []

    for key in ['SML', 'ads', 'coads']:
        if key not in ads_model:
            ads_model[key] = ads_default_dict[key]

    for key in ['ads', 'coads']:
        for j in ads_model[key]:
            if not j:
                continue
            for k in substrates:
                ads_init_dict = {'SML': ads_model['SML'], 'type': key, 'value': j, 'substrate': k}
                ads.append(ads_from_input(ads_init_dict))

    return ads


def get_templator():
    templator = \
        """#Templator for HTMACat-kit input
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
        print(print_str + " are construct with %d configurations!" % (i + 1))
    else:
        print(print_str + " are construct with !No! configurations!")
