# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 11:47:30 2023

@author: YuxiaoLan
"""
import yaml
from HTMACat.model.Substrate import substrate_from_input
from HTMACat.model.Ads import ads_from_input
from ase.io.vasp import write_vasp
from HTMACat.model.Structure import Structure

def Input(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        result = yaml.load(f.read(), Loader=yaml.FullLoader)

    # A substrate is one facet with one dop element with on dop_type
    struct_Info = result['StrucInfo']
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
    substrates = []
    for i in range(len(dope_init_list)):
        for j in range(len(surface_init_list)):
            init_dict = {**struct_init_dict, **dope_init_list[i], **surface_init_list[j]}
            substrates.append(substrate_from_input(init_dict))

    ## A adsorption is
    ads_model = result['Model']
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

    return substrates, ads


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
