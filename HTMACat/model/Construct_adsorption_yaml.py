#!/data/jqyang/miniconda3/bin/python
# -*- coding: UTF-8 -*-
from HTMACat.IO import Input, out_vasp


def Construct_adsorption_yaml(filename):
    substrates, species_dict, adsorptions = Input(filename)
    for ads in adsorptions:
        out_vasp(ads)
