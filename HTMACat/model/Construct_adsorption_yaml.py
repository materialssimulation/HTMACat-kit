#!/data/jqyang/miniconda3/bin/python
from HTMACat.IO import read, write


def Construct_adsorption_yaml(filename):
    adsorptions = read(filename)
    for ads in adsorptions:
        write(ads)
