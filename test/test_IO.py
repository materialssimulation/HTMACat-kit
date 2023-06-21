import shutil
from ruamel.yaml import YAML
from HTMACat.model.Substrate import substrate_from_input
from HTMACat.model.Ads import ads_from_input
from HTMACat.model.Species import species_from_input
from ase.io.vasp import write_vasp
from HTMACat.model.Structure import Structure
from HTMACat.IO import yaml2dict, dict2object
from pathlib import Path
from rich import print
import pytest


def test_yaml2dict():
    f = "./test/test_IO.yaml"

def test_dict2object():
    pass

