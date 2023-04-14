from HTMACat.model.Construct_adsorption_yaml import Construct_adsorption_yaml


def test_dope_type():
    Construct_adsorption_yaml('config_dope_type.yaml')


def test_ads_type():
    Construct_adsorption_yaml('config_ads_type.yaml')

def test_b_dope():
    Construct_adsorption_yaml('config_b1.yaml')

def test_sml():
    Construct_adsorption_yaml('config_sml.yaml')

#test_dope_type()
test_b_dope()