from HTMACat.api import construct_adsorption

StrucInfo = {
    "struct": {
        "element": "Pt",
        "lattice_type": "fcc",
        "lattice_constant": 4.16,
        "facet": ["100"],
        "dope": {"Cu": ["0"]}
    }
}

Model = {
    "ads": [
        [{'s': "O"}, '1O', {"settings": {"site_coords": [[1, 1, 2]], "direction": "asphericity"}}],
        [{'s': "O"}, '1O', {"settings": {"direction": "asphericity"}}]
    ]
}

construct_adsorption(StrucInfo=StrucInfo, Model=Model, workdir=".")
