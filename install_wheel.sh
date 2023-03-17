#!/bin/bash
pip install ase numpy==1.23.5 scikit-learn
pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/FireWorks-2.0.3-py3-none-any.whl
pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/CatKit-0.5.4-py3-none-any.whl
pip install --upgrade $(curl -s https://api.github.com/repos/stanfordbshan/HTMACat-kit/releases/latest | grep "browser_download_url.*\.whl" | cut -d : -f 2,3 | tr -d \")