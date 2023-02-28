#!/bin/bash
pip install ase numpy
mkdir tmp
cd tmp
git clone https://github.com/materialsproject/fireworks.git
cd fireworks
rm README.md
touch README.md
python setup.py install
cd ../../
rm -rf tmp
pip install git+https://github.com/SUNCAT-Center/CatKit.git
pip install ./requires_wheel/HTMACat-1.0.0-py3-none-any.whl