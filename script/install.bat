rd /S /Q tmp
pip install ase numpy
md tmp
cd tmp
git clone https://github.com/materialsproject/fireworks.git
cd fireworks
echo " " > README.md
python setup.py install
cd ../../
rd /S /Q tmp
pip install git+https://github.com/SUNCAT-Center/CatKit.git
pip install .\requires_wheel\HTMACat-1.0.0-py3-none-any.whl