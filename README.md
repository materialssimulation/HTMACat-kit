# HTMACat

# 1. Introduction

A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes, it provides key tools
for high-throughput design and screening of catalytic materials. The software mainly includes functional modules such as
surface structure analysis and information extraction, catalytic surface and various adsorption model construction,
automatic construction of primitive reaction processes, automatic extraction of computational data and automatic
extraction and construction of descriptors.The software can perform the following computational workflows: adsorption
energy calculation and analysis workflow, primitive reaction calculation and analysis workflow, high throughput
calculation and automated analysis of adsorption energy and reaction potential of catalytic primitive reaction
processes, etc.

# 2. Installation Guide

## 2.1 Environment

```requirements.txt
Python >= 3.6, <=3.9
ase >= 3.19
numpy >= 1.20.0
networkx >=2.1
scipy >= 0.1
spglib >= 0.1
ruamel.yaml > 0.15
rdkit >= 2022.3.3
typer[all] >= 0.6
```
## 2.2  Installation
Install from PyPI

```shell
pip install HTMACat

```
or install from source code
```shell
# master branch
pip install -U git+https://github.com/stanfordbshan/HTMACat-kit.git@master
# or dev branch
pip install -U git+https://github.com/stanfordbshan/HTMACat-kit.git@dev
```

# 3. Getting started

Please Visit [HTMACat-kit’s documentation](https://stanfordbshan.github.io/HTMACat-kit/) for more information.


# Author

* Jiaqiang Yang email:[jqyang_hust@hust.edu.cn](mailto:jqyang_hust@hust.edu.cn)
* Feifeng Wu email:[wufeifeng_hust@163.com](wufeifeng_hust@163.com)
* Bin Shan email:[bshan@mail.hust.edu.cn](bshan@mail.hust.edu.cn)
* Zhaojie Wang email:[yczgwangzhaojie@163.com](yczgwangzhaojie@163.com)
* Yuxiao Lan email:[husterlanxxt@163.com](husterlanxxt@163.com)
* Haojie Li email:[1197946404@qq.com](1197946404@qq.com)
* Zhang Liu email:[zhangliu@hust.edu.cn](zhangliu@hust.edu.cn)
* Zhihong Zhang email:[zhihongzh_chem@126.com](zhihongzh_chem@126.com)

# Links

* Materials Design and Nano-Manufacturing Center@HUST:http://www.materialssimulation.com/
* Pypi homepage for HTMACat project:https://pypi.org/project/HTMACat/
* HTMACat-kit’s documentation:https://stanfordbshan.github.io/HTMACat-kit/