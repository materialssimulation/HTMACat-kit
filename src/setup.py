from setuptools import setup, find_packages

with open('requirements.txt', 'r') as f:
    requirements = f.readlines()

setup(
    name='HTMACat',
    version='1.0.4',
    description=
    'A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes.',
    long_description=
    'This kit it provides key tools for high-throughput design and screening of catalytic materials. The software mainly includes functional modules such as surface structure analysis and information extraction, catalytic surface and various adsorption model construction, automatic construction of primitive reaction processes, automatic extraction of computational data and automatic extraction and construction of descriptors.',
    author='Jiaqiang Yang, Feifeng Wu, Bin Shan',
    author_email='bshan@mail.hust.edu.cn',
    url='http://www.materialssimulation.com/',
    keywords='high-throughput',
    python_requires=">=3.6, <3.10",
    install_requires=requirements,
    packages=find_packages(),
    entry_points={
        'console_scripts': [  # 命令的入口
            'htmat=HTMACat.command:main',
        ]
    })
