from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.readlines()

with open("README.md", encoding="utf-8") as f:
    readme = f.read()

setup(
    name="HTMACat",
    version="1.0.5",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Operating System :: OS Independent",
    ],
    license="GPLv3",
    description="A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes.",
    long_description_content_type="text/markdown",
    long_description=readme,
    author="Jiaqiang Yang, Feifeng Wu, Zhaojie Wang, Yuxiao Lan, Haojie Li, Zhang Liu, Zhihong Zhang, Bin Shan*",
    author_email="bshan@mail.hust.edu.cn",
    url="https://stanfordbshan.github.io/HTMACat-kit/",
    keywords=["high-throughput", "python", "catalysis", "modeling"],
    python_requires=">=3.6",
    install_requires=requirements,
    packages=find_packages(),
    include_package_data=True,
    package_data={"": ["*.db", "*.txt", "*.json", "*.yaml", "*.md", "*.rst", "Element_Info"]},
    entry_points={
        "console_scripts": [
            "htmat=HTMACat.command:main",
        ]
    },
)
