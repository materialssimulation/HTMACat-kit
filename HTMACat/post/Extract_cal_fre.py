from ase import io
from ase.constraints import FixAtoms
import os
from ase.io.vasp import write_vasp, read_vasp


def get_file_name(reaction):
    specie_f = reaction.split("=")[0].strip()
    specie_b = reaction.split("=")[1].strip()
    ### Construct the reaction name
    ##1.Extract the reactant molecule and type (a g s)
    specie_f_list = []
    for i in range(len(specie_f.split("+"))):
        specie_f_list += [specie_f.split("+")[i].strip()]
    specie_f_mol = []
    specie_f_typ = []
    for j, specie in enumerate(specie_f_list):
        specie_f_mol += [specie.split("(", 1)[0].strip()]
        specie_f_typ += [specie.split("(", 1)[1].split(")")[0].strip()]
    ##2.Extract the product molecule and type (a g s)
    specie_b_list = []
    for i in range(len(specie_b.split("+"))):
        specie_b_list += [specie_b.split("+")[i].strip()]
    specie_b_mol = []
    specie_b_typ = []
    for j, specie in enumerate(specie_b_list):
        specie_b_mol += [specie.split("(", 1)[0].strip()]
        specie_b_typ += specie.split("(", 1)[1].split(")")[0].strip()
    ##3.construct the file name for every reaction
    File = "+".join(specie_f_mol) + "=" + "+".join(specie_b_mol)
    return File


def get_adspecie_dir(react="react.log"):
    ### Extract reaction info
    ##1.Extract reaction type,reactant info,product info,barrier & enthalpy
    react_info = open(react, "r+")
    reaction_info = react_info.readlines()
    if len(reaction_info) > 2:
        reaction_type = reaction_info[0]
        reactant_info = reaction_info[1].split(",")[0].split("=")[0]
        product_info = reaction_info[1].split(",")[1].split("=")[0]
        # print(reactant_info,product_info)
        return reactant_info, product_info
    react_info.close()


if __name__ == "__main__":
    ### Extract the dir
    ReaInfo = open("reaction-bk", "r+")
    dir_list = []
    for index, line in enumerate(ReaInfo):
        File = get_file_name(line)
        Dir = os.getcwd()
        os.chdir(File)
        reactant_info, product_info = get_adspecie_dir(react="react.log")
        dir_list += [reactant_info, product_info]
        os.chdir(Dir)
    ReaInfo.close()
    # print(len(dir_list))
    # print(len(list(set(dir_list))))
    dir_list = list(set(dir_list))

    ### output ###
    Flist = open("file_list", "w+")
    for i in dir_list:
        Flist.write(f"{i}\n")
    Flist.close()
    ### construct cal dir
    dcut = 0.68
    Dir_temp = os.getcwd()
    for j in dir_list:
        print(j)
        os.chdir(f"../test/{j}")
        os.system("rm -rf fre-cal")
        # os.chdir(f'../surface-adsorption/{Dir}')
        os.system("cp -r optmk fre-cal")
        os.chdir("fre-cal")
        print(os.getcwd())
        os.system("mv CONTCAR POSCAR")
        pos = io.read("POSCAR", format="vasp")
        Mask = []
        for j in pos.get_scaled_positions():
            z = j[2]
            if float(z) > float(dcut):
                Mask += [False]
            else:
                Mask += [True]
        constraint = FixAtoms(mask=Mask)
        pos.set_constraint(constraint)
        write_vasp("POSCAR", pos, direct=True, sort=[""], vasp5=True)
        os.chdir(Dir_temp)
