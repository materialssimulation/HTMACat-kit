# import math
import os
from ase import io
from ase.io.vasp import write_vasp, read_vasp
from HTMACat.NEB.Construct_neb import Construct_neb
from HTMACat.NEB.Post_struc import mig_atoms
import numpy as np
from HTMACat.Extract_info import *

# Construct reaction network
ReaInfo = open("reaction", "r+")
Ads_Ener = open("EnerInfo", "a+")
Efile = "energy_list_2.csv"
# to get species involved in the reaction network
rspe1, rspe2 = Extract_reaction(Efile="reaction")
rspe = np.hstack((rspe1, rspe2))

for index, line in enumerate(ReaInfo):
    specie_f = line.split("=")[0].strip()
    specie_b = line.split("=")[1].strip()
    print("--------------------")
    print(line.strip())
    ### Operation on reactant
    """
    1.substract the reactant molecule and type (a g s)
    2.substract the dir and energy of stable configuration
    """
    ##1.substract the product molecule and type (a g s)
    specie_f_list = []
    for i in range(len(specie_f.split("+"))):
        specie_f_list += [specie_f.split("+")[i].strip()]
    specie_f_mol = []
    specie_f_typ = []
    for j, specie in enumerate(specie_f_list):
        specie_f_mol += [specie.split("(", 1)[0].strip()]
        specie_f_typ += [specie.split("(", 1)[1].split(")")[0].strip()]
    print(specie_f_mol)
    ##2.substract the dir and energy of stable configuration
    ads_dir_f = []
    energy_f = []
    EnerInfo = open(f"../surface-adsorption/{Efile}", "r+")
    for i, ads_ener in enumerate(EnerInfo):
        ads = ads_ener.split(",", 1)[0]
        ener = ads_ener.split(",", 1)[1].strip()
        # to get the reaction species info from the name of poscar
        for k in range(2, len(ads.split("_"))):
            tmp = sorted(ads.split("_")[k:-1])
            if np.array([spe in rspe for spe in tmp]).all():
                ads1 = tmp
                break
            else:
                continue

        # if np.array([spe in ads1 for spe in  specie_f_mol]).all():
        if ads1 == sorted(specie_f_mol):
            ads_dir_f += [ads]
            energy_f += [ener]
    ener_max_f = max(energy_f)
    ads_dir_max_f = ads_dir_f[energy_f.index(max(energy_f))]
    EnerInfo.close()

    ### Operation on product
    """
    1.substract the product molecule and type (a g s)
    2.substract the dir and energy of stable configuration
    """
    ##1.substract the product molecule and type (a g s)
    specie_b_list = []
    for i in range(len(specie_b.split("+"))):
        specie_b_list += [specie_b.split("+")[i].strip()]
    specie_b_mol = []
    specie_b_typ = []
    for j, specie in enumerate(specie_b_list):
        specie_b_mol += [specie.split("(", 1)[0].strip()]
        specie_b_typ += specie.split("(", 1)[1].split(")")[0].strip()
    # print(specie_b_mol)
    # print(specie_b_typ)
    energy_b = []
    ads_dir_b = []
    EnerInfo = open(f"../surface-adsorption/{Efile}", "r+")
    for i, ads_ener in enumerate(EnerInfo):
        ads = ads_ener.split(",", 1)[0]
        ener = ads_ener.split(",", 1)[1].strip()
        # to get the reaction species info from the name of poscar
        for k in range(2, len(ads.split("_"))):
            tmp = sorted(ads.split("_")[k:-1])
            if np.array([spe in rspe for spe in tmp]).all():
                ads1 = tmp
                break
            else:
                continue

        # if np.array([spe in ads1 for spe in  specie_b_mol]).all():
        if ads1 == sorted(specie_b_mol):
            ads_dir_b += [ads]
            energy_b += [ener]
    ener_max_b = max(energy_b)
    ads_dir_max_b = ads_dir_b[energy_b.index(max(energy_b))]
    EnerInfo.close()
    print(ads_dir_max_f, ads_dir_max_b)
    ## output ener info
    Rea = line.strip()
    Spef = ads_dir_max_f
    Ef = str(round(float(ener_max_f), 3))
    Speb = ads_dir_max_b
    Eb = str(round(float(ener_max_b), 3))
    E = str(round((float(ener_max_b) - float(ener_max_f)), 2))
    Ads_Ener.write(f"{Rea}:\n{Spef}={Ef} eV,{Speb}={Eb} eV,ReaEner={E} eV\n")

    ### construct dir
    Dir = os.getcwd()
    File = "+".join(specie_f_mol) + "=" + "+".join(specie_b_mol)
    path = os.path.join(Dir, File)
    os.mkdir(path, 0o777)
    Dir1 = os.path.split(Dir)[0]
    # os.system(f'cp {Dir1}/surface-adsorption/{ads_dir_max_f}/optmk/CONTCAR ./{File}/POSstart')
    # os.system(f'cp {Dir1}/surface-adsorption/{ads_dir_max_b}/optmk/CONTCAR ./{File}/POSend')
    os.system(f"cp {Dir1}/surface-adsorption/{ads_dir_max_f}/CONTCAR ./{File}/POSstart")
    os.system(f"cp {Dir1}/surface-adsorption/{ads_dir_max_b}/CONTCAR ./{File}/POSend")
    os.chdir(path)
    Reactlog = open("react.log", "w")
    Reactlog.write(f"{Rea}:\n{Spef}={Ef} eV,{Speb}={Eb} eV,ReaEner={E} eV\n")
    Reactlog.close
    ### format chansfer and delete the constraint
    initial = read_vasp("POSstart")
    final = read_vasp("POSend")
    if initial.constraints or final.constraints:
        write_vasp(
            "POSstart", initial, direct=True, sort=[""], vasp5=True, ignore_constraints=True
        )
        write_vasp("POSend", final, direct=True, sort=[""], vasp5=True, ignore_constraints=True)

    ### construct the cal file
    os.system("cp POSstart POSCAR")
    os.system("jointPOTCAR.sh")
    os.system("rm POSCAR")
    os.system("cp /data/jqyang/high-throughput/script/infile-neb/* .")
    print("calc dir is done!")
    ### construct N images for neb calc
    mig_atoms(["POSstart", "POSend"])
    Images = Construct_neb(
        struct=["POSstart", "POSend"], method_inter="linear", nimages=6, dcut=0.75
    )
    os.mkdir(f"pos", 0o777)
    for i in range(len(Images)):
        os.mkdir(f"0{i}", 0o777)
        write_vasp("0%s/POSCAR" % (i), Images[i], direct=True, sort=[""], vasp5=True)
        write_vasp("pos/POSCAR-0%s" % (i), Images[i], direct=True, sort=[""], vasp5=True)
    print(f"{int(i)-1} images are constructed")
    # print('--------------------')
    # os.system('qsub job.txt')
    os.chdir(Dir)

Ads_Ener.close()
ReaInfo.close()
print("NEB construction is done")
