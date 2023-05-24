from HTMACat.Extract_info import *
import sys


## calculate the energy of radical
def cal_Erad(FErad, Radical):
    Erad = 0
    with open(FErad, "r+") as Eradicals:
        for i, Eradical in enumerate(Eradicals):
            radical = Eradical.split(",")[0]
            E = Eradical.split(",")[1].strip()
            # find the atom energy
            if Radical == radical:
                # print(i,radical)
                Erad = float(Erad) + float(E)
    return Erad


## calculate the energy of radical relative to atom energy: ExCyHzO=xEC+yEH+zEO
def cal_Erad_atom(FErad, Radical):
    Erad = 0
    rad = molecule(Radical)
    for i in rad.get_chemical_symbols():
        with open(FErad, "r+") as Eradicals:
            for j, Eradical in enumerate(Eradicals):
                radical = Eradical.split(",")[0]
                E = Eradical.split(",")[1].strip()
                # find the atom energy
                if i == radical:
                    # print(i,radical)
                    Erad = float(Erad) + float(E)
    return Erad


def cal_Eslab(FEslab, facet):
    E_slab = 0
    with open(FEslab, "r+") as Eslabs:
        for i, Eslab in enumerate(Eslabs):
            slab = Eslab.split(",")[0].split("_")[0:-1]
            E = Eslab.split(",")[1]
            # find the facet energy
            if facet == slab:
                # print(slab)
                E_slab = float(E)
    return E_slab


## calculate the adE with  atom energy: E(xCyHzO)ad=ECHOsurf-Esurf-xEC-yEH-zEO
def cal_Eads(Flist, FErad, FEslab, radicals, Erad_property="radical", Facet_property="all"):
    EnerInfo = open(Flist, "r+")
    Foutput = open(f"adsE_{Erad_property}_{Facet_property}", "w+")
    # num=0
    for i, ads_ener in enumerate(EnerInfo):
        # extract the facet and radical
        conf = ads_ener.split(",", 1)[0]
        conf_facet = conf.split("_")[0:-2]
        conf_radical = conf.split("_")[-2]
        # extract the energy of adsorption configuration
        Eads = ads_ener.split(",", 1)[1].strip()
        print(conf)
        if Eads:
            E_rad, E_slab = 0, 0
            if conf_radical in radicals:
                # if conf.split('_')[-3] == 'b1':
                #   num=num+1
                # else:
                #   num=num
                ## extract the energy of radical based atom energy
                if Erad_property == "atom":
                    E_rad = cal_Erad_atom(FErad, conf_radical)
                elif Erad_property == "radical":
                    E_rad = cal_Erad(FErad, conf_radical)
                else:
                    print("Erad has not the property!")
                    break

                ##  extract the energy of facet
                with open(FEslab, "r+") as Eslabs:
                    for j, Eslab in enumerate(Eslabs):
                        line1 = Eslab.split("[")[0]
                        Ftmp = line1.split(",")
                        conf_slab = Ftmp[0].split("_")[0:-1]
                        line2 = [eval(i) for i in Eslab.split("[")[1].split("]")[0].split(",")]
                        # print(Ftmp[2])
                        # print(Ftmp[2]== 'True')
                        if conf_slab == conf_facet:
                            poscar = f"./poscar/{conf}.vasp"
                            (
                                adatoms,
                                adatoms_symb,
                                surfatoms,
                                surfatoms_symb,
                                subsurfatoms,
                                subsurfatoms_symb,
                            ) = distinguish_atom_binding(
                                poscar,
                                tol=0.05,
                                base_layer=int(Ftmp[3]),
                                atoms_layer=int(float(Ftmp[4])),
                            )
                            if line2 == surfatoms_symb:
                                if Facet_property == "all":
                                    E_slab = Ftmp[1]
                                elif (Facet_property == "stable") & (Ftmp[2] == "True"):
                                    E_slab = Ftmp[1]
                                else:
                                    continue
                            else:
                                continue
                        else:
                            continue
            else:
                continue

            if E_slab:
                Eads = Extract_adsE(slab_E=E_slab, radical_E=E_rad, tot_E=Eads)
                # print(conf,Eads)
                Foutput.write(f"{conf},{Eads}\n")
            else:
                print("NO Eslab")
        else:
            print(f"{conf} is not calculated!")
    # print(f'Species: {num}')
    EnerInfo.close()
    Foutput.close()


def cal_adE_coad(Flist, FErad, FEslab, Erad_property="radical"):
    EnerInfo = open(Flist, "r+")
    Foutput = open(f"adsE_coad_{Erad_property}", "w+")

    for i, ads_ener in enumerate(EnerInfo):
        conf = ads_ener.split(",", 1)[0]
        Eads = ads_ener.split(",", 1)[1].strip()
        conf_facet = conf.split("_")[0:-3]
        conf_radicals = conf.split("_")[-3:-1]

        # if Eads and len(conf.split('_')) > 4:
        # if Eads and len(conf.split('_')) > 3:
        if Eads:
            # print(Eads,len(conf.split('_')))
            E_rad, E_slab = 0, 0
            E_slab = cal_Eslab(FEslab, conf_facet)
            if Erad_property == "atom":
                for i, conf_radical in enumerate(conf_radicals):
                    ## extract the energy of radical based atom energy
                    E_rad = float(E_rad) + float(cal_Erad_atom(FErad, conf_radical))
            elif Erad_property == "radical":
                for i, conf_radical in enumerate(conf_radicals):
                    E_rad = float(E_rad) + float(cal_Erad(FErad, conf_radical))
            else:
                print("Erad has not the property!")
                break

            Eads = Extract_adsE(slab_E=E_slab, radical_E=E_rad, tot_E=Eads)
            print(conf, Eads)
            Foutput.write(f"{conf},{Eads}\n")
        else:
            print(f"{conf} is not calculated!")
    EnerInfo.close()
    Foutput.close()


if __name__ == "__main__":
    """poscar=f'./poscar/Au_Cu_100_b1_N_7.vasp' adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsu
    rfatoms,subsurfatoms_symb=distinguish_atom_binding(poscar,tol=0.05,base_layer=4,atoms_layer=8)
    print(surfatoms_symb)"""

    Flist = "energy_list.csv"
    FErad = "/data3/home/jqyang/general-script/energy_radical"
    FEslab = "/data3/home/jqyang/general-script/energy_facet_f"
    # radicals=['NH3','NH2','NH','N','O','OH','H2O','O2','NO','N2O','N2','H']
    # radicals=['N','O']
    # radicals=['NH3','NH2','NH','N','O','OH','H','H2O','O2','NO','N2O','N2']
    radicals = ["NH3", "NH2", "NH", "N", "O", "OH", "H", "NO", "N2O", "N2"]
    # cal_Eads(Flist,FErad,FEslab,radicals,Erad_property='radical',Facet_property='all')
    cal_Eads(Flist, FErad, FEslab, radicals, Erad_property="radical", Facet_property="stable")
    """Flist="energy_list_coad.csv" FErad="/data3/home/jqyang/general-script/energy_radical"
    FEslab="/data3/home/jqyang/general-script/energy_facet"
    cal_adE_coad(Flist,FErad,FEslab,Erad_property='radical')"""
