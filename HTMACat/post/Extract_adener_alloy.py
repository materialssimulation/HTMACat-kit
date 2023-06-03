from HTMACat.Extract_info import *
import os
import csv

if __name__ == "__main__":
    # species=['N','O']
    # species=['NH3','NH2','NH','N','O','OH','H','H2O','O2','NO','N2O','N2']
    species = ["NH3", "NH2", "NH", "N", "O", "OH", "H"]
    # dop_typs=['1','1L','2','4']
    dop_typs = ["1", "1L", "2", "3"]
    Efile = "adsE_radical_all"
    facets = ["111"]
    ele = [["Ag", "Au", "Cu", "Pd", "Pt", "Ir", "Rh"]]
    # facets=['100','111']
    # ele=[['Ag','Au','Ir','Cu','Pd','Pt','Rh'],['Ag','Au','Ir','Cu','Pd','Pt','Rh','Ru']]

    for dop_typ in dop_typs:
        ### extract the adsorption energy of the most stable configuration
        for i, facet in enumerate(facets):
            for j in ele[i]:
                if j == "Ru":
                    continue
                else:
                    for k in ele[i]:
                        if j == k:
                            continue
                        else:
                            Dir = os.getcwd()
                            Dir1 = f"/data3/home/jqyang/surface/{j}{k}/{facet}/surface-adsorption/"
                            os.chdir(Dir1)
                            print(f"------{j}{k}--{facet}------")
                            fads = open(f"ads_final_{dop_typ}", "w+")
                            for specie in species:
                                # print(dop_typ)
                                print(f"{specie}")
                                spe, ene = get_adsorption_energy_stable(Efile, specie, dop_typ)
                                # poscar=f'{spe}/CONTCAR'
                                # bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(poscar)
                                fads.write(f"{spe},{ene}\n")
                                # print(spe,ene,bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb)
                            fads.close()
                            os.chdir(Dir)
        ### extract the adsorption energy of the most stable configuration
        fads0 = open(f"ads_all_alloy_{facets[0]}_{dop_typ}.csv", "w+")
        writer = csv.writer(fads0)
        writer.writerow(["sys", "Eads"])
        for specie in species:
            for i, facet in enumerate(facets):
                for j in ele[i]:
                    if j == "Ru":
                        continue
                    else:
                        for k in ele[i]:
                            if j == k:
                                continue
                            else:
                                Dir = os.getcwd()
                                Dir1 = f"/data3/home/jqyang/surface/{j}{k}/{facet}/surface-adsorption/"
                                os.chdir(Dir1)
                                fads1 = open(f"ads_final_{dop_typ}", "r+")
                                for m, ads_ene in enumerate(fads1):
                                    spe = ads_ene.split(",")[0]
                                    ene = ads_ene.split(",")[1].strip()
                                    if len(spe.split("_")) < 2:
                                        spe1 = spe.split("_")[0]
                                    else:
                                        spe1 = spe.split("_")[-2]
                                        typ = spe.split("_")[-3]

                                    if (spe1 == specie) & (typ == dop_typ):
                                        tmp = [spe, ene]
                                        writer.writerow(tmp)
                                fads1.close()
                                os.chdir(Dir)
        fads0.close()

        ### to arrange the binding energy info
        fads2 = open(f"ads_all_alloy_{facets[0]}_{dop_typ}_arr.csv", "w+")
        writer2 = csv.writer(fads2)
        writer2.writerow(["sys", "E(NH3)", "E(NH2)", "E(NH)", "E(N)", "E(O)", "E(OH)", "E(H)"])

        for i, facet in enumerate(facets):
            for j in ele[i]:
                if j == "Ru":
                    continue
                else:
                    for k in ele[i]:
                        if j == k:
                            continue
                        else:
                            Dir = os.getcwd()
                            Dir1 = f"/data3/home/jqyang/surface/{j}{k}/{facet}/surface-adsorption/"
                            os.chdir(Dir1)
                            name = "_".join([j, k, facet, dop_typ])
                            row = []
                            row += [name]

                            for specie in species:
                                fads1 = open(f"ads_final_{dop_typ}", "r+")
                                for m, ads_ene in enumerate(fads1):
                                    spe = ads_ene.split(",")[0]
                                    ene = ads_ene.split(",")[1].strip()
                                    if len(spe.split("_")) < 2:
                                        spe1 = spe.split("_")[0]
                                    else:
                                        spe1 = spe.split("_")[-2]
                                        typ = spe.split("_")[-3]
                                    if (spe1 == specie) & (typ == dop_typ):
                                        row += [ene]
                                fads1.close()
                            os.chdir(Dir)
                            # row +=['\n']
                            writer2.writerow(row)
        fads2.close()
