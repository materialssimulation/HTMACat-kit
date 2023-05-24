from HTMACat.Extract_info import *
import os
import csv

if __name__ == "__main__":
    species = ["NH3", "NH2", "NH", "N", "O", "OH", "H2O", "O2", "H", "NO", "N2", "N2O"]
    Efile = "adsE_radical"
    facets = ["100", "111"]
    ele = [
        ["Ag", "Au", "Ir", "Cu", "Pd", "Pt", "Rh"],
        ["Ag", "Au", "Ir", "Cu", "Pd", "Pt", "Rh", "Ru"],
    ]
    ### extract the adsorption energy of the most stable configuration
    for i, facet in enumerate(facets):
        for j in ele[i]:
            Dir = os.getcwd()
            Dir1 = f"/data3/home/jqyang/surface/{j}/{facet}/surface-adsorption/"
            os.chdir(Dir1)
            print(f"------{j}--{facet}------")
            fads = open("ads_final", "w+")
            for specie in species:
                print(f"{specie}")
                spe, ene = get_adsorption_energy_stable(Efile, specie)
                fads.write(f"{spe},{ene}\n")
            fads.close()
            os.chdir(Dir)
    ### extract the adsorption energy of the most stable configuration
    fads0 = open("ads_all.csv", "w+")
    writer = csv.writer(fads0)
    writer.writerow(["sys", "Eads"])
    for specie in species:
        for i, facet in enumerate(facets):
            for j in ele[i]:
                Dir = os.getcwd()
                Dir1 = f"/data3/home/jqyang/surface/{j}/{facet}/surface-adsorption/"
                os.chdir(Dir1)
                fads1 = open("ads_final", "r+")
                for k, ads_ene in enumerate(fads1):
                    spe = ads_ene.split(",")[0]
                    ene = ads_ene.split(",")[1].strip()
                    if len(spe.split("_")) < 2:
                        spe1 = spe.split("_")[0]
                    else:
                        spe1 = spe.split("_")[-2]
                    # ene=ads_ene.split(',')[1].strip()
                    print(spe1, ene)
                    if spe1 == specie:
                        tmp = [spe, ene]
                        writer.writerow(tmp)
                fads1.close()
                os.chdir(Dir)
    fads0.close()
