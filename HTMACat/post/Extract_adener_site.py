from HTMACat.Extract_info import *
import os
import csv

if __name__ == "__main__":
    species = ["N", "O"]
    dop_typs = ["b1"]
    Efile = "ads_final_b1"

    facets = ["100"]
    ele = [["Ag", "Au", "Cu", "Pd", "Pt", "Ir", "Rh", "Ru"]]

    for dop_typ in dop_typs:
        ### extract the adsorption energy of the most stable configuration
        fads0 = open(f"ads_all_alloy_{facets[0]}_{dop_typ}_site.csv", "w+")
        writer = csv.writer(fads0)
        writer.writerow(["sys", "site", "surfatoms"])
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
                                print(f"{j}{k}-{facet}")
                                (
                                    spec_ads,
                                    spec_ads_stable,
                                    dir_list_final,
                                    spec_ads_stable_surfa,
                                ) = get_site_stable(Efile, Ecut=-0.1)
                                for k, sys in enumerate(dir_list_final):
                                    spe = sys.split("_")[-2]
                                    dop = sys.split("_")[-3]
                                    if (spe == specie) & (dop in dop_typ):
                                        tmp = [
                                            sys,
                                            spec_ads_stable.get(spe),
                                            spec_ads_stable_surfa.get(spe),
                                        ]
                                        writer.writerow(tmp)

                                os.chdir(Dir)
        fads0.close()
