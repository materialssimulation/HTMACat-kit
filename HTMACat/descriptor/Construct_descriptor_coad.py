from HTMACat.Extract_info import *
from ase.build import molecule
from HTMACat.configuration import base_info


def Construct_descriptor_coad(poscar, feature_surf, feature_ads, feature_site):
    """Construct a descriptor for a catalyst-adsorbate system.

    Parameters
    ----------
    poscar : str
        The POSCAR file path of the catalyst-adsorbate system.
    feature_surf : list
        A list of features of the catalyst surface for constructing the descriptor.
        The features include:

        - surface valence electron
        - surface atomic radius
        - surf+subsurf mean valence electron
        - surf+subsurf mean atomic radius
    feature_ads : list
        A list of features of the adsorbate for constructing the descriptor.
        The features include:

        - mean enegativity
        - mean valence_electron
    feature_site : list
        A list of features of the binding site for constructing the descriptor.
        The features include:

        - mean valence electron
        - mean atomic radius
        - binding type

    Returns
    -------
    descriptor : numpy.ndarray
        A descriptor array for the catalyst-adsorbate system, including the following features:

        - descriptor_surf: the descriptor of the catalyst surface, including the features specified in `feature_surf`.
        - descriptor_ads: the descriptor of the adsorbate, including the features specified in `feature_ads`.
        - descriptor_site: the descriptor of the binding site, including the features specified in `feature_site`.
        - descriptor_distance: the distance of adsorbed species.

    Notes
    -----
    This function first extracts the atom types and coordinates from the POSCAR file of the catalyst-adsorbate system,
    and then distinguishes the adsorbate, surface atoms, and subsurface atoms. The features of the surface and
    subsurface atoms are calculated based on the `feature_surf` input list, while the features of the adsorbate and
    binding site are calculated based on the `feature_ads` and `feature_site` input lists, respectively. If multiple
    adsorbates are found in the system, the descriptors are constructed for each adsorbate and then combined into a
    single descriptor. If no adsorbate is found, an error message will be printed.
    """
    descriptor = []
    descriptor_surf = []
    typ = {None: 0, "top": 1, "bri": 2, "fcc": 3, "hcp": 3, "4-fold": 4}
    typ2 = {None: 0, "top": 0, "bri": 0, "fcc": 0, "hcp": 1, "4-fold": 0}

    ## the features of catalyst surface : surfce valence electron,surface atomic radius, surf+subsurf mean valence electron,surf+subsurf mean atomic radius
    # surfatoms,surfatoms_symb=distinguish_atom_binding(poscar, tol=0.03,layer='surf_atom')
    # subsurfatoms,subsurfatoms_symb=distinguish_atom_binding(poscar, tol=0.03,layer='subsurf_atom')
    (
        adatoms,
        adatoms_symb,
        surfatoms,
        surfatoms_symb,
        subsurfatoms,
        subsurfatoms_symb,
    ) = distinguish_atom_binding(poscar, tol=0.05)
    feature_value_surf = Construct_descriptor_info(base_info, surfatoms_symb, feature_surf)
    feature_value_subsurf = Construct_descriptor_info(base_info, subsurfatoms_symb, feature_surf)

    ## Ignore structures without standard or integrated surface configuration
    if feature_value_surf != []:
        feature_value = np.hstack((feature_value_surf, feature_value_subsurf))
        # print(feature_value)
        descriptor_surf = np.around(np.mean(feature_value, 0), 2)
        ## the feature of adspecies and binding sites
        (
            bind_adatoms,
            bind_adatoms_symb,
            adspecie,
            bind_type_symb,
            bind_surfatoms,
            bind_surfatoms_symb,
        ) = get_binding_adatom(poscar)
        # print(bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb)
        # print(adspecie,bind_type_symb,bind_surfatoms_symb)
        # print(bind_type_symb[0])
        # print(adspecie)
        if adspecie == []:
            print(f"The molecule can not adsorb on the {poscar}!")
        elif len(adspecie) > 1:
            # print(f"More then 1 adspecie are found in {poscar}!")
            sys_ad = poscar.split("/")[1].split("_")[-3:-1]
            if sorted(sys_ad) != sorted(adspecie):
                print(f"Adsorbed configuration reconstruct on {poscar}!")
            else:
                ## construct the descriptor of adsorbate: mean enegativity, mean valence_electron
                descriptor_ads, descriptor_ads_tmp = [], []
                for i, specie in enumerate(adspecie):
                    ads = molecule(specie)
                    ads_symb = ads.get_chemical_symbols()
                    # print(ads_symb)
                    feature_value_ads = Construct_descriptor_info(base_info, ads_symb, feature_ads)
                    # print(feature_value_ads)
                    descriptor_ads_tmp += [np.around(np.mean(feature_value_ads, 0), 2)]
                descriptor_ads = np.hstack((descriptor_ads_tmp[0], descriptor_ads_tmp[1]))

                ## construct the descriptor of site: mean valence electron, mean atomic radius, bind type
                descriptor_site = []
                site_type, descriptor_site_tmp, descriptor_site_tmp2 = [], [], []
                # typ={None:0,'top':1,'bri':2,'fcc':3,'hcp':3,'4-fold':4}
                # typ2={None:0,'top':0,'bri':0,'fcc':0,'hcp':1,'4-fold':0}
                for j, typ_site in enumerate(bind_type_symb):
                    descriptor_site_tmp += [np.hstack((typ.get(typ_site), typ2.get(typ_site)))]
                    feature_value_site = Construct_descriptor_info(
                        base_info, bind_surfatoms_symb[j], feature_site
                    )
                    descriptor_site_tmp += [np.around(np.mean(feature_value_site, 0), 2)]
                    # print(len(descriptor_site_tmp))
                    descriptor_site_tmp2 += [
                        np.hstack((descriptor_site_tmp[0], descriptor_site_tmp[1]))
                    ]
                    # print(descriptor_site)
                descriptor_site = np.hstack((descriptor_site_tmp2[0], descriptor_site_tmp2[1]))
                # print(descriptor_site)

                ## construct the descriptor:distance of adsorbed species
                dis_symb, descriptor_distance = dis_symb_matrix, dis_matrix = get_distance_adatoms(
                    poscar
                )
                descriptor = np.hstack(
                    (descriptor_surf, descriptor_ads, descriptor_site, descriptor_distance)
                )
                # print(descriptor)
                return descriptor

        else:
            ## construct the descriptor of adsorbate: mean enegativity, mean valence_electron
            sys_ad = poscar.split("/")[1].split("_")[-3:-1]
            adspecie_tmp = sys_ad
            # print(adspecie_tmp)
            descriptor_ads, descriptor_ads_tmp = [], []
            for i, specie in enumerate(adspecie_tmp):
                ads = molecule(specie)
                ads_symb = ads.get_chemical_symbols()
                # print(ads_symb)
                feature_value_ads = Construct_descriptor_info(base_info, ads_symb, feature_ads)
                # print(feature_value_ads)
                descriptor_ads_tmp += [np.around(np.mean(feature_value_ads, 0), 2)]
            descriptor_ads = np.hstack((descriptor_ads_tmp[0], descriptor_ads_tmp[1]))

            ## construct the descriptor of site: mean valence electron, mean atomic radius, bind type
            descriptor_site = []
            site_type, descriptor_site_tmp, descriptor_site_tmp2 = [], [], []

            for j, typ_site in enumerate(bind_type_symb):
                descriptor_site_tmp += [np.hstack((typ.get(typ_site), typ2.get(typ_site)))]
                feature_value_site = Construct_descriptor_info(
                    base_info, bind_surfatoms_symb[j], feature_site
                )
                descriptor_site_tmp += [np.around(np.mean(feature_value_site, 0), 2)]
                # print(len(descriptor_site_tmp))
                descriptor_site_tmp2 += [
                    np.hstack((descriptor_site_tmp[0], descriptor_site_tmp[1]))
                ]
            # print(descriptor_site_tmp2)
            descriptor_site1 = descriptor_site_tmp2[0]
            # print(descriptor_site1)

            ## add the descriptor info about the physical adsorbed site
            site2_type = [0, 0]
            site2_surfatom_feature = [0 for i in range(len(feature_site))]
            # print(np.hstack((site2_type,site2_surfatom_feature)))

            if adspecie[0] == adspecie_tmp[0]:
                descriptor_site = np.hstack((descriptor_site1, site2_type, site2_surfatom_feature))
                ## construct the descriptor:distance of adsorbed species
                descriptor_distance = [0]
                descriptor = np.hstack(
                    (descriptor_surf, descriptor_ads, descriptor_site, descriptor_distance)
                )
                # print(descriptor)
                # return descriptor
            elif adspecie[0] == adspecie_tmp[1]:
                descriptor_site = np.hstack((site2_type, site2_surfatom_feature, descriptor_site1))
                ## construct the descriptor:distance of adsorbed species
                descriptor_distance = [0]
                descriptor = np.hstack(
                    (descriptor_surf, descriptor_ads, descriptor_site, descriptor_distance)
                )
                # print(descriptor)
                # return descriptor
            else:
                print(f"Adsorbed configuration changes on {poscar}!")
    else:
        print(f"Surface info of {poscar} can not be obtained ")


from HTMACat.Extract_info import *
from HTMACat.Base_tools import *
import os
import numpy as np
import operator

if __name__ == "__main__":
    feature_surf = ["Valence_electron", "Atomic_radius"]
    # feature_surf=['Atomic_radius']
    feature_ads = ["Enegativity", "Valence_electron"]
    feature_site = ["Valence_electron", "Atomic_radius"]
    # feature_site=['Atomic_radius']
    adspecies = [
        "NH3_OH",
        "NH3_O",
        "NH2_O",
        "NH2_OH",
        "NH2_H2O",
        "NH_O",
        "NH_OH",
        "NH_H2O",
        "N_O",
        "N_OH",
        "N_H2O",
        "N_N",
        "N_NO",
        "O_O",
    ]
    print("----------------------------------------")
    print("Construct descriptor starts:")
    print("1st step: Get the whole 'Descriptor+Ead'")

    file_des = open("descriptor-coad", "w+")
    file_all = open("descriptor-coad-all", "w+")
    file_log = open("descriptor-coad-log", "w+")

    EnerInfo = open("adsE_coad", "r+")
    des, des_all, des_tmp = [], [], []
    for i, Ener in enumerate(EnerInfo):
        sys = Ener.split(",")[0]
        ene = Ener.split(",")[-1].strip()
        poscar = f"./{sys}/optmk/CONTCAR"
        # poscar= f'./{sys}/CONTCAR'

        if adspecies != []:
            sys_all = "_".join(sys.split("_")[-3:-1])

            specie = {sys_all}.intersection(set(adspecies))
            # print(set([sys_all]),set(adspecies),specie)
            if specie:
                print(sys)
                descriptor = Construct_descriptor_coad(
                    poscar, feature_surf, feature_ads, feature_site
                )
                if descriptor is None:
                    tmp = np.hstack(([sys], ["None"], [ene]))
                    for d in tmp:
                        if d == tmp[-1]:
                            file_all.write("%s\n" % d)
                        else:
                            file_all.write("%s\t" % d)
                else:
                    tmp = np.hstack(([sys], descriptor, [ene]))
                    for d in tmp:
                        if d == tmp[-1]:
                            file_all.write("%s\n" % d)
                        else:
                            file_all.write("%s\t" % d)
                    # file_all.writelines('%s\n' %tmp)
                    des_tmp += [np.hstack(([sys], descriptor, [ene]))]
        else:
            print(sys)
            descriptor = Construct_descriptor_coad(poscar, feature_surf, feature_ads, feature_site)
            if descriptor is None:
                tmp = np.hstack(([sys], ["None"], [ene]))
                for d in tmp:
                    if d == tmp[-1]:
                        file_all.write("%s\n" % d)
                    else:
                        file_all.write("%s\t" % d)
            else:
                tmp = np.hstack(([sys], descriptor, [ene]))
                for d in tmp:
                    if d == tmp[-1]:
                        file_all.write("%s\n" % d)
                    else:
                        file_all.write("%s\t" % d)
                # file_all.writelines('%s\n' %tmp)
                des_tmp += [tmp]
    # print(des_tmp[0])
    print("2nd: Extrate the repeated values")
    ###Substrate the repeated items according to the feature list
    m, n = np.array(des_tmp).shape
    # Extrate the feature value
    des_feature = [des_tmp[j][1:-1] for j in range(m)]
    # The feaure list after substrate the repeat
    des_feature_tmp = np.array(list({tuple(t) for t in des_feature}))
    k, l = des_feature_tmp.shape
    # The threshold value Ecut of Ead: if Ead > Ecut, ignore it
    Ecut = 0.25
    # Output the feature and energy
    if m > k:
        print("Subtrate Repeat Values Starts:")
        for i, d1 in enumerate(des_feature_tmp):
            for j, d2 in enumerate(des_tmp):
                if all(d1 == d2[1:-1]):
                    tmp = np.hstack((d1, d2[-1]))
                    tmp2 = np.hstack((d1, [d2[0].split("_")[0]]))
                    ## output
                    if float(tmp[-1]) < float(Ecut):
                        # output des+ene
                        for d in tmp:
                            if d == tmp[-1]:
                                file_des.write("%s\n" % d)
                            else:
                                file_des.write("%s\t" % d)
                        # output des+type: des+'Au'
                        for d in tmp2:
                            if d == tmp2[-1]:
                                file_log.write("%s\n" % d)
                            else:
                                file_log.write("%s\t" % d)
                    else:
                        print(f"Ignore the value {tmp[-1]} >= {Ecut}")

                    break
                else:
                    continue
    else:
        print("No Repeated Value")
    file_all.close()
    file_des.close()
    file_log.close()
    print("Subtrate Repeat Values End !")
    """print('---------Raw data---------') os.system('cat descriptor-all') print('---------Final
    data---------') os.system('cat descriptor')"""
