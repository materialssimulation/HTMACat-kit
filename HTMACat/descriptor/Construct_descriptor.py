from HTMACat.Extract_info import *
from HTMACat.configuration import base_info


def Construct_descriptor(poscar, feature_surf, feature_ads, feature_site, adspecies, facet="100"):
    """Construct a descriptor for a given catalytic system based on its geometry and adsorption
    properties.

    Parameters
    ----------
    poscar : str
        The POSCAR file path of the catalyst-adsorbate system.
    feature_surf : list
        A list of strings specifying the features to extract from the surface of the catalytic system. Possible options
        are 'valence electron', 'atomic radius', 'mean valence electron', and 'mean atomic radius'.
    feature_ads : list
        A list of strings specifying the features to extract from the adsorbate molecule. Possible options are 'mean
        electronegativity' and 'mean valence electron'.
    feature_site : list
        A list of strings specifying the features to extract from the binding site of the catalytic system. Possible
        options are 'mean valence electron', 'mean atomic radius', and 'binding type'.
    adspecies : list
        A list of strings specifying the adsorbate species that can bind to the catalytic system.
    facet : str, optional
        The surface facet of the catalytic system. Possible options are '100' and '111'. The default is '100'.

    Returns
    -------
    descriptor : array_like
        A 1D array containing the descriptor for the catalytic system. The descriptor is constructed by concatenating the
        following features: surface valence electron, surface atomic radius, mean valence electron of surface and subsurface
        atoms, mean atomic radius of surface and subsurface atoms, mean electronegativity of adsorbate molecule, mean valence
        electron of adsorbate molecule, mean valence electron of binding site atoms, mean atomic radius of binding site atoms,
        and binding type of binding site ('top', 'bridge', 'fcc', 'hcp', or '4-fold').
    """
    descriptor = []
    descriptor_surf = []
    ### the features of catalyst surface : surfce valence electron,surface atomic radius,
    ### surf+subsurf mean valence electron,surf+subsurf mean atomic radius
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
    # print(surfatoms_symb)
    ## If NO, The symmetry of surface atoms could be broken
    if get_symmetry_surfatoms(poscar, tol=0.3) == "NO":
        print(f"The symmetry of surface atoms of {poscar} could not be keeped!")
        # return descriptor
    else:
        feature_value_surf = Construct_descriptor_info(base_info, surfatoms_symb, feature_surf)
        feature_value_subsurf = Construct_descriptor_info(
            base_info, subsurfatoms_symb, feature_surf
        )
        # print(feature_value_surf)
        m0, n0 = np.array(feature_value_surf).shape
        m1, n1 = np.array(feature_value_subsurf).shape
        # print(m0,m1)
        ## Ignore structures without standard or integrated surface configuration
        ## if m0 not equal to m1 menas that the large surface reconstruction occurs and causes the broken of surface.
        if (feature_value_surf != []) and (m0 == m1):
            # print(feature_value_surf,feature_value_subsurf)
            feature_value = np.hstack((feature_value_surf, feature_value_subsurf))
            # print(feature_value)
            descriptor_surf_tmp = np.around(np.mean(feature_value, 0), 2)
            facet_coord = {"100": 8, "111": 9}
            descriptor_surf = np.hstack((descriptor_surf_tmp, [facet_coord.get(facet)]))
            ### the feature of adspecies and binding sites
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
            if adspecie == []:
                print(f"The molecule can not adsorb in the {poscar}!")
            elif len(adspecie) > 1:
                print(f"More than 1 adspecie are found in {poscar}!")
            elif set(adspecie).intersection(set(adspecies)):
                ## construct the descriptor of adsorbate: mean enegativity, mean valence_electron
                # print(adspecie)
                descriptor_ads = []
                # print(adspecie)
                ads = molecule(adspecie[0])
                # print(ads)
                ads_symb = ads.get_chemical_symbols()
                feature_value_ads = Construct_descriptor_info(base_info, ads_symb, feature_ads)
                descriptor_ads = np.around(np.mean(feature_value_ads, 0), 2)
                # print(descriptor_ads)

                ## construct the descriptor of site: mean valence electron, mean atomic radius, bind type
                descriptor_site = []
                typ = {None: 0, "top": 1, "bri": 2, "fcc": 3, "hcp": 3, "4-fold": 4}
                typ2 = {None: 0, "top": 0, "bri": 0, "fcc": 0, "hcp": 1, "4-fold": 0}
                site_type = np.hstack((typ.get(bind_type_symb[0]), typ2.get(bind_type_symb[0])))
                feature_value_site = Construct_descriptor_info(
                    base_info, bind_surfatoms_symb[0], feature_site
                )
                descriptor_site_tmp = np.around(np.mean(feature_value_site, 0), 2)
                descriptor_site = np.append(descriptor_site_tmp, site_type)
                # print(descriptor_site)

                descriptor = np.hstack((descriptor_surf, descriptor_ads, descriptor_site))
                # print(descriptor)
                return descriptor
            else:
                print(adspecie)
                print(f"{poscar} can not be identified!")
        else:
            print(f"Surface info of {poscar} can not be obtained ")


from HTMACat.Extract_info import *
from HTMACat.descriptor.Construct_descriptor import *
from HTMACat.Base_tools import *
import os
import numpy as np
import operator

if __name__ == "__main__":
    feature_surf = ["Valence_electron", "Atomic_radius"]
    # feature_ads=['Enegativity','Valence_electron']
    feature_ads = ["Valence_electron", "Atomic_radius"]
    feature_site = ["Valence_electron", "Atomic_radius"]
    adspecies = ["NH3", "NH2", "NH", "N", "O", "OH", "H"]
    dop_typ_all = ["1", "2", "3", "4", "1L"]
    facet = "111"
    print("----------------------------------------")
    print("Construct descriptor starts:")
    print("1st step: Get the whole 'Descriptor+Ead'")
    EnerInfo = open("adsE_radical_all", "r+")
    file_des = open("descriptor", "w+")
    file_all = open("descriptor-all", "w+")
    file_log = open("descriptor-log", "w+")
    des, des_all, des_tmp = [], [], []
    for i, Ener in enumerate(EnerInfo):
        sys = Ener.split(",")[0]
        ene = Ener.split(",")[-1].strip()
        sys_all = sys.split("_")
        dop_typ = sys_all[-3]
        specie = set(sys_all).intersection(set(adspecies))
        if specie:
            if dop_typ in dop_typ_all:
                # poscar= f'./{sys}/optmk/CONTCAR'
                poscar = f"./{sys}/CONTCAR"
                print(sys)
                descriptor = Construct_descriptor(
                    poscar, feature_surf, feature_ads, feature_site, adspecies, facet=facet
                )
                # print(descriptor)
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
    print(m, k)
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
    print("Subtrate Repeat Values End !")
    """print('---------Raw data---------') os.system('cat descriptor-all') print('---------Final
    data---------') os.system('cat descriptor')"""
