from HTMACat.Extract_info import *
from HTMACat.configuration import base_info


def Construct_descriptor(
    poscar, feature_surf, feature_ads, feature_site, adspecies, facet="100", label="all"
):
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

            if label == "surface":
                return descriptor_surf
            else:
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
                    site_type = np.hstack(
                        (typ.get(bind_type_symb[0]), typ2.get(bind_type_symb[0]))
                    )
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


def Construct_descriptor_info_E(ele, facet, specie, facet_dop=True):
    """Construct descriptor information based on element, crystal plane and doping species class.

     Parameters
     ----------
     ele : str
         The chemical symbol of an element.
     facet : str
         crystal face
     specie : str
         The chemical symbol of the dopant.
     facet_dop : bool, optional
         Whether the doping type is included (default True).

     Returns
     -------
     str
         The energy value in the descriptor information.

     Notes
     -----
    If 'facet_dop' is True, the doping type will be included when building descriptor information.
    """
    with open("/data3/home/jqyang/general-script/energy_single", "r+") as Efile:
        f1 = Efile.readline().split()
        EN, EO, label = [], [], []
        f2 = Efile.readlines()
        for i, sys in enumerate(f2):
            EN.append(sys.split()[0])
            EO.append(sys.split()[1])
            label.append(sys.split()[2].strip())
        if facet_dop:
            lab = f"{ele}_{facet}"
            # print(label[0],specie)
            for j, l in enumerate(label):
                if (lab == l) and (f1[0] == specie):
                    return EN[j]
                elif (lab == l) and (f1[1] == specie):
                    return EO[j]
                else:
                    continue
        else:
            lab = f"{ele}"
            for j, ls in enumerate(label):
                l = ls.split("_")[0]
                if (lab == l) and (f1[0] == specie):
                    return EN[j]
                elif (lab == l) and (f1[1] == specie):
                    return EO[j]
                else:
                    continue


from HTMACat.Extract_info import *
import os
import numpy as np
import operator


def Construct_des_module(adspecies, facet, dop_typ_all):
    """Constructs the descriptor module for a given set of adsorbates on a surface with a
    particular facet and doping type.

    Parameters
    ----------
    adspecies : list of str
        A list of adsorbate species to consider.
    facet : str
        The surface facet to consider.
    dop_typ_all : list of str
        A list of doping types to consider.

    Returns
    -------
    None
        The function writes the descriptor and energy information to an output file.
    """
    feature_surf = ["Valence_electron", "Atomic_radius"]
    # feature_ads=['Enegativity','Valence_electron']
    feature_ads = ["Valence_electron", "Atomic_radius"]
    feature_site = ["Valence_electron", "Atomic_radius"]

    # adspecies=['N']
    # facet='111'
    # dop_typ_all = ['1','2','3','1L','b1']

    file_all = open("descriptor-all-2", "w+")
    for typ in dop_typ_all:
        print("----------------------------------------")
        print("Construct descriptor {facet} {dop} starts:")
        print("1st step: Get the whole 'Descriptor+Ead'")

        EnerInfo = open(f"ads_final_{typ}", "r+")

        for i, Ener in enumerate(EnerInfo):
            sys = Ener.split(",")[0]
            ene = Ener.split(",")[-1].strip()
            sys_all = sys.split("_")
            dop_typ = sys_all[-3]
            base = sys_all[0]
            base2 = sys_all[1]
            facet3 = sys_all[2]
            specie = list(set(sys_all).intersection(set(adspecies)))
            if specie:
                if (dop_typ in dop_typ_all) and (facet == facet3):
                    # poscar= f'./{sys}/optmk/CONTCAR'
                    poscar = f"./{sys}/CONTCAR"
                    # print(sys,specie[0])
                    E1 = Construct_descriptor_info_E(base, facet, specie[0])
                    E2 = Construct_descriptor_info_E(base2, facet, specie[0])
                    descriptor = Construct_descriptor(
                        poscar,
                        feature_surf,
                        feature_ads,
                        feature_site,
                        adspecies,
                        facet=facet,
                        label="surface",
                    )
                    # print(descriptor)
                    if dop_typ == "1L":
                        dop_typ = 9
                    if descriptor is None:
                        tmp = np.hstack(([sys], ["None"], [ene], [base]))
                        for d in tmp:
                            if d == tmp[-1]:
                                file_all.write("%s\n" % d)
                            else:
                                file_all.write("%s\t" % d)

                    else:
                        tmp = np.hstack(([sys], [E1], [E2], [dop_typ], descriptor, [ene], [base]))
                        for d in tmp:
                            if d == tmp[-1]:
                                file_all.write("%s\n" % d)
                            else:
                                file_all.write("%s\t" % d)
                        # file_all.writelines('%s\n' %tmp)
        EnerInfo.close()
    file_all.close()


def Construct_des_module_2(adspecies, facet, dop_typ_all):
    """Construct descriptor module for a given adsorbate species, facet, and dopant type(s).

    Parameters
    ----------
    adspecies : str
        The adsorbate species to construct the descriptor for.
    facet : str
        The facet to construct the descriptor for.
    dop_typ_all : list
        A list of dopant types to consider.

    Returns
    ----------
    None
    This function does not return anything. It writes the constructed descriptors to a file.

    Notes
    -----
    This function reads energy information from a file named 'all' and uses the function Construct_descriptor_info_E to
    obtain information on the energy of a given adsorbate on a given facet and dopant type. The function Construct_descriptor
    is then used to construct a descriptor for a given system, and the descriptor is written to a file named
    'descriptor-{facet}-{adspecies[0]}'.
    """
    feature_surf = ["Valence_electron", "Atomic_radius"]
    feature_ads = ["Valence_electron", "Atomic_radius"]
    feature_site = ["Valence_electron", "Atomic_radius"]

    file_all = open(f"descriptor-{facet}-{adspecies[0]}", "w+")
    # for typ in dop_typ_all:

    print("----------------------------------------")
    print("Construct descriptor {facet} {dop} starts:")
    print("1st step: Get the whole 'Descriptor'")

    EnerInfo = open(f"all", "r+")
    for i, Ener in enumerate(EnerInfo):
        sys = Ener.split(".")[0]
        sys_all = sys.split("_")
        dop_typ = sys_all[-2]

        base = sys_all[0]
        base2 = sys_all[1]
        facet3 = sys_all[2]

        if (dop_typ in dop_typ_all) and (facet == facet3):
            # poscar= f'./{sys}/optmk/CONTCAR'
            poscar = f"./{sys}.vasp"
            print(sys)
            E1 = Construct_descriptor_info_E(base, facet, adspecies[0], facet_dop=True)
            E2 = Construct_descriptor_info_E(base2, facet, adspecies[0], facet_dop=False)
            descriptor = Construct_descriptor(
                poscar,
                feature_surf,
                feature_ads,
                feature_site,
                adspecies,
                facet=facet,
                label="surface",
            )

            ### construct descriptor ##
            if dop_typ == "1L":
                dop_typ = 9
            if descriptor is None:
                tmp = np.hstack(([sys], ["None"], [base]))
                for d in tmp:
                    if d == tmp[-1]:
                        file_all.write("%s\n" % d)
                    else:
                        file_all.write("%s\t" % d)

            else:
                tmp = np.hstack(([sys], [E1], [E2], [dop_typ], descriptor, [base]))
                for d in tmp:
                    if d == tmp[-1]:
                        file_all.write("%s\n" % d)
                    else:
                        file_all.write("%s\t" % d)
                # file_all.writelines('%s\n' %tmp)
    EnerInfo.close()
    file_all.close()


if __name__ == "__main__":
    # E=Construct_descriptor_info_E('Co','0001','N',facet_dop=False)
    adspecies = ["N", "O"]
    # facet='100'
    # dop_typ_all=['1','2','4','1L']
    for specie in adspecies:
        # facet='111'
        # dop_typ_all=['1','2','3','1L']
        facet = "100"
        dop_typ_all = ["1", "2", "4", "1L"]

        Construct_des_module_2([specie], facet, dop_typ_all)
