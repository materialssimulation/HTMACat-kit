# -*- coding: UTF-8 -*-
from HTMACat.Extract_info import *
import os
import numpy as np


### To get the file names of elememtary reaction steps
def get_file_name(reaction):
    specie_f = reaction.split('=')[0].strip()
    specie_b = reaction.split('=')[1].strip()
    ### Construct the reaction name
    ##1.Extract the reactant molecule and type (a g s)
    specie_f_list = []
    for i in range(len(specie_f.split('+'))):
        specie_f_list += [specie_f.split('+')[i].strip()]
    specie_f_mol = []
    specie_f_typ = []
    for j, specie in enumerate(specie_f_list):
        specie_f_mol += [specie.split('(', 1)[0].strip()]
        specie_f_typ += [specie.split('(', 1)[1].split(')')[0].strip()]
    ##2.Extract the product molecule and type (a g s)
    specie_b_list = []
    for i in range(len(specie_b.split('+'))):
        specie_b_list += [specie_b.split('+')[i].strip()]
    specie_b_mol = []
    specie_b_typ = []
    for j, specie in enumerate(specie_b_list):
        specie_b_mol += [specie.split('(', 1)[0].strip()]
        specie_b_typ += specie.split('(', 1)[1].split(')')[0].strip()
    ##3.construct the file name for every reaction
    File = '+'.join(specie_f_mol) + '=' + '+'.join(specie_b_mol)
    return File


### To get info of catalytic surface
def get_surface_info(poscar, feature_surf, base_info):
    ### Extract the surface info
    descriptor_surf = []
    #base_info='/data3/home/jqyang/high-throughput/script/construction/base-function/codes/Info/Element_Info'
    #feature_surf=['Valence_electron','Atomic_radius']
    #feature_surf=['Atomic_radius']
    adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
        poscar, tol=0.05)
    #print(surfatoms_symb)
    feature_value_surf = Construct_descriptor_info(base_info, surfatoms_symb, feature_surf)
    ## Ignore structures without standard or integrated surface configuration
    if feature_value_surf != []:
        descriptor_surf = np.around(np.mean(feature_value_surf, 0), 2)
    #print(descriptor_surf)

    return descriptor_surf


### Extract the site info
def get_site_info(poscar):
    descriptor_site = []
    bind_adatoms, bind_adatoms_symb, adspecie, bind_type_symb, bind_surfatoms, bind_surfatoms_symb = get_binding_adatom(
        poscar)
    typ = {None: 0, 'top': 1, 'bri': 2, 'fcc': 3, 'hcp': 3, '4-fold': 4}
    print(bind_type_symb)
    tmp = 0
    for i in bind_type_symb:
        tmp += int(typ.get(i))
    print(tmp)
    descriptor_site = [tmp]
    return descriptor_site


### Extract the adsorbate info
def get_adsorbate_info(poscar, feature_ads, base_info):
    ### Extract adsorbate info
    bind_adatoms, bind_adatoms_symb, adspecie, bind_type_symb, bind_surfatoms, bind_surfatoms_symb = get_binding_adatom(
        poscar)
    descriptor_ads, descriptor_ads_tmp = [], []
    for i, specie in enumerate(adspecie):
        ads = molecule(specie)
        ads_symb = ads.get_chemical_symbols()
        #print(ads_symb)
        feature_value_ads = Construct_descriptor_info(base_info, ads_symb, feature_ads)
        #print(feature_value_ads)
        descriptor_ads_tmp += [np.around(np.mean(feature_value_ads, 0), 2)]
    if len(adspecie) == 2:
        descriptor_ads = np.hstack((descriptor_ads_tmp[0], descriptor_ads_tmp[1]))
    else:
        descriptor_ads = np.hstack((descriptor_ads_tmp[0], [0, 0]))
    return descriptor_ads


def get_reaction_info(react='react.log', File_adsE='adsE_coad_radical'):
    ### Extract reaction info
    ##1.Extract reaction type,reactant info,product info,barrier & enthalpy
    react_info = open(react, 'r+')
    reaction_info = react_info.readlines()
    if len(reaction_info) > 2:
        reaction_type = reaction_info[0]
        reactant_info = reaction_info[1].split(',')[0].split('=')[0]
        product_info = reaction_info[1].split(',')[1].split('=')[0]
        #print(reactant_info,product_info)
        barrier = reaction_info[-3].split(':')[-1].strip()
        enthalpy = reaction_info[-2].split(':')[-1].strip()
        #print(barrier,enthalpy)
        react_info.close()

        ##2.extract the reactant and product adsorbed energy
        Ead_reactant, Ead_product = 0, 0
        #Ead_info=open('../../surface-adsorption/adsE_coad_atom','r+')
        Ead_info = open(f'../../surface-adsorption/{File_adsE}', 'r+')
        for k, sys in enumerate(Ead_info):
            sys_name = sys.split(',')[0]
            sys_ener = sys.split(',')[1].strip()
            if sys_name == reactant_info:
                Ead_reactant = sys_ener
            elif sys_name == product_info:
                Ead_product = sys_ener
            else:
                continue
        #print(Ead_reactant,Ead_product)
        Ead_info.close()
        ##Descriptor of Ead and reaction enthalpy
        descriptor_reaction = np.hstack((Ead_reactant, Ead_product, enthalpy))
        #descriptor_reaction=np.hstack((Ead_reactant,enthalpy))
        return descriptor_reaction, barrier
    else:
        print(f'NOTE: NEB of {File} is not calculated!')
        return None, None
    react_info.close()


if __name__ == '__main__':
    ReaInfo = open("reaction-bk", "r+")
    des_barr = open('descriptor-barr', 'w+')
    des_type = open('descriptor-type', 'w+')
    base_info = '/data3/home/jqyang/high-throughput/script/construction/base-function/codes/Info/Element_Info'
    feature_surf = ['Valence_electron', 'Atomic_radius']
    #feature_surf=['Atomic_radius']
    feature_ads = ['Enegativity', 'Valence_electron']
    #File_adsE="adsE_coad_radical"
    for index, line in enumerate(ReaInfo):
        File = get_file_name(line)
        ## enter the dir
        Dir = os.getcwd()
        Element_type = Dir.split('/')[-3]
        print(f'Construction for {File} on {Element_type} system starts:')
        os.chdir(File)
        ### Extract the surface info
        descriptor_surf = get_surface_info('POSstart', feature_surf, base_info)
        ### Extract reaction info
        descriptor_reaction, barrier = get_reaction_info()
        ### Extract site info
        descriptor_site1 = get_site_info('POSstart')
        descriptor_site2 = get_site_info('POSend')
        ### Extract adsorbate info
        #descriptor_ads=get_adsorbate_info('POSstart',feature_ads,base_info)
       '''
       bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb = get_binding_adatom('POSstart')
       descriptor_ads,descriptor_ads_tmp=[],[]
       for i,specie in enumerate(adspecie):
             ads = molecule(specie)
             ads_symb=ads.get_chemical_symbols()
             #print(ads_symb)
             feature_value_ads=Construct_descriptor_info(base_info,ads_symb,feature_ads)
             #print(feature_value_ads)
             descriptor_ads_tmp += [np.around(np.mean(feature_value_ads,0),2)]
       if len(adspecie) == 2:
          descriptor_ads=np.hstack((descriptor_ads_tmp[0],descriptor_ads_tmp[1]))
       else:
          descriptor_ads=np.hstack((descriptor_ads_tmp[0],[0,0]))
       '''

        if barrier:
            ### Recombination of descriptor info
            ## Integrated descriptor
            descriptor = np.hstack(
                (descriptor_surf, descriptor_reaction, descriptor_site1, descriptor_site2))
            # descriptor=np.hstack((descriptor_surf,descriptor_reaction,descriptor_site,descriptor_ads))
            ## Descriptor+barrier
            descriptor_barrier = np.hstack((descriptor, barrier))
            ## Descriptor+element type
            descriptor_type = np.hstack((descriptor, Element_type))
            print(descriptor_barrier)
            print(descriptor_type)

            ### Output the descriptor
            #des_barr.write(f'{descriptor_barrier}')
            #des_type.write(f'{descriptor_type}')
            len_des = len(descriptor_barrier)
            print(len_des)
            for i, d in enumerate(descriptor_barrier):
                #if d == descriptor_barrier[-1]:
                if i == (int(len_des) - 1):
                    des_barr.write('%s\n' % d)
                else:
                    des_barr.write('%s\t' % d)
            for j, d in enumerate(descriptor_type):
                #if d == descriptor_type[-1]:
                if j == (int(len_des) - 1):
                    des_type.write('%s\n' % d)
                else:
                    des_type.write('%s\t' % d)

        print(f'Construction for {File} on {Element_type} system end!')
        print('-----------------------------')
        os.chdir(Dir)
    ReaInfo.close()
    des_barr.close()
    des_type.close()
