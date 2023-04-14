#!/data/jqyang/miniconda3/bin/python
# -*- coding: UTF-8 -*-
import ase
from ase.build import bulk
from ase.visualize import view
from ase.io.vasp import write_vasp
from catkit.gen.surface import SlabGenerator
from catkit.build import surface
from catkit.gen.adsorption import Builder
from catkit.gratoms import *
#from catkit.build import molecule
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import os
import numpy as np
import math
from HTMACat.Extract_info import *
from catkit.build import molecule
### 1.Construct the slab
def Construct_slab(Path_info,N_dop_bulk=[],super_cell=[3,3,1],output_bulk=False):
  # Construct lattice using the optimized lattice constant
  #if isinstance(Path_info, str):
  #   LatInfo = open(Path_info,"r+")
  #else:
  #   LatInfo = Path_info   
  #LatInfo = open(Path_info,"r+")
    slab=[]
    line=Path_info
    #print(line)
    #for index, line in enumerate(LatInfo):
    mname = line.split()[0].strip()
    mlat = line.split()[1].strip()
    mlatcon = line.split()[2].strip()
    # hcp 
    if mlat == 'hcp':
       mc = line.split(' ')[3].strip()
       mfacet = line.split(' ')[4].strip()
       La,Lb,Lc,Ld=list(mfacet)
       mbulk = bulk(mname, mlat, a=float(mlatcon), covera=float(mc)/float(mlatcon), cubic=False)
       #mbulk=mbulk*[2,1,1]
       #view(mbulk)
       if N_dop_bulk == []:
          pass
       else:
          for k in range(len(N_dop_bulk)):
              mbulk[k].symbol=N_dop_bulk[k]
       gen = SlabGenerator(
             mbulk,
             miller_index=(int(La),int(Lb),int(Lc),int(Ld)),
             layers=4,
             fixed=2,
             layer_type = 'trim',
             vacuum=8,
             standardize_bulk = True)

    else:
       # bcc fcc etc.
       mfacet = line.split(' ')[3].strip()
       La,Lb,Lc=list(mfacet)
       mbulk = bulk(mname, mlat, a=mlatcon, cubic=True)
       #view(mbulk)
       if N_dop_bulk == []:
          pass
          #return mbulk,mname
       else:
          for k in range(len(N_dop_bulk)):
              mbulk[k].symbol=N_dop_bulk[k]
       ##### generate the surfaces #####
       gen = SlabGenerator(
             mbulk,
             miller_index=(int(La),int(Lb),int(Lc)),
             layers=4,
             fixed=2,
             layer_type = 'trim',
             vacuum=8,
             standardize_bulk = True) 
    terminations = gen.get_unique_terminations()
    for i, t in enumerate(terminations):
        slab += [gen.get_slab(iterm=i)*super_cell]
        #print("%s%s_%s facet is finished!" %(mname,mfacet,i))
        #write_vasp("%s%s%s.vasp" %(mname,mfacet,i), slab, direct=True, sort=[''], vasp5=True)
         
  #LatInfo.close()
    if output_bulk :
       return mbulk,mname,N_dop_bulk
    else:
       return slab,mname,mfacet
#slab,mname,mfacet= Construct_slab(Path_info='../Info/StrucInfo')
#view(slab)

### 2.Construct the doped surface configuration
from catkit.gen.adsorption import AdsorptionSites
#from Extract_info import *
def Construct_doped_slab(Path_info,Ele_dop,Natom=1):
    # generate surface adsorption configuration
    slabs,mname,mfacet = Construct_slab(Path_info)
    slabs_dop=[]
    p1=[]
    p1_symb=[]
    for i, slab in enumerate(slabs):
        site = AdsorptionSites(slab)
        site_typ = site.get_connectivity()
        topo = site.get_topology()
        #atom_number_dop =[]
        for i in range(len(site_typ)):
            if site_typ[i] == Natom:
               atom_number_dop = topo[i] 
               p1 += [slab.get_positions()[k] for k in atom_number_dop]
               
               slb=slab.copy()
               symbol = slb.get_chemical_symbols()
               for j,item in enumerate(atom_number_dop):
                   symbol[item]=Ele_dop
               slb.set_chemical_symbols(symbol)
               p1_symb += [slb.get_chemical_symbols()[k] for k in atom_number_dop]
               slabs_dop += [slb]
               #view(slb*(2,2,1))
    return slabs_dop,mname,mfacet,p1,p1_symb
      
#slab_dop,mname,mfacet=Construct_doped_slab(Path_info='./StrucInfo',Ele_dop='Au',Natom=1)
#view(slab_dop)

### 3.Construct the doped surface configuration
def Construct_1stLayer_slab(Path_info,Ele_dop):
    # generate surface adsorption configuration
    slabs,mname,mfacet = Construct_slab(Path_info)
    slabs_dop=[]
    for i, slab in enumerate(slabs):
        surface_atoms=slab.get_surface_atoms()
        slb=slab.copy()
        for j,surf_atom in enumerate(surface_atoms) :
            slb[surf_atom].symbol=Ele_dop
        slabs_dop += [slb]
        #view(slb*(2,2,1))
    return slabs_dop,mname,mfacet,surface_atoms
#Construct_1stLayer_slab(Path_info='./StrucInfo',Ele_dop='Au')

def MolToNXGraph(m):
    '''
    convert molecule object to graph in networkx
    params:
        m: RDKit Mol object
    returns:
        G: networkx Graph object of molecule m
    '''
    G = nx.Graph()
    for i_n in range(m.GetNumAtoms()):
        G.add_node(i_n)
    bonds = [m.GetBondWithIdx(k) for k in range(len(m.GetBonds()))]
    edges = []
    for edge in bonds:
        edges.append((edge.GetBeginAtomIdx(),edge.GetEndAtomIdx()))
    G.add_edges_from(edges)
    return G

### 4.Construct single adsorption configuration 
#from catkit.gen.adsorption import AdsorptionSites
#from Extract_info import *
def Construct_single_adsorption(slabs,ads,SML):
    # generate surface adsorption configuration
    #slabs = Construct_slab()
    slab_ad=[]
    print(ads)
    for i, slab in enumerate(slabs):
        site = AdsorptionSites(slab)
        coordinates = site.get_coordinates()
        #print(site.get_symmetric_sites(unique=False, screen=True))
        #print(site.get_periodic_sites(screen=True))
        #print(len(site.get_coordinates(unique=False)))
        builder = Builder(slab)
        if SML:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads))
            '''
            AllChem.EmbedMolecule(mole)
            AllChem.MMFFOptimizeMolecule(mole)
            conf = mole.GetConformer()
            conf_coords = [conf.GetAtomPosition(k) for k in range(mole.GetNumAtoms())]
            '''
            G = MolToNXGraph(mole)
            ads_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads_use = ads_list[0]
            for ads_ in ads_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads_use = ads_
                    break
        else:
            ads_use = molecule(ads)[0]
        print(ads_use)
        for i, coord in enumerate(coordinates):
            slab_ad += [builder._single_adsorption(ads_use,bond=0,site_index=i)] ### wzj
    return slab_ad

### 4.Construct double sites adsorption configuration
def Construct_double_adsorption(slabs,ads,SML):
    # generate surface adsorption configuration
    slab_ad=[]
    for i, slab in enumerate(slabs):
        site = AdsorptionSites(slab)
        builder=Builder(slab)
        if SML:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads))
            G = MolToNXGraph(mole)
            ads_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads_use = ads_list[0]
            for ads_ in ads_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads_use = ads_
                    break
        else:
            ads_use = molecule(ads)[0]
        edges = site.get_adsorption_edges()
        for i, edge01 in enumerate(edges):
            slab_ad += [builder._double_adsorption(ads_use,bonds=[0,1],edge_index=i)]
    return slab_ad
### 5.Construct coadsoprtion configuration with mono+mono adsorption
def Construct_coadsorption_11_pristine(slabs,ads,dis_inter,ads_type,SML):
    slab_ad=[]
    #dis_inter=[3.0,5.5]
    for i, slab in enumerate(slabs):
        site01 = AdsorptionSites(slab)
        #print(site01.get_symmetric_sites())
        builder01=Builder(slab)
        coordinate01 = site01.get_coordinates()
        ##generate surface adsorption configuration
        ads1 = ads.split(' ')[0].strip()
        ads2 = ads.split(' ')[1].strip()
        if SML:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
            G = MolToNXGraph(mole)
            ads1_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads1_use = ads1_list[0]
            for ads_ in ads1_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads1_use = ads_
                    break
            mole = Chem.AddHs(Chem.MolFromSmiles(ads2))
            G = MolToNXGraph(mole)
            ads2_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads2_use = ads2_list[0]
            for ads_ in ads2_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads2_use = ads_
                    break
        else:
            ads1_use = molecule(ads1)[0]
            ads2_use = molecule(ads2)[0]
        for i, sitetype in enumerate(site01.get_symmetric_sites()):
            slab = builder01._single_adsorption(ads1_use,bond=0,site_index=i,auto_construct=True)
            coord01=site01.get_coordinates()[i]
            #view(slab)
            site02 = AdsorptionSites(slab)
            coordinates02 = site02.get_coordinates()
            for j, coord02 in enumerate(coordinates02):
                dis = np.linalg.norm(coord01-coord02, ord=2)
                if dis < float(dis_inter[0]):
                   continue
                elif dis > float(dis_inter[1]):
                   continue
                else:
                   builder02=Builder(slab)
                   slab_ad += [builder02._single_adsorption(ads2_use,bond=0,site_index=j,auto_construct=True)]
        return slab_ad
### 5.Construct coadsoprtion configuration with mono+mono adsorption
def Construct_coadsorption_11(slabs,ads,dis_inter,ads_type,SML):
    slab_ad=[]
    #dis_inter=[3.0,5.5]
    for i, slab in enumerate(slabs):
        site01 = AdsorptionSites(slab)
        #print(site01.get_symmetric_sites())
        builder01=Builder(slab)
        coordinate01 = site01.get_coordinates()
        ##generate surface adsorption configuration
        ads1 = ads.split(' ')[0].strip()
        ads2 = ads.split(' ')[1].strip()
        if SML:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
            G = MolToNXGraph(mole)
            ads1_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads1_use = ads1_list[0]
            for ads_ in ads1_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads1_use = ads_
                    break
            mole = Chem.AddHs(Chem.MolFromSmiles(ads2))
            G = MolToNXGraph(mole)
            ads2_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads2_use = ads2_list[0]
            for ads_ in ads2_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads2_use = ads_
                    break
        else:
            ads1_use = molecule(ads1)[0]
            ads2_use = molecule(ads2)[0]
        for i, sitetype in enumerate(site01.get_symmetric_sites()):
            slab = builder01._single_adsorption(ads1_use,bond=0,site_index=i,auto_construct=True)
            coord01=site01.get_coordinates()[i]
            #view(slab)
            site02 = AdsorptionSites(slab)
            coordinates02 = site02.get_coordinates()
            for j, coord02 in enumerate(coordinates02):
                dis = np.linalg.norm(coord01-coord02)
                if dis < float(dis_inter[0]):
                   continue
                elif dis > float(dis_inter[1]):
                   continue
                else:
                   builder02=Builder(slab)
                   slab_ad += [builder02._single_adsorption(ads2_use,bond=0,site_index=j,auto_construct=True)]
    #view(slab_ad[0]*(2,2,1))
    #ads_type={'NH3':[1],'NH2':[2],'NH':[2,3],'N':[2,3],'O':[2,3],'OH':[1,2,3],'NO':[1,2,3]}
    #ads_type=ads_type
    #ads_type={'NH2': [2], 'NH': [3], 'NO': [3], 'NH3': [1], 'N': [3], 'O': [3], 'OH': [3]}
    typ={None:0,'top':1,'bri':2,'fcc':3,'hcp':3,'4-fold':4}
    #view(slab_ad)
    slab_ad_final=[]
    for j, adslab in enumerate(slab_ad):
       bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(adslab)
       adspecie_tmp,bind_type_symb_tmp=[],[]
       #print(adspecie)
       for k,spe in enumerate(adspecie):
           if spe in ads_type.keys():
              adspecie_tmp += [spe]
              bind_type_symb_tmp +=[bind_type_symb[k]]
       #print(adspecie_tmp,bind_type_symb_tmp)
       if len(adspecie_tmp) < 2:
           print('Can not identify the config!')
           slab_ad_final += [adslab]
       elif typ.get(bind_type_symb_tmp[0]) in ads_type.get(adspecie_tmp[0]) and typ.get(bind_type_symb_tmp[1]) in ads_type.get(adspecie_tmp[1]):
             slab_ad_final += [adslab]
    return slab_ad_final   
### 6.Construct coadsoprtion configuration with mono+mono adsorption
def Construct_coadsorption_12(slabs,ads,dis_inter,SML):
    slab_ad=[]
    for i, slab in enumerate(slabs):
        site01 = AdsorptionSites(slab)
        #print(site01.get_symmetric_sites())
        builder01=Builder(slab)
        coordinate01 = site01.get_coordinates()
        #generate surface adsorption configuration
        ads1 = ads.split(' ')[0].strip()
        ads2 = ads.split(' ')[1].strip()
        if SML:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
            G = MolToNXGraph(mole)
            ads1_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads1_use = ads1_list[0]
            for ads_ in ads1_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads1_use = ads_
                    break
            mole = Chem.AddHs(Chem.MolFromSmiles(ads2))
            G = MolToNXGraph(mole)
            ads2_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads2_use = ads2_list[0]
            for ads_ in ads2_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads2_use = ads_
                    break
        else:
            ads1_use = molecule(ads1)[0]
            ads2_use = molecule(ads2)[0]
        for i, coord01 in enumerate(coordinate01):
            slab=builder01._single_adsorption(ads1_use,bond=0,site_index=i)
            ## after adsorbing an atoms
            #site analysis and surface builder
            site02 = AdsorptionSites(slab)
            builder02 = Builder(slab)
            #obtain the edges and ID and coordinates of all sites
            edge02 = site02.get_adsorption_edges()
            site_coord02 = site02.get_coordinates(unique=False)
            site_type02 = site02.get_periodic_sites(screen=True)
            dic_site = {}
            keys=[site_type02[i] for i in range(len(site_type02))]
            values=[site_coord02[i] for i in range(len(site_type02))]
            dic_site=dict(zip(keys,values))
            # emiliate the too near site
            for j, edge in enumerate(edge02):
                  coord10 = dic_site.get(edge[0])
                  coord11 = dic_site.get(edge[1])
                  if not coord10 is None:
                     dis1 = np.linalg.norm(coord01-coord10)
                  else:
                     dis1=dis_inter[0]+0.1
                  if not coord11 is None:
                     dis2 = np.linalg.norm(coord01-coord11)
                  else:
                     dis2 = dis_inter[1]+0.1
                  dis=min(dis1,dis2)

                  if  dis < dis_inter[0]:
                      continue
                  elif dis > dis_inter[1]:
                      continue
                  else:
                      slab_ad += [builder02._double_adsorption(ads2_use,bonds=[0,1],edge_index=j)]
    return slab_ad
### 7.Construct coadsoprtion configuration with bi+bi adsorption
def Construct_coadsorption_22(slabs,ads,dis_inter,SML):
    slab_ad=[]
    for i, slab in enumerate(slabs):
        site01 = AdsorptionSites(slab)
        #print(site01.get_symmetric_sites())
        builder01=Builder(slab)
        #coordinate01 = site01.get_coordinates()
        edge01 = site01.get_adsorption_edges()
        site_coord01 = site01.get_coordinates(unique=False)
        site_type01 = site01.get_periodic_sites(screen=True)
        dic_site01 = {}
        keys=[site_type01[i] for i in range(len(site_type01))]
        values=[site_coord01[i] for i in range(len(site_type01))]
        dic_site01=dict(zip(keys,values))
        #generate surface adsorption configuration
        ads1 = ads.split(' ')[0].strip()
        ads2 = ads.split(' ')[1].strip()
        if SML:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
            G = MolToNXGraph(mole)
            ads1_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads1_use = ads1_list[0]
            for ads_ in ads1_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads1_use = ads_
                    break
            mole = Chem.AddHs(Chem.MolFromSmiles(ads2))
            G = MolToNXGraph(mole)
            ads2_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
            ads2_use = ads2_list[0]
            for ads_ in ads2_list:
                if nx.is_isomorphic(ads_._graph, G):
                    ads2_use = ads_
                    break
        else:
            ads1_use = molecule(ads1)[0]
            ads2_use = molecule(ads2)[0]
        for i, edge01 in enumerate(edge01):
            slab = builder01._double_adsorption(ads1_use,bonds=[0,1],edge_index=i)
            coord00 = dic_site01.get(edge01[0])
            coord01 = dic_site01.get(edge01[1])
            ## after adsorbing an atoms
            #site analysis and surface builder
            site02 = AdsorptionSites(slab)
            builder02 = Builder(slab)
            #obtain the edges and ID and coordinates of all sites
            edge02 = site02.get_adsorption_edges()
            site_coord02 = site02.get_coordinates(unique=False)
            site_type02 = site02.get_periodic_sites(screen=True)
            dic_site02 = {}
            keys=[site_type02[i] for i in range(len(site_type02))]
            values=[site_coord02[i] for i in range(len(site_type02))]
            dic_site02=dict(zip(keys,values))
            ## emiliate the too near site
            for j, edge in enumerate(edge02):
                  coord10 = dic_site02.get(edge[0])
                  coord11 = dic_site02.get(edge[1])
                  if not coord10 is None:
                     dis1 = np.linalg.norm(coord01-coord10)
                  else:
                     dis1=dis_inter[0]+0.1
                  if not coord11 is None:
                     dis2 = math.hypot(coord00[0]-coord11[0],coord00[1]-coord11[1],coord00[2]-coord11[2])
                  else:
                     dis2 = dis_inter[1]+0.1
                  dis=min(dis1,dis2)

                  if  dis < dis_inter[0]:
                      continue
                  elif dis > dis_inter[1]:
                      continue
                  else:
                      ads_slab += [builder02._double_adsorption(ads2_use,bonds=[0,1],edge_index=j)]
    return slab_ad
