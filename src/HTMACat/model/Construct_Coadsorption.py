#!/data/jqyang/miniconda3/bin/python
# -*- coding: UTF-8 -*-
from ase.visualize import view
from ase.io.vasp import write_vasp
from HTMACat.model.Construct_model import *
from HTMACat.Extract_info import *
def Construct_coadsorption(Path_Info,Model,spec_ads_stable):
    model = open(Model,'r+')
    Ele_dop=model.readline().split()
    Natom_dop=model.readline().split()
    Ads=model.readline().split(';')
    Index=model.readline().split(';')
    model.close()
    LatInfo = open(Path_Info,"r+")
    para=LatInfo.readline().split(' ')
    latcon = para[2]
    facet = para[-1].strip()
    if facet == '111':
       dis_min = float(latcon)*math.sqrt(2)/2*math.sqrt(3)/2 
       dis_max= float(latcon)*math.sqrt(2)/2*math.sqrt(3)
       dis_inter=[dis_min,dis_max]
    elif (facet == '100') or (facet == '0001'):
       #dis_min = float(latcon)*math.sqrt(2)/2*math.sqrt(3)/2
       #dis_max= float(latcon)*math.sqrt(2)/2*math.sqrt(3)
       #dis_inter=[dis_min,dis_max]
       #print(dis_inter)
       dis_min = float(latcon)*1
       dis_max = float(latcon)*2
       dis_inter=[dis_min,dis_max]
       #print(dis_inter)
    LatInfo.close()
    
    for i,ele_dop in enumerate(Ele_dop):
        for j,natom_dop in enumerate(Natom_dop):
         with open(Path_Info,'r+') as Path:
          for i,path in enumerate(Path):
            if natom_dop == '0':
               for k,ads in enumerate(Ads):
                   slabs,mname,mfacet=Construct_slab(path)
                   #dis_inter=[1.50,3.50]
                   ads=ads.strip() 
                   ads_symb='_'.join(ads.split(' '))
                   if Index[k].strip() == '1 1':
                      slabs_ads=Construct_coadsorption_11(slabs,ads,dis_inter,spec_ads_stable) 
                   elif Index[k].strip() == '1 2':
                      slabs_ads=Construct_coadsorption_12(slabs,ads,dis_inter)
                   elif Index[k].strip() == '2 2':
                      slabs_ads=Construct_coadsorption_22(slabs,ads,dis_inter)
                   #view(slabs_ads) 
                   for m,slab_ad in enumerate(slabs_ads):
                      #write_vasp("%s_%s_%s_%s_%s.vasp" %(mname,mfacet,natom_dop,ads_symb,m), slab_ad, direct=True, sort=[''], vasp5=True)
                      write_vasp("%s_%s_%s_%s.vasp" %(mname,mfacet,ads_symb,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print("%s %s adsorption on pristine system are finished!" %(ads_symb,m+1))
            elif natom_dop == 'b1':
               # Doped system with a ratio
               for k,ads in enumerate(Ads):
                   slabs,mname,mfacet=Construct_slab(path,N_dop_bulk=['Cu'],super_cell=[2,2,1])
                   #dis_inter=[2.50,4.50]
                   ads=ads.strip()
                   ads_symb='_'.join(ads.split(' '))
                   if Index[k].strip() == '1 1':
                      slabs_ads=Construct_coadsorption_11(slabs,ads,dis_inter,spec_ads_stable)
                   elif Index[k].strip() == '1 2':
                      slabs_ads=Construct_coadsorption_12(slabs,ads,dis_inter)
                   elif Index[k].strip() == '2 2':
                      slabs_ads=Construct_coadsorption_22(slabs,ads,dis_inter)
                   for m,slab_ad in enumerate(slabs_ads):
                      write_vasp("%s_%s_%s_b1_%s.vasp" %(mname,mfacet,ads_symb,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print("%s %s adsorption on b1 doped systemare finished!" %(ads_symb,m+1))
            elif natom_dop == '1L':
               #spec_ads_stable={'NH2': [2], 'NH': [3], 'NO': [3], 'NH3': [1], 'N': [3], 'O': [3], 'OH': [3]}
               for k,ads in enumerate(Ads):
                   slabs_dop,mname,mfacet,surface_atoms=Construct_1stLayer_slab(path,ele_dop)
                   #dis_inter=[1.50,5.50]
                   #print(dis_inter)
                   ads=ads.strip() 
                   ads_symb='_'.join(ads.split(' '))
                   if Index[k].strip() == '1 1':
                      slabs_ads=Construct_coadsorption_11(slabs_dop,ads,dis_inter,spec_ads_stable)
                   elif Index[k].strip() == '1 2':
                      slabs_ads=Construct_coadsorption_12(slabs_dop,ads,dis_inter)
                   elif Index[k].strip() == '2 2':
                      slabs_ads=Construct_coadsorption_22(slabs_dop,ads,dis_inter)
                   #view(slabs_ads)
                   for m,slab_ad in enumerate(slabs_ads):
                      write_vasp("%s_%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,ads_symb,m), slab_ad, direct=True, sort=[''], vasp5=True)
                    
                   print("%s %s adsorption on 1 layer doped system are finished!" %(ads_symb,m+1))
            else:
               natom_dop=int(natom_dop)
               for k,ads in enumerate(Ads):
                   slabs_dop,mname,mfacet,p1,p1_symb=Construct_doped_slab(path,ele_dop,Natom=natom_dop) 
                   #dis_inter=[1.50,3.50]
                   ads=ads.strip()
                   ads_symb='_'.join(ads.split(' '))
                   if Index[k].strip() == '1 1':
                      slabs_ads=Construct_coadsorption_11(slabs_dop,ads,dis_inter,spec_ads_stable)
                   elif Index[k].strip() == '1 2':
                      slabs_ads=Construct_coadsorption_12(slabs_dop,ads,dis_inter)
                   elif Index[k].strip() == '2 2':
                      slabs_ads=Construct_coadsorption_22(slabs_dop,ads,dis_inter)
                   ## To choose the config whose neighbor list includes the doped atoms
                   slabs_ads_near=[]
                   for slb in slabs_ads:
                      bind_adatoms,bind_adatoms_symb,bind_type_symb,adspecie,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(slb)
                      bind_surfatoms_symb_all=sum(bind_surfatoms_symb,[]) 
                      #print(bind_surfatoms_symb_all)
                      if set(p1_symb) & set(bind_surfatoms_symb_all):
                         slabs_ads_near += [slb]
                   #view(slabs_ads)
                   #view(slabs_ads_near)
                   for m,slab_ad in enumerate(slabs_ads_near):
                      write_vasp("%s_%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,ads_symb,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print("%s %s adsorption on %s atom doped system are finished!" %(ads_symb,m+1,natom_dop))
     
if __name__ == '__main__':
   Path_Info='StrucInfo'
   Model='Model_coad'
   Efile='../adsE'
   #spec_ads,spec_ads_stable=get_site_stable(Efile,Ecut=-0.1)
   #spec_ads_stable={'NH3':[1,2],'NH2':[2],'NH':[2,4],'N':[2,4],'O':[2,4],'OH':[2,4],'NO':[2,4],'H2O':[1],'H':[2,4]}
   spec_ads_stable={'NH3':[1],'NH2':[2],'NH':[2,4],'N':[2,4],'O':[2,4],'OH':[2,4],'NO':[2,4],'H2O':[1],'H':[2,4]}
   #spec_ads_stable={'NH2':[2],'NH':[3],'NO':[3],'NH3':[1],'N':[3],'O':[3],'OH':[3]} 
   Construct_coadsorption(Path_Info,Model,spec_ads_stable)
    
