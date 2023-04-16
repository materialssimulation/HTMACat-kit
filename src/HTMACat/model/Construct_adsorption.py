#!/data/jqyang/miniconda3/bin/python
# -*- coding: UTF-8 -*-
from ase.visualize import view
from ase.io.vasp import write_vasp
from HTMACat.model.Construct_model import *
from HTMACat.Extract_info import *
def Construct_adsorption(Path_Info,Model):
    """
    Automated batch enumeration and construction of all possible configurations

    Parameters
    ----------
    Path_Info:str
       a filename,containing the crystallographic information for generating the slabs.
    Model:str
       a filename,containing surface doping and adsorption modeling parameters
   
    Returns
    -------
       None
          No return value, generating a vasp structure file
   """
    model = open(Model,'r')
    Ele_dop=model.readline().split()
    Natom_dop=model.readline().split()
    Ads=model.readline().split()
    SML = False
    if 'SML' == Ads[0]:
       SML = True
       Ads = Ads[1:]
    else:
       pass
    Index=model.readline().split()
    model.close()

    for i,ele_dop in enumerate(Ele_dop):
        for j,natom_dop in enumerate(Natom_dop):
            if natom_dop == '0':
               for k,ads in enumerate(Ads):
                   slabs,mname,mfacet=Construct_slab(Path_Info)
                   if Index[k] == '1':
                      slabs_ads=Construct_single_adsorption(slabs,ads,SML)
                   elif Index[k] == '2':
                      slabs_ads=Construct_double_adsorption(slabs,ads,SML)
                   for m,slab_ad in enumerate(slabs_ads):
                      write_vasp("%s_%s_%s_%s.vasp" %(mname,mfacet,ads,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print(f"{ads} {m+1} adsorption on {mname} {natom_dop} doped system are finished!")
            elif natom_dop == 'b1':
               #N_dop_bulk=[ele_dop]
               for k,ads in enumerate(Ads):
                   slabs,mname,mfacet=Construct_slab(Path_Info,N_dop_bulk=[ele_dop],super_cell=[2,2,1]) ### wzj note: Formic_acid
                   if Index[k] == '1':
                      slabs_ads=Construct_single_adsorption(slabs,ads,SML)
                   elif Index[k] == '2':
                      slabs_ads=Construct_double_adsorption(slabs,ads,SML)
                   #view(slabs_ads)
                   for m,slab_ad in enumerate(slabs_ads):
                      write_vasp("%s_%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,ads,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print(f"{ads} {m+1} adsorption on {mname}_{ele_dop} {natom_dop} doped system are finished!")
            elif natom_dop == '1L':
               for k,ads in enumerate(Ads):
                   slabs_dop,mname,mfacet,surface_atoms=Construct_1stLayer_slab(Path_Info,ele_dop)
                   if Index[k] == '1':
                      slabs_ads=Construct_single_adsorption(slabs_dop,ads,SML)
                   elif Index[k] == '2':
                      slabs_ads=Construct_double_adsorption(slabs_dop,ads,SML)
                   for m,slab_ad in enumerate(slabs_ads):
                      write_vasp("%s_%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,ads,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print(f"{ads} {m+1} adsorption on {mname}_{ele_dop} {natom_dop} doped system are finished!")
            else:
               natom_dop=int(natom_dop)
               for k,ads in enumerate(Ads):
                   slabs_dop,mname,mfacet,p1,p1_symb=Construct_doped_slab(Path_Info,ele_dop,Natom=natom_dop)
                   #print(p1,p1_symb)
                   if Index[k] == '1':
                      slabs_ads=Construct_single_adsorption(slabs_dop,ads,SML)
                      ## To choose the config whose neighbor list includes the doped atoms
                      slabs_ads_near=[]
                      #view(slabs_ads)
                      for slb in slabs_ads:
                          bind_adatoms,bind_adatoms_symb,bind_type_symb,adspecie,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(slb) 
                          #slabs_ads_near += [slb]
                          if (set(p1_symb) & set(bind_surfatoms_symb[0])):
                             slabs_ads_near += [slb]
                          
                   elif Index[k] == '2':
                      slabs_ads=Construct_double_adsorption(slabs_dop,ads,SML)
                      ## To choose the config whose neighbor list includes the doped atoms
                      slabs_ads_near=[]
                      for slb in slabs_ads:
                          bind_adatoms,bind_adatoms_symb,bind_type_symb,adspecie,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(slb)
                          if (set(p1_symb) & set(bind_surfatoms_symb[0])) or (set(p1_symb) & set(bind_surfatoms_symb[0])):
                             slabs_ads_near += [slb]
                   #view(slabs_ads)
                   #view(slabs_ads_near)
                   for m,slab_ad in enumerate(slabs_ads_near):
                      write_vasp("%s_%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,ads,m), slab_ad, direct=True, sort=[''], vasp5=True)
                   print(f"{ads} {m+1} adsorption on {mname}_{ele_dop} {natom_dop} doped system are finished!")
if __name__ == '__main__':
   Path_Info='StrucInfo'
   Model='Model' 
   with open(Path_Info,'r+') as Path:
      for i,path in enumerate(Path):
          Construct_adsorption(path,Model)

   #Construct_adsorption(Path_Info,Model)
  
