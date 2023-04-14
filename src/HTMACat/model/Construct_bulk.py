#!/data/jqyang/miniconda3/bin/python
# -*- coding: UTF-8 -*-
from ase.visualize import view
from ase.io.vasp import write_vasp
from HTMACat.model.Construct_model import *
from HTMACat.Extract_info import *
def Construct_bulks(Path_Info,Model):
    """
    Construct bulk systems write the generated structures in VASP format.

    Parameters
    ---------- 
    Path_Info:str
       a filename,containing the crystallographic information for generating the slabs.
    Model:str
       a filename,containing surface doping and adsorption modeling parameters

    Returns
    -------
    None

    """
    model = open(Model,'r')
    Ele_dop=model.readline().split()
    Natom_dop=model.readline().split()
    model.close()
    
    struct=Path_Info.split()[0]
    print(f'{struct} system construction')
    for i,ele_dop in enumerate(Ele_dop):
        for j,natom_dop in enumerate(Natom_dop):
            if natom_dop == '0':
               slabs,mname,dop=Construct_slab(Path_Info,output_bulk=True)
               #view(slabs)
               write_vasp("%s.vasp" %(mname), slabs, direct=True, sort=[''], vasp5=True)
               print("pristine system are finished!")
          
            elif natom_dop == 'b1':
               if struct == ele_dop:
                  pass
               else:
                  slabs,mname,dop=Construct_slab(Path_Info,N_dop_bulk=[ele_dop],super_cell=[3,3,1],output_bulk=True)
                  print(dop)
                  write_vasp("%s_%s_%s.vasp" %(mname,ele_dop,natom_dop), slabs, direct=True, sort=[''], vasp5=True)
                  print(f"{struct}_{ele_dop} {natom_dop} doped systemare finished!")

if __name__ == '__main__':
   Path_Info='StrucInfo'
   Model='Model' 
   with open(Path_Info,'r+') as Path:
      for i,path in enumerate(Path):
          Construct_bulks(path,Model)    
   #Construct_slabs(Path_Info,Model)
 
