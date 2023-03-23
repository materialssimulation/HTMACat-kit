#!/data/jqyang/miniconda3/bin/python
# -*- coding: UTF-8 -*-
from ase.visualize import view
from ase.io.vasp import write_vasp
from HTMACat.model.Construct_model import *
from HTMACat.Extract_info import *
def Construct_slabs(Path_Info,Model):
    model = open(Model,'r')
    Ele_dop=model.readline().split()
    Natom_dop=model.readline().split()
    model.close()
    #Struct=open(Path_Info,'r+')
    #struct=Struct.readline().split()[0]
    #Struct.close()
    struct=Path_Info.split()[0]
    print(f'{struct} system construction')
    for i,ele_dop in enumerate(Ele_dop):
        for j,natom_dop in enumerate(Natom_dop):
            if natom_dop == '0':
               slabs,mname,mfacet=Construct_slab(Path_Info)
               for m,slab in enumerate(slabs):
                   write_vasp("%s_%s_%s.vasp" %(mname,mfacet,m), slab, direct=True, sort=[''], vasp5=True, long_format=False)
               print("pristine system are finished!")
          
            elif natom_dop == 'b1':
               if struct == ele_dop:
                  pass
               else:
                  slabs,mname,mfacet=Construct_slab(Path_Info,N_dop_bulk=[ele_dop],super_cell=[3,3,1])
                  for m,slab in enumerate(slabs):
                     write_vasp("%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,m), slab, direct=True, sort=[''], vasp5=True, long_format=False)
                  print(f"{struct}_{ele_dop} {natom_dop} doped systemare finished!")

            elif natom_dop == '1L':
               if struct == ele_dop:
                  pass
               else:
                  slabs,mname,mfacet,surface_atoms=Construct_1stLayer_slab(Path_Info,ele_dop)
                  for m,slab in enumerate(slabs):
                      write_vasp("%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,m), slab, direct=True, sort=[''], vasp5=True, long_format=False)
                  print(f"{struct}_{ele_dop} {natom_dop} doped system are finished!" )
          
            else:
              if struct == ele_dop:
                  pass
              else:    
                  natom_dop=int(natom_dop)
                  slabs,mname,mfacet,p1,p1_symb=Construct_doped_slab(Path_Info,ele_dop,Natom=natom_dop)
                  for m,slab in enumerate(slabs):
                      write_vasp("%s_%s_%s_%s_%s.vasp" %(mname,ele_dop,mfacet,natom_dop,m), slab, direct=True, sort=[''], vasp5=True, long_format=False)
                  print(f"{struct}_{ele_dop} {natom_dop} atom doped system are finished!")
            #view(slabs)
if __name__ == '__main__':
   Path_Info='StrucInfo'
   Model='Model' 
   with open(Path_Info,'r+') as Path:
      for i,path in enumerate(Path):
          Construct_slabs(path,Model)    
   #Construct_slabs(Path_Info,Model)
  
