# -*- coding: UTF-8 -*-
from ase.io import read
from catkit.gen.utils import to_gratoms
from collections import Counter
import os
#from catkit.build import molecule
from ase.build import molecule
from itertools import chain


### 1.Substract the reaction species
def Extract_reaction(Efile):
    ReaInfo = open(f"{Efile}", "r+")
    ad_species = []
    gas_species = []
    specie_f_mol = []
    specie_f_typ = []
    specie_b_mol = []
    specie_b_typ = []

    for index, line in enumerate(ReaInfo):
        specie_f = line.split('=')[0].strip()
        specie_b = line.split('=')[1].strip()

        ### Operation on reactant
        ##Extract the product molecule and type (a g s)
        specie_f_list = []
        for i in range(len(specie_f.split('+'))):
            specie_f_list += [specie_f.split('+')[i].strip()]
        for j, specie in enumerate(specie_f_list):
            mol = specie.split('(', 1)[0].strip()
            typ = specie.split('(', 1)[1].split(')')[0].strip()
            if typ == 'a':
                if mol not in ad_species:
                    ad_species += [mol]
            elif typ == 'g':
                if mol not in gas_species:
                    gas_species += [mol]
            specie_f_mol += [mol]
            specie_f_typ += [typ]
        #print(ad_species)
        #print(gas_species)

        ### Operation on product
        ##Exstract the product molecule and type (a g s)
        specie_b_list = []
        for i in range(len(specie_b.split('+'))):
            specie_b_list += [specie_b.split('+')[i].strip()]
        for j, specie in enumerate(specie_b_list):
            mol = specie.split('(', 1)[0].strip()
            typ = specie.split('(', 1)[1].split(')')[0].strip()
            if typ == 'a':
                if mol not in ad_species:
                    ad_species += [mol]
            elif typ == 'g':
                if mol not in gas_species:
                    gas_species += [mol]
            specie_b_mol += [mol]
            specie_b_typ += [typ]
        #print(ad_species)
        #print(gas_species)

    ReaInfo.close()
    return ad_species, gas_species


#print(Extract_reaction())


### 2.Substract energy of adsoprtion configuration
def Extract_energy(Efile, struc):
    EnerInfo = open(Efile, "r+")
    order = []
    Ener = []
    for i, ads_ener in enumerate(EnerInfo):
        ads = ads_ener.split(',', 1)[0]
        ener = ads_ener.split(',', 1)[1].strip()
        ads1 = sorted(ads.split('_', 3)[1:-1])
        num = ads.split('_', 3)[-1]
        if ads1 == sorted(struc):
            order += [num]
            Ener += [round(float(ener), 3)]
    EnerInfo.close()
    return Ener, order


#Dir = "../Info/energy_list_coad.csv"
#Ener,order=Extract_energy(Efile=Dir,struc=["NH3","OH"])
#print(Ener,order)


### 3.Extract energy of radical & slab
def Extract_energy_single(Efile, struc):
    EnerInfo = open(Efile, "r+")
    #Ener=[]
    for i, ads_ener in enumerate(EnerInfo):
        ads = ads_ener.split(',', 1)[0]
        ener = ads_ener.split(',', 1)[1].strip()
        #conf=ads.s[lit('_')[0:-2]
        if ads == struc:
            #Ener += [round(float(ener),2)]
            Ener = round(float(ener), 3)

    EnerInfo.close()
    return Ener, struc


#Dir = "../Info/energy_radical"
#Ener,struc=Exstract_energy_single(Efile=Dir,struc="NH3")
#print(Ener,struc)


### 4.Extract adsorption energy
def Extract_adsE(slab_E, radical_E, tot_E):
    ads_E = float(tot_E) - float(slab_E) - float(radical_E)
    return round(ads_E, 3)


### 5.Get the potcar file
def get_potcar(poscar, path='/data/jqyang/src/mypps/potpaw_PBE/'):
    # extract the element species:
    list_elem = []
    for i, pos in enumerate(poscar):
        struc = read(pos, format='vasp')
        symbols = struc.get_chemical_symbols()
        dic_sym = dict(Counter(symbols))
        list_elem.append(list(dic_sym.keys()))
    # if num of file, compare them
    if i > 0:
        if list_elem[0] == list_elem[1]:
            list_pot = list_elem[0]
        else:
            print("Elements in structure files are different && please correct! ")
    else:
        list_pot = list_elem[0]
    # extract the potcars file and combines them
    if os.path.exists('POTCAR'):
        os.system('rm POTCAR')
    potcar = open('POTCAR', 'w')
    for j in range(len(list_pot)):
        path_tmp = os.path.join(path, list_pot[j], 'POTCAR')
        for line in open(path_tmp):
            potcar.writelines(line)
    potcar.close()


### 6. Get the site info
from ase.neighborlist import NeighborList
from ase.neighborlist import natural_cutoffs


#from ase.build import molecule
def get_site(poscar, mole):
    struc = read(poscar, format='vasp')
    ### molecule and atom info
    mol = molecule(mole[0])
    atom_b = mole[1]

    ### POSACR structure analysis
    ### neighbor list of every atoms
    cutOff = natural_cutoffs(struc)
    nl = NeighborList(cutOff, self_interaction=False, bothways=True)
    nl.update(struc)
    matrix = nl.get_connectivity_matrix()
    ### list of every atoms : the atom order in poscar, number of binding atoms, the symbol of the atom
    adatom_num = []
    adatom_conn = []
    adatom_symb = []
    for i in range(len(struc.get_positions())):
        indices, offsets = nl.get_neighbors(i)
        #print(indices)
        if len(indices) < 5:
            #print('atoms %s is adatom'%(i))
            adatom_num += [str(i)]
            adatom_conn += [len(indices)]
            adatom_symb += [struc.get_chemical_symbols()[i]]
            #for j, offset in zip(indices, offsets):
            #print(struc.positions[j] + offset @ struc.get_cell())
    ### the dict: key= the atom order in poscar, the  number of binding atom
    dic_adatom = dict(zip(adatom_num, adatom_conn))
    # get connectivity of the atom we studied in adsorbed configuration
    if atom_b in adatom_symb:
        # the atom order of the atom we studied
        atom_b_num = struc.get_chemical_symbols().index(atom_b)

    # the binding number of the atom,list of binding atoms of the atom
    # corresponding chemical symbol atom in catalyst
    atom_b_conn_a = dic_adatom.get(str(atom_b_num))
    atom_b_conn, offsets = nl.get_neighbors(atom_b_num)
    atom_b_conn_symb = [struc.get_chemical_symbols()[i] for i in atom_b_conn if i < 36]
    #print(atom_b,atom_b_conn_a,atom_b_conn,atom_b_conn_symb)

    ### adsorbed molecule analysis
    cutOff_a = natural_cutoffs(mol)
    nl_a = NeighborList(cutOff_a, self_interaction=False, bothways=True)
    nl_a.update(mol)
    matrix = nl_a.get_connectivity_matrix()
    mol_num = []
    mol_conn = []
    mol_symb = []
    for j in range(len(mol.get_chemical_symbols())):
        indices, offsets = nl_a.get_neighbors(j)
        #atom number
        mol_num += [str(j)]
        #atom connectivity
        mol_conn += [len(indices)]
        #atom symbol
        mol_symb += [mol.get_chemical_symbols()[j]]
    dic_mol = dict(zip(mol_num, mol_conn))
    # get connectivity of atom in molecule
    if atom_b in mol_symb:
        atom_b_num = mol.get_chemical_symbols().index(atom_b)
    atom_b_conn_m = dic_mol.get(str(atom_b_num))

    ### calculate binding site type
    typ = {'0': None, '1': 'top', '2': 'bri', '3': 'hol', '4': '4-fold'}
    bind_type = typ.get(str(int(atom_b_conn_a) - int(atom_b_conn_m)))
    #print(bind_type)
    return bind_type, atom_b_conn_symb


### 7. To distinguish the surface and bulk atoms,atom binded with or not the surface
from catkit.gen.utils import get_unique_coordinates
import numpy as np


def distinguish_atom_binding(poscar, tol_layer=0.01, tol=0.05, base_layer=4, atoms_layer=9):
    '''
    To distinguish the atoms at different layer
    Parameters:
    ----------------------
    poscar: 
       atoms objects
    tol: 
       tolerance whether atoms in the same layer 
    layer: 
       layer number; eg. surface atoms:'surf_atom',subsurface atom:'subsurf_atom',
       adsorb layer:'adatom',others could be represented by the corresponding number
    base_layer:
       the layer of the substrate 
    ----------------------
    '''
    if isinstance(poscar, str):
        struct = read(poscar, format='vasp')
    else:
        struct = poscar
    #the representative Z coords of every layer
    coordinates_layer = get_unique_coordinates(struct, axis=2, tag=True, tol=tol_layer)
    coordinates = struct.get_scaled_positions()[:, 2]
    #the dict about the atoms and corresponding layers
    key_atom = []
    value_layer = []
    #print(coordinates_layer)
    for i, coord in enumerate(coordinates):
        for j, coord_layer in enumerate(coordinates_layer):
            dis = np.abs(coord - coord_layer)
            if dis < tol:
                #print(f'atom {i+1} belongs to layer {j+1}')
                key_atom += [i]
                value_layer += [j + 1]
                break
            elif j >= base_layer:
                #print(f'atom {i+1} belongs to layer {j+1}')
                key_atom += [i]
                value_layer += [j + 1]

    #if len(key_atom) != (i+1):
    #   raise ValueError('tol is too large and not to distinguish the layers; Please reduce tol!')
    dict_atom_layer = dict(zip(key_atom, value_layer))
    #print(dict_atom_layer)
    adatoms, adatoms_symb = [], []
    surfatoms, surfatoms_symb = [], []
    subsurfatoms, subsurfatoms_symb = [], []
    for k in range(len(coordinates)):
        #print(dict_atom_layer.get(k))
        if int(dict_atom_layer.get(k)) > base_layer:
            dict_atom_layer[k] = 'adatom'
            adatoms += [k]
            adatoms_symb += [struct.get_chemical_symbols()[k]]
        elif dict_atom_layer.get(k) == base_layer:
            dict_atom_layer[k] = 'surf_atom'
            surfatoms += [k]
            surfatoms_symb += [struct.get_chemical_symbols()[k]]
        elif dict_atom_layer.get(k) == (base_layer - 1):
            dict_atom_layer[k] = 'subsurf_atom'
            subsurfatoms += [k]
            subsurfatoms_symb += [struct.get_chemical_symbols()[k]]
        else:
            continue
    #print(dict_atom_layer)
    base_element = ['Au', 'Ag', 'Pd', 'Pt', 'Rh', 'Ru', 'Ir', 'Cu', 'Fe', 'Co', 'Ni', 'Zn']
    if len(surfatoms) > int(atoms_layer):
        #print('surface atoms >9')
        #print(surfatoms_symb)
        surfatoms_final, surfatoms_symb_final = [], []
        tmp = []
        for n, atom_symb in enumerate(surfatoms_symb):
            if atom_symb not in base_element:
                #print(atom_symb)
                adatoms += [surfatoms[n]]
                adatoms_symb += [atom_symb]
                #surfatoms.pop(n)
                #surfatoms_symb.pop(n)
                tmp += [n]
        for m, atom_symb in enumerate(surfatoms_symb):
            if m not in tmp:
                surfatoms_final += [surfatoms[m]]
                surfatoms_symb_final += [surfatoms_symb[m]]
        surfatoms = surfatoms_final
        surfatoms_symb = surfatoms_symb_final
        #print(surfatoms_symb)
        ### Ignore structures where Z coord of adatoms < Z coord of surfatoms
        Z_mean = np.mean([struct.get_scaled_positions()[i][2] for i in surfatoms])
        Z_adatom = np.mean([struct.get_scaled_positions()[i][2] for i in adatoms])
        if Z_adatom <= Z_mean:
            surfatoms = []
            surfatoms_symb = []
    ### Ignore the structures without standard and integrated surface configurations
    elif len(surfatoms) < int(atoms_layer):
        #raise ValueError(f'{poscar} can not been analyzed!')
        surfatoms = []
        surfatoms_symb = []
        #return adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsurfatoms,subsurfatoms_symb
    '''
    #extract the atoms and corresponding chemical symbols of specific layer
    layer_atoms=[]
    layer_atoms_symb=[] 
    for m in range(len(coordinates)):
        if dict_atom_layer.get(m) == layer:           
           layer_atoms += [m]
           layer_atoms_symb += [struct.get_chemical_symbols()[m]]
    return layer_atoms,layer_atoms_symb
    '''

    return adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb


### get the neighboring atoms fo specific adatoms
def get_atom_neigh(poscar, atom):
    if isinstance(poscar, str):
        struct = read(poscar, format='vasp')
    else:
        struct = poscar
    adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
        poscar, tol=0.03)
    struct_symb = struct.symbols
    atom_symb = chemical_symbols[atom]

    atom_index, atom_neighs = [], []
    # whether input atom occur in the bulk
    if atom_symb in struct.get_chemical_symbols():
        cutOff = natural_cutoffs(struct, mult=1.0)
        nl = NeighborList(cutOff, self_interaction=False, bothways=True)
        nl.update(struct)
        # whether input atom occur in the adatoms
        for i, adatom_symb in enumerate(adatoms_symb):
            if atom_symb == adatom_symb:
                indices, offsets = nl.get_neighbors(adatoms[i])
                for index in indices:
                    atom_index += [index]
                    atom_neighs += [struct.get_chemical_symbols()[index]]
            else:
                print(f'No {atom_symb} adatom')
    else:
        print(f'No {atom_symb} in {struct_symb}')
    return atom_index, atom_neighs


### 8. TO get atoms binding with surface among adatoms
def get_binding_adatom(poscar):
    # extract surface atoms and adsorbed atoms
    #adatoms,adatoms_symb=distinguish_atom_binding(poscar,tol=0.05,layer='adatom')
    #surf_atoms,surf_atom_symb=distinguish_atom_binding(poscar,tol=0.05,layer='surf_atom')

    # neighbor list of atoms in struct
    if isinstance(poscar, str):
        struct = read(poscar, format='vasp')
    else:
        struct = poscar
    adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
        poscar, tol=0.05)
    #print(adatoms_symb,surfatoms_symb)
    #print(struct.symbols)
    cutOff = natural_cutoffs(struct, mult=1.0)
    #print(cutOff)
    nl = NeighborList(cutOff, self_interaction=False, bothways=True)
    nl.update(struct)
    #### extract the adatoms binded with surface and corresponding surface atoms
    bind_adatoms = []
    bind_adatoms_symb = []
    bind_surfatoms = []
    bind_surfatoms_symb = []
    site_type = []
    site_type2 = []
    site_type_symb = []
    ### Extract the binded surface atoms and binded adsorbed atoms
    for i, atom in enumerate(adatoms):
        indices, offsets = nl.get_neighbors(atom)
        tmp = []
        tmp2 = []
        #print(indices)
        for index in indices:
            #print(index)
            if index in surfatoms:
                bind_adatoms += [atom]
                tmp += [index]
                tmp2 += [struct.get_chemical_symbols()[index]]
            else:
                continue
        #tmp3_symb=[struct.get_chemical_symbols()[i] for i in tmp3]
        if tmp != []:
            bind_surfatoms += [tmp]
            bind_surfatoms_symb += [tmp2]
            site_type += [len(tmp2)]
        #print(struct.get_scaled_positions()[atom][0:-1])
    bind_adatoms = list(set(bind_adatoms))
    bind_adatoms_symb = [struct.get_chemical_symbols()[i] for i in bind_adatoms]
    ### Extract the bind type
    item = []
    bind_type_symb = []
    for j, adatom in enumerate(bind_adatoms):
        p1 = struct.get_scaled_positions()[adatom][0:-1]
        item_tmp = []
        for k, atom in enumerate(subsurfatoms):
            p2 = struct.get_scaled_positions()[atom][0:-1]
            if abs(p1[0] - p2[0]) < 0.025 and abs(p1[1] - p2[1]) < 0.025:
                item_tmp += [True]
            else:
                item_tmp += [False]
        if any(item_tmp):
            item += [int(1)]
        else:
            item += [int(0)]
    #print(item)
    typ = {0: None, 1: 'top', 2: 'bri', 3: 'hol', 4: '4-fold'}
    for m, bind in enumerate(site_type):
        bind_type = typ.get(bind)
        if bind_type == 'hol' and item[m] == 0:
            bind_type_symb += ['fcc']
        elif bind_type == 'hol' and item[m] == 1:
            bind_type_symb += ['hcp']
        else:
            bind_type_symb += [bind_type]

    ### Extract the adsorbed species
    adspecie = []
    for i, atom in enumerate(bind_adatoms):
        indices, offsets = nl.get_neighbors(atom)
        tmp3 = [atom]
        #print(indices)
        for index in indices:
            if (index in surfatoms) or (index in subsurfatoms):
                continue
            else:
                tmp3 += [index]
        #print(tmp3)
        tmp3_symb = ''.join(list([struct.get_chemical_symbols()[i] for i in tmp3]))
        Ele = list(dict(Counter(tmp3_symb)).keys())
        Num = list(dict(Counter(tmp3_symb)).values())
        #print(tmp3_symb)
        #print(Ele)
        mol = '*'
        for j, E in enumerate(Ele):
            mol = mol + E
            if Num[j] == 1:
                continue
            else:
                mol = mol + str(Num[j])
        if mol == '*OH2':
            mol = '*H2O'
        else:
            adspecie += [mol.split('*')[1]]
        #adspecie+=[''.join(list(chain.from_iterable(zip(Ele,Num))))]
        #mol=molecule(tmp3_symb)[0]
        #print(mol.symbol)
    return bind_adatoms, bind_adatoms_symb, adspecie, bind_type_symb, bind_surfatoms, bind_surfatoms_symb


### 9. To get the distance between adatoms
from ase.geometry import get_distances
from catkit.gen import defaults


def get_distance_adatoms(poscar, tol=0.1):
    if isinstance(poscar, str):
        struct = read(poscar, format='vasp')
    else:
        struct = poscar
    #struct=read(poscar,format='vasp')
    bind_adatoms, bind_adatoms_symb, adspecie, bind_type_symb, bind_surfatoms, bind_surfatoms_symb = get_binding_adatom(
        poscar)
    radii = defaults.get('radii')
    radii_adatom = [radii[struct.numbers[i]] for i in bind_adatoms]
    dis_matrix = []  # distcance
    dis_symb_matrix = []  # symbol corresponding distance
    if len(bind_adatoms) < 2:
        raise ValueError('Only 1 atom and No distance')
    ## distance when there are 2 adatoms
    #elif len(bind_adatoms) == 2:
    #bond_distance=sum(radii_adatom)
    #print(bond_distance)
    #dis = get_distances(struct.get_positions()[bind_adatoms[0]],struct.get_positions()[bind_adatoms[1]])[1][0][0]
    #print(dis)
    #if abs(dis-bond_distance) < tol:
    #   print('The dis may below the possible bond length')
    #else:
    #print(round(dis,3))
    #return round(dis,3)
    ## distance when there are above 2 adatoms
    else:
        for i in range(len(bind_adatoms)):
            if (i + 1) > (len(bind_adatoms) - 1):
                continue
            else:
                dis_symb = '-'.join([bind_adatoms_symb[i], bind_adatoms_symb[i + 1]])
                bond_distance = radii_adatom[i] + radii_adatom[i + 1]
                ## get minmum distance of adsorbed atoms
                p_all = struct.get_scaled_positions()
                p1_scaled = struct.get_scaled_positions()[bind_adatoms[i]]
                p2_scaled = struct.get_scaled_positions()[bind_adatoms[i + 1]]
                # replace the minus coord
                struct.set_scaled_positions(p_all)
                #print(p1_scaled[0],p2_scaled[0])
                if p1_scaled[0] - p2_scaled[0] >= 0.50:
                    p_all[bind_adatoms[i + 1]][0] = p2_scaled[0] + 1
                    struct.set_scaled_positions(p_all)
                elif p1_scaled[0] - p2_scaled[0] <= -0.50:
                    p_all[bind_adatoms[i]][0] = p1_scaled[0] + 1
                    struct.set_scaled_positions(p_all)
                if p1_scaled[1] - p2_scaled[1] >= 0.50:
                    p_all[bind_adatoms[i + 1]][1] = p2_scaled[1] + 1
                    struct.set_scaled_positions(p_all)
                elif p1_scaled[1] - p2_scaled[1] <= -0.50:
                    p_all[bind_adatoms[i]][1] = p1_scaled[1] + 1
                    struct.set_scaled_positions(p_all)
                dis = get_distances(struct.get_positions()[bind_adatoms[i]],
                                    struct.get_positions()[bind_adatoms[i + 1]])[1][0][0]
                #print(struct.get_positions()[bind_adatoms[i]],struct.get_positions()[bind_adatoms[i+1]])
                #print(get_distances(struct.get_positions()[bind_adatoms[i]],struct.get_positions()[bind_adatoms[i+1]]))

                ## ignore the configuraton with too near distance
                if abs(dis - bond_distance) < tol:
                    print('The dis may below the possible bond length')
                else:
                    dis_symb_matrix += [dis_symb]
                    dis_matrix += [round(dis, 3)]
    return dis_symb_matrix, dis_matrix


### 10. To get the minimum distance between two group
from ase.geometry import get_distances


def get_min_distance_group(struct, group1, group2):
    p1_scaled = [struct.get_scaled_positions()[i] for i in group1]
    p2_scaled = [struct.get_scaled_positions()[i] for i in group2]
    dis_matrix = get_distances(p1_scaled, p2_scaled)
    return dis_matrix


### 11. To peprocess the raw info for descritor construction
def Construct_descriptor_info(raw_file, atoms, feature):
    comp = open(raw_file, "r+")
    index_line, index_row, matrix_info = [], [], []
    for i, com in enumerate(comp):
        line_tmp = com.split()
        line_tmp[-1] = com.split()[-1].strip()
        if i == 0:
            tmp1, tmp2 = [], []
            for j, ltmp in enumerate(line_tmp):
                tmp1 += [ltmp.strip()]
                tmp2 += [j - 1]
            index_row = dict(zip(tmp1, tmp2))
        else:
            index_line += [line_tmp[0]]
            matrix_info += [line_tmp[1:]]
    comp.close()
    feature_value = []
    for i, ele in enumerate(atoms):
        line = index_line.index(ele)
        tmp = []
        for j, des in enumerate(feature):
            row = index_row.get(des)
            tmp += [float(matrix_info[line][row])]
        feature_value += [tmp]
    #return index_line,index_row,matrix_info,feature_value
    return feature_value


### 12. Extrate the based info:atomic_numbers, atomic_names, atomic_masses, covalent_radii
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii, vdw_radii


def Extract_atomic_info(atoms):
    numbers_atoms = [atomic_numbers[i] for i in atoms]
    names_atoms = [atomic_names[i] for i in numbers_atoms]
    radii_atoms = [covalent_radii[i] for i in numbers_atoms]
    mass_atoms = [atomic_masses[i] for i in numbers_atoms]
    vdw_radii_atoms = [vdw_radii[i] for i in numbers_atoms]
    return names_atoms, radii_atoms, mass_atoms, vdw_radii_atoms


### 13. Extract the stable adsorption type for single molecule or radical adsorption based on adsorption energy
def get_site_stable(Efile, Ecut=-0.1):
    EnerInfo = open(Efile, "r+")
    species_ads = []
    typ = {None: 0, 'top': 1, 'bri': 2, 'fcc': 3, 'hcp': 3, '4-fold': 4}
    for i, ads_ener in enumerate(EnerInfo):
        species_ads += [ads_ener.split(',')[0].split('_')[-2]]
    specie_ads = set(species_ads)
    EnerInfo.close()

    spec, ads_type, ads_type_stable = [], [], []
    dir_list_final = []
    ads_type_surfa, ads_type_stable_surfa = [], []

    for specie in specie_ads:
        ads_type_tmp, energy, dir_list = [], [], []
        #ads_type_surfa_tmp=[]
        spec += [specie]
        EnerInfo = open(Efile, "r+")
        for j, ads_ener in enumerate(EnerInfo):
            ads = ads_ener.split(',')[0]
            ads1 = ads.split('_')[-2]
            ener = ads_ener.split(',')[-1].strip()
            #print(ads1,specie,ads1==specie)
            if ads1 == specie and float(ener) <= float(Ecut):
                bind_adatoms, bind_adatoms_symb, adspecie, bind_type_symb, bind_surfatoms, bind_surfatoms_symb = get_binding_adatom(
                    f'{ads}/CONTCAR')
                #bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(f'{ads}/optmk/CONTCAR')
                #print(bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb)
                if bind_type_symb == []:
                    bind_type_symb = ['top']
                ads_type_tmp += [typ.get(bind_type_symb[0])]
                energy += [ener]
                dir_list += [ads_ener.split(',')[0]]
                ads_type_surfa += [bind_surfatoms_symb[0]]
        EnerInfo.close()

        ads_type += [list(set(ads_type_tmp))]
        ads_type_stable += [[ads_type_tmp[energy.index(max(energy))]]]
        dir_list_final += [dir_list[energy.index(max(energy))]]
        ads_type_stable_surfa += [ads_type_surfa[energy.index(max(energy))]]
    spec_ads = dict(zip(spec, ads_type))
    spec_ads_stable = dict(zip(spec, ads_type_stable))
    spec_ads_stable_surfa = dict(zip(spec, ads_type_stable_surfa))
    #print(ads_type_stable_surfa)
    return spec_ads, spec_ads_stable, dir_list_final, spec_ads_stable_surfa


###13-2. Extract the stable adsorption type for single molecule or radical adsorption based on the calculated energy
def get_site_stable_energy(Efile, Ecut=0.0):
    EnerInfo = open(Efile, "r+")
    species_ads = []
    typ = {None: 0, 'top': 1, 'bri': 2, 'fcc': 3, 'hcp': 3, '4-fold': 4}
    for i, ads_ener in enumerate(EnerInfo):
        species_ads += [ads_ener.split(',')[0].split('_')[-2]]
    specie_ads = set(species_ads)
    EnerInfo.close()

    spec, ads_type, ads_type_stable = [], [], []
    dir_list_final = []
    for specie in specie_ads:
        ads_type_tmp, energy, dir_list = [], [], []
        spec += [specie]
        EnerInfo = open(Efile, "r+")
        for j, ads_ener in enumerate(EnerInfo):
            ads = ads_ener.split(',')[0]
            ads1 = ads.split('_')[-2]
            ener = ads_ener.split(',')[-1].strip()
            #print(ads1,specie,ads1==specie)
            if ads1 == specie and float(ener) < float(Ecut):
                bind_adatoms, bind_adatoms_symb, adspecie, bind_type_symb, bind_surfatoms, bind_surfatoms_symb = get_binding_adatom(
                    f'{ads}/CONTCAR')

                if bind_type_symb == []:
                    bind_type_symb = ['top']
                ads_type_tmp += [typ.get(bind_type_symb[0])]
                energy += [ener]
                dir_list += [ads_ener.split(',')[0]]
        EnerInfo.close()

        ads_type += [list(set(ads_type_tmp))]
        ads_type_stable += [[ads_type_tmp[energy.index(max(energy))]]]
        dir_list_final += [dir_list[energy.index(max(energy))]]
    spec_ads = dict(zip(spec, ads_type))
    spec_ads_stable = dict(zip(spec, ads_type_stable))
    return spec_ads, spec_ads_stable, dir_list_final


### get the energy of  most stable configurationi for single radical
def get_adsorption_energy_stable(Efile, specie, dop_typ):
    EnerInfo = open(Efile, "r+")
    dir_spe, dir_ene = [], []
    for i, ads_ener in enumerate(EnerInfo):
        ads = ads_ener.split(',')[0]
        ads1 = ads.split('_')[-2]
        dop = ads.split('_')[-3]
        ener = ads_ener.split(',')[-1].strip()
        #if (ads1 == specie) & (dop in dop_typ):
        if (ads1 == specie) & (dop == dop_typ):
            dir_spe += [ads]
            dir_ene += [float(ener)]
    EnerInfo.close()
    if dir_ene:
        spe = dir_spe[dir_ene.index(min(dir_ene))]
        ene = dir_ene[dir_ene.index(min(dir_ene))]
        return spe, ene
    else:
        return specie, None


### get the file name of reaction dir
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


### calculate the energy of radicals
def cal_Erad(FErad, Radical):
    Erad = 0
    with open(FErad, 'r+') as Eradicals:
        for i, Eradical in enumerate(Eradicals):
            radical = Eradical.split(',')[0]
            E = Eradical.split(',')[1].strip()
            # find the atom energy
            if Radical == radical:
                #print(i,radical)
                Erad = float(Erad) + float(E)
    return Erad


## calculate the energy of radical relative to atom energy: ExCyHzO=xEC+yEH+zEO
def cal_Erad_atom(FErad, Radical):
    Erad = 0
    rad = molecule(Radical)
    for i in rad.get_chemical_symbols():
        with open(FErad, 'r+') as Eradicals:
            for j, Eradical in enumerate(Eradicals):
                radical = Eradical.split(',')[0]
                E = Eradical.split(',')[1].strip()
                # find the atom energy
                if i == radical:
                    #print(i,radical)
                    Erad = float(Erad) + float(E)
    return Erad


def cal_Eslab(FEslab, facet):
    E_slab = 0
    with open(FEslab, 'r+') as Eslabs:
        for i, Eslab in enumerate(Eslabs):
            slab = Eslab.split(',')[0].split('_')[0:-1]
            E = Eslab.split(',')[1]
            # find the facet energy
            if facet == slab:
                #print(slab)
                E_slab = float(E)
    return E_slab


## calculate the adE with  atom energy: E(xCyHzO)ad=ECHOsurf-Esurf-xEC-yEH-zEO
def cal_Eads(Flist, FErad, FEslab, radicals, Erad_property='radical', Facet_property='all'):
    EnerInfo = open(Flist, 'r+')
    Foutput = open(f'adsE_{Erad_property}_{Facet_property}', 'w+')
    #num=0
    for i, ads_ener in enumerate(EnerInfo):
        # extract the facet and radical
        conf = ads_ener.split(',', 1)[0]
        conf_facet = conf.split('_')[0:-2]
        conf_radical = conf.split('_')[-2]
        # extract the energy of adsorption configuration
        Eads = ads_ener.split(',', 1)[1].strip()
        #print(conf)
        if Eads:
            E_rad, E_slab = 0, 0
            if conf_radical in radicals:
                #if conf.split('_')[-3] == 'b1':
                #   num=num+1
                #else:
                #   num=num
                ## extract the energy of radical based atom energy
                if Erad_property == 'atom':
                    E_rad = cal_Erad_atom(FErad, conf_radical)
                elif Erad_property == 'radical':
                    E_rad = cal_Erad(FErad, conf_radical)
                else:
                    print('Erad has not the property!')
                    break

                ##  extract the energy of facet
                with open(FEslab, 'r+') as Eslabs:
                    for j, Eslab in enumerate(Eslabs):
                        line1 = Eslab.split('[')[0]
                        Ftmp = line1.split(',')
                        conf_slab = Ftmp[0].split('_')[0:-1]
                        line2 = [eval(i) for i in Eslab.split('[')[1].split(']')[0].split(',')]
                        #print(Ftmp[2])
                        #print(Ftmp[2]== 'True')
                        if conf_slab == conf_facet:
                            poscar = f'./poscar/{conf}.vasp'
                            adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
                                poscar,
                                tol=0.05,
                                base_layer=int(Ftmp[3]),
                                atoms_layer=int(float(Ftmp[4])))
                            if line2 == surfatoms_symb:
                                if Facet_property == 'all':
                                    E_slab = Ftmp[1]
                                elif (Facet_property == 'stable') & (Ftmp[2] == 'True'):
                                    E_slab = Ftmp[1]
                                else:
                                    continue
                            else:
                                continue
                        else:
                            continue
            else:
                continue

            if E_slab:
                Eads = Extract_adsE(slab_E=E_slab, radical_E=E_rad, tot_E=Eads)
                #print(conf,Eads)
                Foutput.write(f'{conf},{Eads}\n')
            else:
                print('NO Eslab')
        else:
            print(f'{conf} is not calculated!')
    #print(f'Species: {num}')
    EnerInfo.close()
    Foutput.close()


def cal_adE_coad(Flist, FErad, FEslab, Erad_property='radical'):
    EnerInfo = open(Flist, 'r+')
    Foutput = open(f'adsE_coad_{Erad_property}', 'w+')

    for i, ads_ener in enumerate(EnerInfo):
        conf = ads_ener.split(',', 1)[0]
        Eads = ads_ener.split(',', 1)[1].strip()
        conf_facet = conf.split('_')[0:-3]
        conf_radicals = conf.split('_')[-3:-1]

        #if Eads and len(conf.split('_')) > 4:
        #if Eads and len(conf.split('_')) > 3:
        if Eads:
            #print(Eads,len(conf.split('_')))
            E_rad, E_slab = 0, 0
            E_slab = cal_Eslab(FEslab, conf_facet)
            if Erad_property == 'atom':
                for i, conf_radical in enumerate(conf_radicals):
                    ## extract the energy of radical based atom energy
                    E_rad = float(E_rad) + float(cal_Erad_atom(FErad, conf_radical))
            elif Erad_property == 'radical':
                for i, conf_radical in enumerate(conf_radicals):
                    E_rad = float(E_rad) + float(cal_Erad(FErad, conf_radical))
            else:
                print('Erad has not the property!')
                break

            Eads = Extract_adsE(slab_E=E_slab, radical_E=E_rad, tot_E=Eads)
            print(conf, Eads)
            Foutput.write(f'{conf},{Eads}\n')
        else:
            print(f'{conf} is not calculated!')
    EnerInfo.close()
    Foutput.close()


### TO extract the struct info and energy of specific slab and
def Extract_slab_info_1(Flist, facet):
    E_slab, surf = [], []
    with open(Flist, 'r+') as Eslabs:
        for i, Eslab in enumerate(Eslabs):
            conf = Eslab.split(',')[0]
            slab = Eslab.split(',')[0].split('_')[0:-1]
            E = Eslab.split(',')[1]
            ## TO extract the energy and surface info
            if facet == slab:
                E_slab += [float(E)]
                # TO extract struct info: layers and atom numbers of layer
                poscar = f'./poscar/{conf}.vasp'
                struct = read(poscar, format='vasp')
                layer = len(get_unique_coordinates(struct, axis=2, tag=True, tol=0.01))
                atoms_layer = int(len(struct.get_chemical_symbols())) / int(layer)
                adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
                    poscar, tol=0.05, base_layer=layer, atoms_layer=atoms_layer)
                surf += [surfatoms_symb]
    ## TO find out the most stable facet
    if len(E_slab) == 1:
        pass
    elif len(E_slab) > 1:
        E_stable = min(E_slab)
        surf_stable = surf[E_slab.index(E_stable)]

    return E_slab, surf, [E_stable], [surf_stable], [layer, atoms_layer]


###TO Extract the struct info and energy of slabs
def Extract_slab_info_2(Flist, layer=4):
    Fslab = open('energy_facet_f', 'w+')
    with open(Flist, 'r+') as Eslabs:
        for i, Eslab in enumerate(Eslabs):
            conf = Eslab.split(',')[0]
            #TO extract energy
            E = Eslab.split(',')[1]
            E_slab = float(E)
            # TO extract struct info: layers and atom numbers of layer
            poscar = f'./poscar/{conf}.vasp'
            struct = read(poscar, format='vasp')
            layer = len(get_unique_coordinates(struct, axis=2, tag=True, tol=0.01))
            atoms_layer = int(len(struct.get_chemical_symbols())) / int(layer)
            adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
                poscar, tol=0.05, base_layer=layer, atoms_layer=atoms_layer)
            surf = surfatoms_symb

            Fslab.write(f'{conf},\t{E_slab},\t{surf}\n')
    Fslab.close()


###To identify whether the symmetry is keeped###
def get_symmetry_surfatoms(poscar, tol=0.3):
    if isinstance(poscar, str):
        struct = read(poscar, format='vasp')
    else:
        struct = poscar
    #struct=read(poscar,format='vasp')
    adatoms, adatoms_symb, surfatoms, surfatoms_symb, subsurfatoms, subsurfatoms_symb = distinguish_atom_binding(
        poscar, tol=0.05)
    #print(surfatoms)

    dis_matrix = []  # distcance
    if len(surfatoms) < 2:
        #raise ValueError('Only or less than 1 atom and No distance')
        symmetry_prop = 'NO'
        return symmetry_prop
    else:
        # p_all=struct.get_scaled_positions()
        # struct.set_scaled_positions(p_all)

        for i, atom1 in enumerate(surfatoms):
            dis = []
            #print(atom1)
            for j, atom2 in enumerate(surfatoms):
                #print(atom2)
                if atom1 == atom2:
                    continue
                else:
                    if isinstance(poscar, str):
                        struct = read(poscar, format='vasp')
                    else:
                        struct = poscar
                    p_all = struct.get_scaled_positions()
                    struct.set_scaled_positions(p_all)

                    p1_scaled = struct.get_scaled_positions()[surfatoms[i]]
                    p2_scaled = struct.get_scaled_positions()[surfatoms[j]]
                    if p1_scaled[0] - p2_scaled[0] >= 0.50:
                        p_all[surfatoms[j]][0] = p2_scaled[0] + 1
                        struct.set_scaled_positions(p_all)
                    elif p1_scaled[0] - p2_scaled[0] <= -0.50:
                        p_all[surfatoms[i]][0] = p1_scaled[0] + 1
                        struct.set_scaled_positions(p_all)
                    if p1_scaled[1] - p2_scaled[1] >= 0.50:
                        p_all[surfatoms[j]][1] = p2_scaled[1] + 1
                        struct.set_scaled_positions(p_all)
                    elif p1_scaled[1] - p2_scaled[1] <= -0.50:
                        p_all[surfatoms[i]][1] = p1_scaled[1] + 1
                        struct.set_scaled_positions(p_all)
                    dis += [
                        get_distances(struct.get_positions()[surfatoms[i]],
                                      struct.get_positions()[surfatoms[j]])[1][0][0]
                    ]
            #print(dis)
            dis_matrix += [min(dis)]
        #print(min(dis_matrix),max(dis_matrix),np.mean(dis_matrix))
        #if max(dis_matrix)-min(dis_matrix) > tol:
        if abs(np.mean(dis_matrix) - min(dis_matrix)) > tol:
            symmetry_prop = 'NO'
            return symmetry_prop
        else:
            symmetry_prop = 'YES'
        return symmetry_prop


from catkit.gen import defaults
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii, chemical_symbols
import numpy as np
if __name__ == '__main__':
    #get_potcar(['opt/POSstart','opt/POSend'])
    #bind_type,atom_b_conn_symb=get_site('opt/POSend',['NH2','N'])
    #Sub_repeat('opt/descriptor-all')
    #poscar='./opt/CONTCAR-Ag111'
    #dis_symb_matrix,dis_matrix=get_distance_adatoms(poscar,tol=0.1)
    #print(dis_symb_matrix,dis_matrix)
    #poscar='CONTCAR'
    #adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsurfatoms,subsurfatoms_symb=distinguish_atom_binding(poscar,tol=0.05)
    #print(adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsurfatoms,subsurfatoms_symb)
    #bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom(poscar)
    #print(bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_symb)
    #print(get_symmetry_surfatoms(poscar,tol=0.2))
    ''' 
   Efile='ads_final_b1' 
   spec_ads,spec_ads_stable,dir_list_final, spec_ads_stable_surfa=get_site_stable(Efile,Ecut=-0.1)
   #print(spec_ads)
   print(spec_ads_stable, spec_ads_stable_surfa)
   #print(dir_list_final)
   '''
    ''' 
   poscar='CONTCAR'
   atom = atomic_numbers['N']
   atom_index,atom_neighs=get_atom_neigh(poscar,atom)
   print(atom_index,atom_neighs)
   '''

    #poscar='Pt_Au_100_b1_N_3/CONTCAR'
    #bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_sym = get_binding_adatom(poscar)
    #print(bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_sym)
    #dis_symb_matrix,dis_matrix=get_distance_adatoms('opt/CONTCAR1')
    #print(dis_symb_matrix,dis_matrix)
    #adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsurfatoms,subsurfatoms_symb=distinguish_atom_binding(poscar,tol_layer=0.01,tol=0.02,base_layer=4,atoms_layer=9)
    #print( adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsurfatoms,subsurfatoms_symb)

    #atoms=['Au','Ag','Pt','Pd','Rh','Ru','Ir','Cu','N','O','H','Zn','Fe','Co','Ni']
    atoms = [
        'Os', 'C', 'Cu', 'Zn', 'Al', 'Sc', 'Ti', 'V', 'Cr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Hf', 'Ta',
        'W', 'Re'
    ]
    names_atom, radii_atoms, mass_atoms, vdw_radii_atoms = Extract_atomic_info(atoms)
    print(f'{atoms}\n{radii_atoms}\n{vdw_radii_atoms}')
    ''' 
   layer_adatoms,layer_adatoms_symb=distinguish_atom_binding('opt/CONTCAR', tol=0.05,layer='adatom')
   print("adatom:%s%s" %(layer_adatoms,layer_adatoms_symb))
   layer_surfatoms,layer_surfatoms_symb=distinguish_atom_binding('opt/CONTCAR', tol=0.05,layer='surf_atom')
   print("surf_atom:%s%s" %(layer_surfatoms,layer_surfatoms_symb)) 
   #print("surf_atom:%s" %(layer_surfatoms_symb))
   layer_subsurfatoms,layer_subsurfatoms_symb=distinguish_atom_binding('opt/CONTCAR', tol=0.05,layer='subsurf_atom')
   print("subsurf_atom:%s%s" %(layer_subsurfatoms,layer_subsurfatoms_symb))
   '''
    '''
   #radii = defaults.get('radii')
   #radii_adatom=[radii[i.numbers] for i in layer_atoms_symb]
   #number_adatom=[atomic_numbers[i] for i in layer_atoms_symb]
   radii_adatom=[covalent_radii[i] for i in number_adatom]
   radii_mean = np.mean(radii_adatom)
   print(radii_mean)
   '''
    ''' 
   adatoms,adatoms_symb,surfatoms,surfatoms_symb,subsurfatoms,subsurfatoms_symb=distinguish_atom_binding('opt/CONTCAR-N',tol=0.01)  
   print("adatom:%s%s" %(adatoms,adatoms_symb))
   print("surf_atom:%s%s" %(surfatoms,surfatoms_symb))
   bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_sym = get_binding_adatom('opt/CONTCAR-N')
   print(bind_adatoms,bind_adatoms_symb,adspecie,bind_type_symb,bind_surfatoms,bind_surfatoms_sym)
   '''
    '''
   descriptor=[]
   #feature=['Enegativity','Valence_electron','Atomic_radius','Atomic_mass']
   feature=['Valence_electron','Atomic_radius']
   feature_value_surf=Construct_descriptor_info('Info/Element_Info',layer_surfatoms_symb,feature)
   feature_value_subsurf=Construct_descriptor_info('Info/Element_Info',layer_subsurfatoms_symb,feature)
   feature_value=np.hstack((feature_value_surf,feature_value_subsurf))
   print(feature_value)
   descriptor += [np.around(np.mean(feature_value,0),2)]
   print(descriptor)
   '''
    #bind_adatoms,bind_adatoms_symb,bind_type_symb,adspecie,bind_surfatoms,bind_surfatoms_symb=get_binding_adatom('opt/POSend')
    #print(bind_adatoms,bind_adatoms_symb,bind_type_symb,adspecie,bind_surfatoms,bind_surfatoms_symb)

    #bind_type,atom_b_conn_symb=get_site('opt/CONTCAR',['NH2','N'])
    #print(bind_type)

    #dis_symb_matrix,dis_matrix=get_distance_adatoms('opt/CONTCAR')
    #print(dis_symb_matrix,dis_matrix)
    #s='opt/CONTCAR'
    #if isinstance(s, str):
    #print("string")
    #else:
    #print("not string")
    #Sub_repeat('opt/descriptor-all')
