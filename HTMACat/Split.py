#lbxent
import numpy as np
from ase import Atoms
from HTMACat.catkit.gen import utils
from collections import defaultdict
from collections import Counter,OrderedDict
from ase.data import atomic_numbers, atomic_names, covalent_radii





def coads_split(filename,key_atom):
    print('Dealing with vasp file:', filename,key_atom)
    contcar_file_path = filename
    with open(contcar_file_path, 'r') as f:
        contcar_content = f.readlines()
    lattice_matrix = np.array([list(map(float, line.split())) for line in contcar_content[2:5]])
    coords_start_line = 9 
    atomic_coords = [line.split()[:3] for line in contcar_content[coords_start_line:]]
    atomic_coords = [[float(coord) for coord in coords] for coords in atomic_coords]
    threshold_z = 0.1 #需要改动（改为真实值？）

    grouped_coords = defaultdict(list)
    
    # 将 Z 轴坐标按照有效数字进行分组
    for coord in atomic_coords:
        z_value = round(coord[2], 2)  # 取两位有效数字
        grouped_coords[z_value].append(coord)
    threshold_z_values = [max(coords, key=lambda c: c[2]) for z_value, coords in grouped_coords.items() if len(coords) >= 6]
    threshold_z0 = max(threshold_z_values, key=lambda c: c[2])[2]
   


    substrate_atoms = [coord for coord in atomic_coords if coord[2] < threshold_z + threshold_z0]
    molecule_atoms = [coord for coord in atomic_coords if coord[2] >= threshold_z + threshold_z0]
    molecule_atoms_indices = [index for index, coord in enumerate(atomic_coords) if coord in molecule_atoms]
    #print(molecule_atoms_indices)



    cartesian_coordinates = np.dot(molecule_atoms, lattice_matrix)
    cartesian_coordinates_list = cartesian_coordinates.tolist()
    
    


    #原子数及元素符号
    element_symbols = contcar_content[5].split()
    atom_counts = [int(count) for count in contcar_content[6].split()]
    all_atoms = [(element, count) for element, count in zip(element_symbols, atom_counts)]
    all_atoms_list = [char * count if len(char) == 1 else [char] * count for char , count in all_atoms]
    all_atoms_list0 = [char for sublist in all_atoms_list for char in sublist]#所有元素包含基底与分子
    #print(all_atoms_list0)
    molecule_elements = [all_atoms_list0[index] for index in molecule_atoms_indices]
    unique_atoms = list(dict.fromkeys(molecule_elements))#获得分子元素
    #print(molecule_elements,unique_atoms)
    atom_counts1 = Counter(molecule_elements)
    atom_counts1_dict = dict(atom_counts1)
    repeat_counts = list(atom_counts1_dict.values())#每个元素对应的数量
    #print(repeat_counts)


    selected_elements = unique_atoms #输入
    #print(selected_elements)
    number_list = [atomic_numbers[key] for key in selected_elements]
    r_list = [covalent_radii[k] for k in number_list]
    

    #total_atoms = sum(count for element, count in zip(element_symbols, atom_counts) if element in selected_elements)
    #molecule_atoms_copy = list(molecule_atoms)
    selected_molecule_elements = [(element, count) for element, count in zip(element_symbols, atom_counts) if element in selected_elements]





    unzipped = list(zip(*selected_molecule_elements))
    #counts_list = list(unzipped[1])
    counts_list = repeat_counts
    #print(counts_list,number_list)
    number_counts = list(zip(number_list,counts_list))
    #print(number_counts)


    target_key = atomic_numbers[key_atom]
    #number_counts = [(number, 1) if number == target_key else (number, count) for number, count in number_counts] #基底含中心原子元素
    #sum_H = len(molecule_atoms) - sum(value for key, value in number_counts if key != 1) #基底被羟基化
    #index_of_key_1 = next((index for index, (key, _) in enumerate(number_counts) if key == 1), None)
    #number_counts[index_of_key_1] = (1,sum_H)
    #print(number_counts)

    atom_coordinates_list = []

    for element, count in number_counts:
        for _ in range(count):
            atom_coords = cartesian_coordinates_list.pop(0)
            atom_coordinates_list.append((element, atom_coords))

    result_rows = [index for index , row in enumerate(atom_coordinates_list) if row[0] == atomic_numbers[key_atom]]
    molecule_coords = np.array(cartesian_coordinates)
    
    
    for i in range(len(molecule_coords)):
        for j in range(i+1, len(molecule_coords)):
            distance = np.linalg.norm(molecule_coords[i] - molecule_coords[j])
    
    from scipy.spatial.distance import pdist, squareform

    threshold_factor = 1.3
    distances = squareform(pdist(np.array([coord[1] for coord in atom_coordinates_list])))
    bond_matrix = np.zeros_like(distances, dtype=int)
    for i in range(len(atom_coordinates_list)):
        for j in range(i + 1, len(atom_coordinates_list)):
            r1 = covalent_radii[atom_coordinates_list[i][0]]
            r2 = covalent_radii[atom_coordinates_list[j][0]]
            threshold_distance = threshold_factor * (r1 + r2)

            if distances[i, j] < threshold_distance:
                bond_matrix[i, j] = bond_matrix[j, i] = 1
    #函数调用
    def find_columns_0(matrix, selected_rows):
        return {row: np.where(matrix[row] == 1)[0] for row in selected_rows}
    def find_columns_1(matrix, selected_rows, excluded_column):
        return {row: np.where(np.logical_and(matrix[row] == 1, np.arange(len(matrix[row])) != excluded_column))[0] for row in selected_rows}

    result0 = result_rows
    result = find_columns_0(bond_matrix, result0)
    key = list(result.keys())
    for i in range(len(key)):
        result_list = result[key[i]].tolist()
    result1 = find_columns_1(bond_matrix,result_list,result0)
    atoms_with = result1#
    key0 = list(result1.keys())

    #print(bond_matrix)



    for i in range(len(key0)):
        result_list1 = result1[key0[i]].tolist()
        result22 = find_columns_1(bond_matrix, result_list1, result_list[i])
        key1 = list(result22.keys())
        for x in range(len(key1)):
            result_list2 = result22[key1[x]].tolist()
            atoms_with[key0[i]] = np.append(atoms_with[key0[i]],result_list2)
            result33 = find_columns_1(bond_matrix, result_list2, result_list1[x])
            key2 = list(result33.keys())
            for z in range(len(key2)):
                result_list3 = result33[key2[z]].tolist()
                atoms_with[key0[i]] = np.append(atoms_with[key0[i]], result_list3)
                result44 = find_columns_1(bond_matrix, result_list3, result_list2[z])
                key3 = list(result44.keys())
                for y in range(len(key3)):
                    result_list4 = result44[key3[y]].tolist()
                    atoms_with[key0[i]] = np.append(atoms_with[key0[i]], result_list4)
                    result55 = find_columns_1(bond_matrix, result_list4, result_list3[z])
                    key4 = list(result55.keys())
                    
    final_list = []
    for key, values in atoms_with.items():
        final_list.append([key]+values.tolist())
    final_list = [tuple(item) for item in final_list]
    #print(final_list)

    from ase import Atoms

    def move_point_away_xy(a, b, lattice_matrix, distance=3):
        a_xy = np.array(a[:2])
        b_xy = np.array(b[:2])
        a_actual_xy = np.dot(a_xy, lattice_matrix[:2, :2])
        b_actual_xy = np.dot(b_xy, lattice_matrix[:2, :2])
        ab_vector_xy = b_actual_xy - a_actual_xy
        ab_length_xy = np.linalg.norm(ab_vector_xy)
        unit_vector = ab_vector_xy / ab_length_xy
        new_b_actual_xy = b_actual_xy + distance * unit_vector
        new_b_original_xy = np.linalg.solve(lattice_matrix[:2, :2].T, new_b_actual_xy)
        diff_xy = new_b_original_xy - b_xy[:2]

        return diff_xy
    

    #分离操作
    def separate_rows(molecule_coords, selected_rows):
        selected_list = []
        other_list = []

        for i, atom_coord in enumerate(molecule_coords):
            if i in selected_rows:
                selected_list.append(atom_coord)
            else:
                other_list.append(atom_coord)

        return np.array(selected_list), np.array(other_list)




    atom_key = result0[0]
    for tup in final_list:#各个基团对应的原始位置，如返回contcar需要+9,每一步重新读取，生成不同的结构
        #处理行，转为真实坐标
        contcar_file_path = filename
        with open(contcar_file_path, 'r') as f:
            contcar_content = f.readlines()
        lattice_matrix = np.array([list(map(float, line.split())) for line in contcar_content[2:5]])
        coords_start_line = 9  # 如果从第十行开始
        delta_xy = move_point_away_xy(molecule_atoms[atom_key], molecule_atoms[tup[0]], lattice_matrix, distance=3.8)
        #print(delta_xy)
        atomic_coords = [list(map(float, line.split()[:3])) for line in contcar_content[coords_start_line:]]
        atomic_coords1 = atomic_coords
        np.set_printoptions(precision=16)
        atomic_coords = np.array(atomic_coords)

        #new_contcar_file_path1 = f"origin_CONTCAR"
        ## 将原子坐标写入文件
        #with open(new_contcar_file_path1, 'w') as f:
        #    f.writelines(contcar_content[:coords_start_line])
        #    for coord in atomic_coords1:
        #        f.write(f"  {coord[0]:.16f}  {coord[1]:.16f}  {coord[2]:.16f}   T   T   T\n")
#
        #print("原子坐标文件已生成", new_contcar_file_path1)

        #流程1分离
        selected_list, other_list = separate_rows(molecule_atoms, tup)
        first_atom = tup[0]
        
 
        
        ############################################################################################################################
        #ase处理，旋转以及移动

        v01 = np.array(cartesian_coordinates[first_atom])-np.array(cartesian_coordinates[-1])
        

        
        
        other_list2 = np.dot(other_list, lattice_matrix)
        other_list2 = other_list2 - np.array(cartesian_coordinates[atom_key])    #处理到原点
        selected_list2 = np.dot(selected_list, lattice_matrix)
        selected_list2 = selected_list2 - np.array(cartesian_coordinates[first_atom]) 
        if len(selected_list) >= 1 :#

            #中心分子部分旋转
            symbols = ['H'] * len(other_list2)
            ase_atoms1 = Atoms(symbols=symbols, positions=other_list2)
            ase_atoms_test = Atoms('H', positions=[(v01[0],v01[1],v01[2])])
            x_0 = np.degrees(np.arctan2(v01[1],v01[2]))
            ase_atoms1.rotate(x_0,'x')
            ase_atoms_test.rotate(x_0,'x') #相对原子
            atom_rotation2 = ase_atoms_test.positions[0]
            y_0 = np.degrees(-np.arctan2(atom_rotation2[0],atom_rotation2[2]))
            ase_atoms1.rotate(y_0+180,'y')
            
            #基团旋转
            symbols = ['H'] * len(selected_list2)
            ase_atoms2 = Atoms(symbols=symbols, positions=selected_list2)
            ase_atoms_test1 = Atoms('H', positions=[(-v01[0],-v01[1],-v01[2])])
            x_1 = np.degrees(np.arctan2(-v01[1],-v01[2]))
            ase_atoms_test1.rotate(x_1,'x')
            ase_atoms2.rotate(x_1,'x')
            atom_rotation3 = ase_atoms_test1.positions[0]
            y_1 = np.degrees(-np.arctan2(atom_rotation3[0],atom_rotation3[2]))
            ase_atoms2.rotate(y_1+180,'y')

            #处理到基底上方一定距离处
            substrate_atoms0 = np.dot(substrate_atoms, lattice_matrix)
            z_0 = substrate_atoms0[:, 2]
            max_z = np.max(z_0)
            z_1 = (ase_atoms1.positions+cartesian_coordinates[atom_key])[:, 2]
            z_2 = (ase_atoms2.positions+cartesian_coordinates[first_atom])[:, 2]
            min_z1 = np.min(z_1)
            min_z2 = np.min(z_2)
            dez1 = max_z - min_z1 + 2.0  #2埃处
            dez2 = max_z - min_z2 + 2.0
            dez1 = np.array([dez1]).reshape(1,-1)
            dez2 = np.array([dez2]).reshape(1,-1)
            ase_atoms1.positions[:, 2] = ase_atoms1.positions[:, 2] + dez1
            ase_atoms2.positions[:, 2] = ase_atoms2.positions[:, 2] + dez2


            inverse_lattice_matrix = np.linalg.inv(lattice_matrix.T)
            other_list3 = np.dot(ase_atoms1.positions, inverse_lattice_matrix)
            selected_list3 = np.dot(ase_atoms2.positions, inverse_lattice_matrix)

            other_list3 = other_list3 + np.array(molecule_atoms[atom_key])
            selected_list3 = selected_list3 + np.array(molecule_atoms[first_atom])

        #else:
        #    other_list3 = other_list
        #    selected_list3 = selected_list

        #合并操作
        #if len(other_list) <= 5:
        #    other_list3 = other_list
        #    selected_list = selected_list3
        
        restored_result = np.vstack([selected_list, other_list])
        num_rows = restored_result.shape[0]
        a = 0
        int_tuple = tuple(int(x) for x in tup)
        sorted_tup = tuple(sorted(int_tuple))
        #print(int_tuple)
        for index in sorted_tup:
            restored_result[index] = selected_list3[a]
            a = a+1
        z = 0
        y = 0
        for i in range(num_rows):
            x = 0
            for t in range(0, len(sorted_tup)):
                if i == sorted_tup[t]:
                    x = 1
            if x == 0:
                restored_result[i] = other_list3[y]
                y = y + 1
            i = i + 1
    

 
        for p1, p2 in enumerate(molecule_atoms):
            index_in_original_list = np.where((atomic_coords[:, :3] == p2).all(axis=1))[0][0]
            atomic_coords[index_in_original_list] = restored_result[p1]

        atomic_coords = atomic_coords.tolist()
        restored_result = restored_result.tolist()

        new_contcar_file_path0 = f"{tup[0]}_CONTCAR"
        # 将原子坐标写入文件
        with open(new_contcar_file_path0, 'w') as f:
            f.writelines(contcar_content[:coords_start_line])
            for coord in atomic_coords:
                f.write(f"  {coord[0]:.16f}  {coord[1]:.16f}  {coord[2]:.16f}   T   T   T\n")

        print("原子坐标过渡文件已生成", new_contcar_file_path0)
       
        new_contcar_file_path = f"CONTCAR_{tup}"
        #print(new_contcar_file_path0)
        with open(f"{tup[0]}_CONTCAR", 'r') as f:
            contcar_content1 = f.readlines()
        atomic_coords00 = [list(map(float, line.split()[:3])) for line in contcar_content1[coords_start_line:]]

       

        for i,elem in enumerate(tup):
            molecule_atom_to_find = restored_result[int(elem)]#找到基团位置
            for index, coords in enumerate(atomic_coords00):
            # 将坐标值乘以1000，然后取整数部分
                rounded_coords = [round(coord * 1000) for coord in coords]
                if rounded_coords == [round(value * 1000) for value in molecule_atom_to_find]:
                    index_in_atomic_coords = index
                    break
            #index_in_atomic_coords = atomic_coords00.index(molecule_atom_to_find)
            #移动坐标同时输出文件，这一步是改变坐标，需要改进,根据中心坐标以及基团中心坐标移动，z轴再移动
            modified_coord = np.array(atomic_coords00[index_in_atomic_coords])
            modified_coord[:2] += delta_xy
            # 替换特定行的坐标
            contcar_content1[coords_start_line + index_in_atomic_coords] = f"  {modified_coord[0]:.16f}  {modified_coord[1]:.16f}  {modified_coord[2]:.16f}   T   T   T\n"
            # 将新内容写入新的CONTCAR文件
        with open(new_contcar_file_path, 'w') as f:
            f.writelines(contcar_content1)
        print("基团分离后的CONTCAR文件已生成:", new_contcar_file_path)

