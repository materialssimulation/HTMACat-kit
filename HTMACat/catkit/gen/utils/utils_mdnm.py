### 本文件中是一些与主流程无关的功能函数，可能调用的位置不确定，请暂不要整合到任何类中(zjwang 20230519)
import copy
import numpy as np
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from rdkit import Chem
from ase.data import atomic_numbers, covalent_radii

def mol_to_graph(self, m):
    """Convert a molecule object to a graph.
    
    Parameters
    ----------
    m : mol
        The RDKit molecule object to be converted into a networkx graph.
    
    Returns
    -------
    G : Graph
        The networkx Graph object derived from m.
    """
    G = nx.Graph()
    for i_n in range(m.GetNumAtoms()):
        G.add_nodes_from([(i_n, {'number': m.GetAtomWithIdx(i_n).GetAtomicNum()})])
    bonds = [m.GetBondWithIdx(k) for k in range(len(m.GetBonds()))]
    edges = []
    for edge in bonds:
        edges.append((edge.GetBeginAtomIdx(), edge.GetEndAtomIdx()))
    G.add_edges_from(edges)
    return G

def solve_normal_vector_linearsvc(coords, bond_idx, take_into_account_idx=None):
    """Solve the adsorption direction of the given species using SVM. Returns the direction
    that should be rotated into [001] of the slab model.

    Parameters
    ----------
    coords : numpy.ndarray
        The coordinates of all atoms in the species.
    bond_idx : int
        The index of the atom to be placed on the adsorption sites.
    take_into_account_idx : list
        The list of atom indices should be as close as possible to the surface. If it is not
        None, the geometric center of these atoms will act as the symmetry center.

    Returns
    -------
    vec : list
        The normal vector of the decision boundary.
    flag_linearly_separable : list
        Record whether the extended coordinates are linearly separable.
    """
    if len(coords) == 1:
        return [0, 0, 1], True
    # confirm the symmetry center for generating mirror image atoms, and
    # label original atoms (except the bonding atom)
    if take_into_account_idx is None:
        coord_sym_center = coords[bond_idx]
        coords_ext = np.array(copy.deepcopy(coords))
        coords_ext = np.delete(coords_ext, bond_idx, 0)
        labels = np.array([1 for k in range(len(coords)-1)]) # 前面删掉了作为镜像中心的bonding atom
    else: # 以多位点几何中心为镜像中心，暂未启用
        tmp_vec = np.array([.0 for k in coords[0]])
        for i,idx in enumerate(take_into_account_idx):
            tmp_vec += np.array(coords[idx])
        coord_sym_center = tmp_vec / len(take_into_account_idx)
        print(coord_sym_center)
        coords_ext = np.array(copy.deepcopy(coords))
        labels = np.array([1 for k in range(len(coords))])
    # extend: mirror image atoms (labels are opposite)
    for i,coord in enumerate(coords_ext):
        coords_ext = np.concatenate((coords_ext, [2*coord_sym_center-coord]), axis=0)
    labels = np.concatenate((labels, 1-labels), axis=0)
    # find the normal vector of the SVC decision boundary
    svc = SVC(kernel='linear').fit(coords_ext, labels)
    k = np.sqrt(np.sum(svc.coef_[0]*svc.coef_[0]))
    vec = svc.coef_[0]/k
    flag_linearly_separable = (svc.predict(coords_ext) == labels).all()
    ### print(labels, svc.predict(coords_ext))
    return vec, flag_linearly_separable

def solve_normal_vector_pca(coords, bond_idx):
    """Solve the adsorption direction of the given species using PCA. Returns the direction
    that should be rotated into [001] of the slab model.
    3D to 3D PCA: [x, y, z] to [M1, M2, M3]

    Parameters
    ----------
    coords : numpy.ndarray
        The coordinates of all atoms in the species.
    bond_idx : int
        The index of the atom to be placed on the adsorption sites.

    Returns
    -------
    vec : list
        The shortest component M3.
    """
    # 1. 找方向
    if len(coords) == 1:
        vec = [0, 0, 1]
    else:
        X = np.array(coords)[:,:3]
        pca = PCA(n_components=3)
        pca.fit(X)
        vec = pca.components_[2]
    # 2. 定“上下”
    Center_coord = np.mean(coords, axis=0) # 形心坐标
    if np.dot(vec,(Center_coord-coords[bond_idx])) < 0:
        vec = -vec
    return vec # return M3

def solve_normal_vector_hetero(coords, symbols):
    centroid = np.mean(coords, axis=0)
    coords_hetero = []
    for i,coord in enumerate(coords):
        if not symbols[i] in ['C', 'H']:
            coords_hetero.append(coord)
    centroid_hetero = np.mean(coords_hetero, axis=0)
    vec = centroid - centroid_hetero
    return vec, centroid, centroid_hetero

def solve_principle_axe_pca(coords):
    """PCA: 2D to 1D, using x & y coords of atoms in an adsorbate
    """
    if len(coords) == 1:
        return [1, 1, 0]
    X = np.array(coords)[:,:2]
    pca = PCA(n_components=1)
    pca.fit(X)
    return pca.components_[0]

date = '20230905'

if __name__ == '__main__' and date == '20230510':
    coords_NH3 = 4 * np.array([ # NH3+.xyz
        [ 0.00000,  0.00000,  0.11649], # N
        [ 0.00000,  0.93973,  0.40800], # H
        [ 0.81383, -0.46986,  0.40808], # H
        [-0.81383, -0.46986,  0.40808], # H
    ])
    vec, flag = solve_normal_vector_linearsvc(coords_NH3, 0)
    print(vec, np.sum(vec*vec), flag)
    vec, flag = solve_normal_vector_linearsvc(coords_NH3, 0, [0,1])
    print(vec, np.sum(vec*vec), flag)
    vec, flag = solve_normal_vector_linearsvc(coords_NH3, 0, [1,2,3])
    print(vec, np.sum(vec*vec), flag)
    vec, flag = solve_normal_vector_linearsvc(coords_NH3, 0, [0,1,2,3])
    print(vec, np.sum(vec*vec), flag)

if __name__ == '__main__' and date == '20230905':
    coords_NH3 = 4 * np.array([ # NH3+.xyz
        [ 0.00000,  0.00000,  0.11649], # N
        [ 0.00000,  0.93973,  0.40800], # H
        [ 0.81383, -0.46986,  0.40808], # H
        [-0.81383, -0.46986,  0.40808], # H
    ])
    vec = solve_normal_vector_pca(coords_NH3, 0)
    print(vec, np.sum(vec*vec))
    coords_2 = 4 * np.array([
        [ 0.00000,  0.00000,  0.11649],
        [ 0.00000,  0.00000,  0.40800],
    ])
    vec = solve_normal_vector_pca(coords_NH3, 1)
    print(vec, np.sum(vec*vec))

if __name__ == '__main__' and date == '20230825':
    print(solve_principle_axe_pca(np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])))
    print(solve_principle_axe_pca(np.array([[-1,  0], [-2,  0], [-3, -1], [1, 2], [2, 2], [3, 3]])))
    print(solve_principle_axe_pca(np.array([[0,  0, 1], [2,  0, 3], [0, 1, 2], [2, 1, 5]])))

# 以下与分子生成相关，Species.py get_molecule()中使用
def Check_treatable__HTMATver(sml): # Zhaojie Wang 20230829   改动：兼容非中性物种
    # 若多段SMILES化合物中阴离子or阳离子只有一种，且各段都只有一个带电位点，该物质中的离子键是可以唯一正确连接的
    # 要求：1）必须所有段有且仅有一个带电位点；2）整个分子中只有一个正电位点或只有一个负电位点（作为中心与所有异号位点成离子键）
    smls = sml.split('.') # 按弱连接分段
    moles = [Chem.AddHs(Chem.MolFromSmiles(k)) for k in smls]
    treatable = True
    #
    poschg = [[] for k in range(len(moles))] # 各段正电荷原子下标
    negchg = [[] for k in range(len(moles))] # 各段负电荷原子下标
    for i,m in enumerate(moles):
        for atom in m.GetAtoms():
            if atom.GetFormalCharge() > 0:
                poschg[i].append(atom.GetIdx())
            elif atom.GetFormalCharge() < 0:
                negchg[i].append(atom.GetIdx())
            else:
                pass
    # 确定带正/负电的位点总数
    num_havepos = 0
    num_haveneg = 0
    for i in range(len(moles)):
        if (len(poschg[i]) + len(negchg[i]) != 1): # 条件1
            treatable = False
            break
        num_havepos += len(poschg[i])
        num_haveneg += len(negchg[i])
    if (num_havepos > 1) and (num_haveneg > 1): # 分子中正、负电位点都不止一个，不行
        treatable = False
    return treatable

def Gen_conn_mole(sml): # Zhaojie Wang 20230829   补充离子键使多段SMILES物种连通
    mole = Chem.AddHs(Chem.MolFromSmiles(sml))
    poschg = []
    negchg = []
    for atom in mole.GetAtoms():
        if atom.GetFormalCharge() > 0:
            poschg.append(atom.GetIdx())
        elif atom.GetFormalCharge() < 0:
            negchg.append(atom.GetIdx())
    if 1 == len(poschg):
        center = poschg[0]
        attach = negchg
    else:
        center = negchg[0]
        attach = poschg
    molew = Chem.RWMol(mole)
    btype = Chem.BondType.IONIC
    for idx in attach:
        molew.AddBond(center,idx,btype)
    return molew



def center_molecule(rotated_positions):
    # 计算x轴和y轴的最大值和最小值
    min_x = np.min(rotated_positions[:, 0])
    max_x = np.max(rotated_positions[:, 0])
    min_y = np.min(rotated_positions[:, 1])
    max_y = np.max(rotated_positions[:, 1])

    # 获取x轴最小值对应的y坐标和y轴最大值对应的x坐标
    min_x_y = rotated_positions[np.argmin(rotated_positions[:, 0])][1]
    max_y_x = rotated_positions[np.argmax(rotated_positions[:, 1])][0]
      # 计算原子的中心坐标
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    return center_x, center_y

def center_slab(final_positions):
    # 计算x轴和y轴的最大值和最小值
    min_x = np.min(final_positions[:, 0])
    max_x = np.max(final_positions[:, 0])
    min_y = np.min(final_positions[:, 1])
    max_y = np.max(final_positions[:, 1])

    # 获取x轴最小值对应的y坐标和y轴最大值对应的x坐标
    min_x_y = final_positions[np.argmin(final_positions[:, 0])][1]
    max_y_x = final_positions[np.argmax(final_positions[:, 1])][0]

    # 计算原子的中心坐标
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    return center_x, center_y

# 吸附构型评价相关

# Pauling 电负性数据
pauling_en = {
    'H': 2.20, #
    'Li': 0.98,
    'Be': 1.57,
    'B': 2.04,
    'C': 2.55,
    'N': 3.04,
    'O': 3.44,
    'F': 3.98, #
    'Na': 0.93,
    'Mg': 1.31,
    'Al': 1.61,
    'Si': 1.90,
    'P': 2.19,
    'S': 2.58,
    'Cl': 3.16, #
    'K': 0.82,
    'Ca': 1.00,
    'Sc': 1.36,
    'Ti': 1.54,
    'V': 1.63,
    'Cr': 1.66,
    'Mn': 1.55,
    'Fe': 1.83,
    'Co': 1.88,
    'Ni': 1.91,
    'Cu': 1.90,
    'Zn': 1.65,
    'Ga': 1.81,
    'Ge': 2.01,
    'As': 2.18,
    'Se': 2.55,
    'Br': 2.96,
    'Kr': 3.00,
    'Rb': 0.82, #
    'Sr': 0.95,
    'Y': 1.22,
    'Zr': 1.33,
    'Nb': 1.60,
    'Mo': 2.16,
    'Tc': 2.10,
    'Ru': 2.20,
    'Rh': 2.28,
    'Pd': 2.20,
    'Ag': 1.93,
    'Cd': 1.69,
    'In': 1.78,
    'Sn': 1.96,
    'Sb': 2.05,
    'Te': 2.10,
    'I': 2.66,
    'Xe': 2.60,
    'Cs': 0.79,
    'Ba': 0.89,
    'La': 1.10,
    'Ce': 1.12,
    'Pr': 1.13,
    'Nd': 1.14,
    'Pm': 1.13,
    'Sm': 1.17,
    'Eu': 1.20,
    'Gd': 1.20,
    'Tb': 1.20,
    'Dy': 1.22,
    'Ho': 1.23,
    'Er': 1.24,
    'Tm': 1.25,
    'Yb': 1.10,
    'Lu': 1.27,
    'Hf': 1.30,
    'Ta': 1.50,
    'W': 2.36,
    'Re': 2.20,
    'Os': 2.20,
    'Ir': 2.20,
    'Pt': 2.28,
    'Au': 2.54,
    'Hg': 2.00,
    'Tl': 1.82,
    'Pb': 2.33,
    'Bi': 2.02,
    'Po': 2.00,
    'At': 2.20,
    'Rn': 2.20,
    'Fr': 0.70,
    'Ra': 0.90,
    'Ac': 1.10,
    'Th': 1.30,
    'Pa': 1.50,
    'U': 1.38,
    'Np': 1.36,
    'Pu': 1.28,
    'Am': 1.30,
    'Cm': 1.30,
    'Bk': 1.30,
    'Cf': 1.30,
    'Es': 1.30,
    'Fm': 1.30,
    'Md': 1.30,
    'No': 1.30,
    'Lr': 1.30,
    'Rf': 1.50,
    'Db': 1.50,
    'Sg': 1.50,
    'Bh': 1.50,
    'Hs': 1.50,
    'Mt': 1.50,
    'Ds': 1.50,
    'Rg': 1.50,
    'Cn': 1.50,
    'Nh': 1.50,
    'Fl': 1.50,
    'Mc': 1.50,
    'Lv': 1.50,
    'Ts': 1.50,
    'Og': 1.50,
}

def fitness_01(distance, rsum, EN): # hetero <-> metallic_site
    if distance < 1e-4:
        d = 1e-4
    else:
        d = distance
    return EN * ( (2*d/rsum)**(-1) - (2*d/rsum)**-2 ) * 4

def score_configuration_hetero(coords, symbols, z_surf):
    # print('*** test', fitness_01(1.8, 1.8, 3.0)) # = 3
    # print('*** test', fitness_01(0.9, 1.8, 3.0)) # = 0
    # print('*** test', fitness_01(1.8*1.5, 1.8, 3.0)) # = 2.667
    idx_hetero = [] # 记录各个杂原子在coords中的index
    neigh_of_hetero = [] # 记录idx_hetero中对应的各个杂原子的‘邻近表面原子’在coords中的index列表
    neigh_dist_of_hetero = [] # 与neigh_of_hetero同形，记录对应的距离
    neigh_rsum_of_hetero = [] # 与neigh_of_hetero同形，记录对应的共价半径之和rsum
    for idx,coord in enumerate(coords):
        if (coord[2] > z_surf) and (not symbols[idx] in ['C', 'H']):
            idx_hetero.append(idx)
            tmp_neighl = []
            tmp_distl = []
            tmp_rsuml = []
            for idxn,coordn in enumerate(coords):
                if abs(coordn[2] - z_surf) < 1.0: # 识别表层原子
                    d = np.sqrt( np.sum( np.dot(coord-coordn,coord-coordn) ) )
                    rsum = covalent_radii[atomic_numbers[symbols[idx]]] + covalent_radii[atomic_numbers[symbols[idxn]]]
                    if d < 2.5*rsum:
                        tmp_neighl.append(idxn)
                        tmp_distl.append(d)
                        tmp_rsuml.append(rsum)
            neigh_of_hetero.append(tmp_neighl)
            neigh_dist_of_hetero.append(tmp_distl)
            neigh_rsum_of_hetero.append(tmp_rsuml)
    # calc score
    score = 0
    for idx_a,atom_i in enumerate(idx_hetero):
        for idx_n,atom_j in enumerate(neigh_of_hetero[idx_a]):
            score = score + fitness_01(distance=neigh_dist_of_hetero[idx_a][idx_n],
                                       rsum=neigh_rsum_of_hetero[idx_a][idx_n],
                                       EN=pauling_en[symbols[atom_i]])
    # 检查是否有与表面原子距离过近的H（比H、?共价半径之和小0.1A以上），如果有，将此构型的整体score乘以一个小于1的惩罚系数（每有一对过近的H-?则产生一项惩罚系数）
    for idxh,coordh in enumerate(coords):
        if (coordh[2] > z_surf) and (symbols[idxh] == 'H'):
            for idxsurf,coordsurf in enumerate(coords):
                if abs(coordsurf[2] - z_surf) < 1.0: # 识别表层原子
                    d = np.sqrt( np.sum( np.dot(coordh-coordsurf,coordh-coordsurf) ) )
                    rsum = covalent_radii[atomic_numbers[symbols[idxh]]] + covalent_radii[atomic_numbers[symbols[idxsurf]]]
                    if d < rsum - 0.1: # 距离太近，触发惩罚项
                        k = d / rsum
                        score = score * k
    return score
