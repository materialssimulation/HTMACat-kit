### 本文件中是一些与主流程无关的功能函数，可能调用的位置不确定，请暂不要整合到任何类中(zjwang 20230519)
import copy
import numpy as np
from sklearn.svm import SVC

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



if __name__ == '__main__':
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
