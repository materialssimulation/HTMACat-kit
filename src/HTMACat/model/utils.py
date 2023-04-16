import numpy as np
import networkx as nx

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

def find_topsurface_atoms(coords, tol_zdiff=0.7, tol_zangle_min=0):
    '''
    generate the list of surface atoms (top surface)
    params:
        coords: aoordinates of all atoms, numpy.ndarray
        tol_zdiff: if the z_coord of an atom is higher than zmax-tol_zdiff,
                    this atom is recognized as a "surface atom"
    returns:
        index_topsurf: list of indices of surface atoms
    '''
    index_topsurf = []
    zmax = np.max(coords[:,2])
    for i,icoord in enumerate(coords):
        if icoord[2] > zmax - tol_zdiff:
            index_topsurf.append(i)
    return index_topsurf
