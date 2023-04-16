import numpy as np
import networkx as nx

def MolToNXGraph(m): # (Last modified: 20230416, zjwang)
    """
    Convert a molecule object to a graph.
    Parameters
    ----------
    m : mol
        The RDKit molecule object to be converted into a networkx graph.
    Returns
    ----------
    G : Graph
        The networkx Graph object derived from m.
    """
    G = nx.Graph()
    for i_n in range(m.GetNumAtoms()):
        G.add_nodes_from([(i_n, {'number':m.GetAtomWithIdx(i_n).GetAtomicNum()})])
    bonds = [m.GetBondWithIdx(k) for k in range(len(m.GetBonds()))]
    edges = []
    for edge in bonds:
        edges.append((edge.GetBeginAtomIdx(),edge.GetEndAtomIdx()))
    G.add_edges_from(edges)
    return G

def find_topsurface_atoms(coords, tol_zdiff=0.7, tol_zangle_min=0): # (Last modified: 20230416, zjwang)
    """
    Generate the list of surface atoms (top surface).
    Parameters
    ----------
    coords: array_like
        Aoordinates of all atoms.
    tol_zdiff: number
        If the z_coord of an atom is higher than zmax-tol_zdiff, this atom is recognized as a "surface atom".
    Returns
    ----------
    index_topsurf: list
        List of the indices of the surface atoms.
    """
    index_topsurf = []
    zmax = np.max(coords[:,2])
    for i,icoord in enumerate(coords):
        if icoord[2] > zmax - tol_zdiff:
            index_topsurf.append(i)
    return index_topsurf
