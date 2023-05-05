from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from HTMACat.catkit.gratoms import *
from HTMACat.catkit.gen.utils.utilities import to_gratoms
from ase import Atoms
from ase.build import molecule
from HTMACat.catkit.gen.adsorption import Builder
from abc import abstractmethod, ABC
from ase.io import read
from ase import neighborlist
import copy


class ABS_Species(ABC):
    def __init__(self, form, formtype='sim', alias_name=None):
        self.__type = formtype
        self.form = form.strip()
        if alias_name is None:
            self.alias_name = form
        else:
            self.alias_name = alias_name

    def get_formular(self):
        return self.form

    def out_print(self):
        return self.alias_name

    def out_file_name(self):
        return self.get_formular()

    @abstractmethod
    def get_molecule(self) -> Gratoms:
        pass

    @classmethod
    def from_input_dict(cls, init_dict):
        return cls(init_dict['name'], init_dict['form'], init_dict['formtype'])

    @classmethod
    def from_input(cls, input_dict):
        species_dict = {}
        for key, value in input_dict.items():
            new_dict = {key: cls(form=value, alias_name=key)}
            species_dict.update(new_dict)
        return species_dict


class Sim_Species(ABS_Species):
    def __init__(self, form, formtype='sim', alias_name=None):
        super().__init__(form, formtype, alias_name)

    def get_molecule(self):
        ads1 = self.get_formular()
        atoms = molecule(ads1)
        cutOff = neighborlist.natural_cutoffs(atoms)
        neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
        neighborList.update(atoms)
        matrix = neighborList.get_connectivity_matrix()
        edges_list = []
        for i in range(matrix.shape[0]):
            for j in range(i):
                if matrix[i, j] == 1:
                    edges_list.append((i, j))
        ads_molecule = to_gratoms(atoms, edges=edges_list)
        return ads_molecule


class File_Species(ABS_Species):
    def __init__(self, form, formtype='file', alias_name=None):
        super().__init__(form, formtype, alias_name)
        if '.' in form:
            str_list = form.split('.')
        self.filetype = str_list[-1]

    def set_filetype(self, typename):
        self.filetype = typename

    def out_file_name(self):
        return self.alias_name

    @property
    def atoms(self) -> Atoms:
        atoms = read(self.get_formular(), format=self.filetype)
        return atoms

    @property
    def edges_list(self):
        atoms = self.atoms
        cutOff = neighborlist.natural_cutoffs(atoms)
        neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
        neighborList.update(atoms)
        matrix = neighborList.get_connectivity_matrix()
        edges_list = []
        for i in range(matrix.shape[0]):
            for j in range(i):
                if matrix[i, j] == 1:
                    edges_list.append((i, j))
        return edges_list

    def get_molecule(self) -> Gratoms:
        atoms = self.atoms
        edges_list = self.edges_list
        ads_molecule = to_gratoms(atoms, edges=edges_list)
        return ads_molecule


class Sml_Species(ABS_Species):
    def __init__(self, form, formtype='sml', alias_name=None):
        super().__init__(form, formtype, alias_name)

    def out_file_name(self):
        ads1 = self.get_formular()
        mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
        ads1 = rdMolDescriptors.CalcMolFormula(mole)
        return ads1

    def MolToNXGraph(self, m):
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
            G.add_nodes_from([(i_n, {'number': m.GetAtomWithIdx(i_n).GetAtomicNum()})])
        bonds = [m.GetBondWithIdx(k) for k in range(len(m.GetBonds()))]
        edges = []
        for edge in bonds:
            edges.append((edge.GetBeginAtomIdx(), edge.GetEndAtomIdx()))
        G.add_edges_from(edges)
        return G

    def get_molecule(self) -> Gratoms:
        '''
        mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
        G = self.MolToNXGraph(mole)
        ads1_list = molecule(rdMolDescriptors.CalcMolFormula(mole))
        ads_molecule = ads1_list[0]
        for ads_ in ads1_list:
            nm = iso.categorical_node_match('number', 6)
            if nx.is_isomorphic(ads_._graph, G, node_match=nm):
                ads_molecule = ads_
                print(ads_molecule._graph)
                break
        '''
        ads1 = self.get_formular()
        mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
        stat = AllChem.EmbedMolecule(mole)
        if stat == -1:
            print(
                '[WARNING]: No 3D conformer of specie %s can be generated, using the 2D version instead! (could be unreasonable)' % ads1)
        conf = mole.GetConformer()
        atomicnums_list = []
        coords_list = []
        for i in range(mole.GetNumAtoms()):
            atomicnums_list.append(mole.GetAtomWithIdx(i).GetAtomicNum())
            coords_list.append(tuple(conf.GetAtomPosition(i)))
        edges_list = []
        for b in mole.GetBonds():
            edges_list.append((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        atoms = Atoms(rdMolDescriptors.CalcMolFormula(mole), coords_list)
        atoms.set_atomic_numbers(atomicnums_list)
        ads_molecule = to_gratoms(atoms, edges=edges_list)
        return ads_molecule


def species_from_input(init_dict):
    # ads_init_dict = {'SML':False,'type':[],'value':[]}
    species_dict = {}
    species_type_list = ['sim', 'sml', 'file']
    for key, value in init_dict.items():
        if key == 'sim':
            new_dict = Sim_Species.from_input(init_dict['sim'])
        elif key == 'sml':
            new_dict = Sml_Species.from_input(init_dict['sml'])
        elif key == 'file':
            new_dict = File_Species.from_input(init_dict['file'])
        else:
            msg = ','.join(species_type_list)
            warn_msg = 'Only support species type: %s, Your input %s part in Species will be dismiss' % (msg, key)
            raise Warning(warn_msg)
        species_dict.update(new_dict)

    return species_dict


def init_from_ads(init_str, species_dict=None):
    if species_dict is None:
        species_dict = {}
    if isinstance(init_str, str):
        if init_str in species_dict:
            return copy.deepcopy(species_dict[init_str])
        else:
            return Sim_Species(init_str)

    assert isinstance(init_str, dict), 'The species in adsorption must be dict or str'

    for key, value in init_str.items():
        if key == 's' or key == 'sml':
            return Sml_Species(value)
        elif key == 'f' or key == 'file':
            return File_Species(value)
        else:
            wrn_msg = ', '.join(['s', 'sml', 'f', 'file'])
            raise Warning('The key of initial dict should be one of the following:\n %s' % wrn_msg)
