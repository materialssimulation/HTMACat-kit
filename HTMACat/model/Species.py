from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from HTMACat.catkit.gratoms import *
from HTMACat.catkit.gen import utils
from HTMACat.catkit.gen.utils.utilities import to_gratoms
from ase import Atoms
from ase.build import molecule
from HTMACat.catkit.gen.adsorption import Builder
from abc import abstractmethod, ABC
from ase.io import read
from ase import neighborlist
import copy


class ABS_Species(ABC):
    def __init__(self, form, formtype="sim", alias_name=None):
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
        return self.alias_name

    @abstractmethod
    def get_molecule(self, randomSeed=0) -> Gratoms:
        pass

    @classmethod
    def from_input_dict(cls, init_dict):
        return cls(init_dict["name"], init_dict["form"], init_dict["formtype"])

    @classmethod
    def from_input(cls, input_dict):
        species_dict = {}
        for key, value in input_dict.items():
            new_dict = {key: cls(form=value, alias_name=key)}
            species_dict.update(new_dict)
        return species_dict


class Sim_Species(ABS_Species):
    def __init__(self, form, formtype="sim", alias_name=None):
        super().__init__(form, formtype, alias_name)

    def get_molecule(self, randomSeed=0):
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
        # todo
        ads_use_charges = -1
        return ads_molecule,ads_use_charges


class File_Species(ABS_Species):
    def __init__(self, form, formtype="file", alias_name=None):
        super().__init__(form, formtype, alias_name)
        if "." in form:
            str_list = form.split(".")
        self.filetype = str_list[-1]

    def set_filetype(self, typename):
        self.filetype = typename

    def out_file_name(self):
        if "." in self.form:
            return self.form.split(".")[0]
        else:
            return self.form

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

    def get_molecule(self, randomSeed=0) -> Gratoms:
        atoms = self.atoms
        edges_list = self.edges_list
        ads_molecule = to_gratoms(atoms, edges=edges_list)
        # todo
        ads_use_charges = -1
        return ads_molecule,ads_use_charges


class Sml_Species(ABS_Species):
    def __init__(self, form, formtype="sml", alias_name=None):
        super().__init__(form, formtype, alias_name)

    def out_file_name(self):
        if self.alias_name == self.form:
            ads1 = self.get_formular()
            mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
            ads1 = rdMolDescriptors.CalcMolFormula(mole)
            return ads1
        else:
            return self.alias_name.split("(")[0]

    def get_molecule(self, randomSeed=0) -> Gratoms:
        ### Changed by ZhaojieWang, 20230829 (<>改进：需能处理离子键可连接的SMILES)
        ads1 = self.get_formular()
        if '.' in ads1:
            if utils.Check_treatable__HTMATver(ads1):
                mole = utils.Gen_conn_mole(ads1)
            else:
                print('[ERROR]: Untreatable SMILES:', ads1)
                return None
        else:
            mole = Chem.AddHs(Chem.MolFromSmiles(ads1))
        stat = AllChem.EmbedMolecule(mole, randomSeed=randomSeed)
        if stat == -1:
            print("[WARNING]: No 3D conformer of specie %s can be generated, using the 2D version instead! (could be unreasonable)" % ads1)
        #try:
        #    AllChem.MMFFOptimizeMolecule(mole) # no randomseed param
        #except:
        #    try:
        #        AllChem.UFFOptimizeMolecule(mole) # no randomseed param
        #    except:
        #        pass
        conf = mole.GetConformer()
        atomicnums_list = []
        coords_list = []
        atomiccharges_list = []
        for i in range(mole.GetNumAtoms()):
            atomicnums_list.append(mole.GetAtomWithIdx(i).GetAtomicNum())
            coords_list.append(tuple(conf.GetAtomPosition(i)))
            atomiccharges_list.append(mole.GetAtomWithIdx(i).GetFormalCharge())
        edges_list = []
        for b in mole.GetBonds():
            edges_list.append((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        _idxtmp = rdMolDescriptors.CalcMolFormula(mole).find('-')
        _idxtmp1 = rdMolDescriptors.CalcMolFormula(mole).find('+')
        # ASE处理不了带电的化学式，计算出的电荷量又必然跟在化学式末尾，故只需截取电荷符号之前的部分
        if -1 == _idxtmp and -1 == _idxtmp1:
            form_str = rdMolDescriptors.CalcMolFormula(mole)
        elif _idxtmp != -1 and _idxtmp1 == -1:
            form_str = rdMolDescriptors.CalcMolFormula(mole)[:_idxtmp]
        elif _idxtmp == -1 and _idxtmp1 != -1:
            form_str = rdMolDescriptors.CalcMolFormula(mole)[:_idxtmp1]
        else:
            raise ValueError('Invalid SMILES') # 不可能的分支
        atoms = Atoms(form_str, coords_list)
        atoms.set_atomic_numbers(atomicnums_list)
        ads_molecule = to_gratoms(atoms, edges=edges_list)
        return ads_molecule, atomiccharges_list


def species_from_input(init_dict):
    # ads_init_dict = {'SML':False,'type':[],'value':[]}
    species_dict = {}
    species_type_list = ["sim", "sml", "file"]
    for key, value in init_dict.items():
        if key == "sim":
            new_dict = Sim_Species.from_input(init_dict["sim"])
        elif key == "sml":
            new_dict = Sml_Species.from_input(init_dict["sml"])
        elif key == "file":
            new_dict = File_Species.from_input(init_dict["file"])
        else:
            msg = ",".join(species_type_list)
            warn_msg = (
                "Only support species type: %s, Your input %s part in Species will be dismiss"
                % (msg, key)
            )
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

    assert isinstance(init_str, dict), "The species in adsorption must be dict or str"

    for key, value in init_str.items():
        if key == "s" or key == "sml":
            return Sml_Species(value)
        elif key == "f" or key == "file":
            return File_Species(value)
        else:
            wrn_msg = ", ".join(["s", "sml", "f", "file"])
            raise Warning("The key of initial dict should be one of the following:\n %s" % wrn_msg)
