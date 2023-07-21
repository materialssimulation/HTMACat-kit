from numpy import array
from rdkit import Chem


# ₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ ₊ ₋ ₌ ₍ ₎ ₐ ₑ ₒ ₓ ₔ ₕ ₖ ₗ ₘ ₙ ₚ ₛ |⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ ⁼ ˂ ˃ ⁽ ⁾ ˙ * ′ ˙ ⁿ º
class Species(object):
    smiles_struc = {"O=C=O": "CO₂",
                    "[H][H]": "H₂",
                    "O=CO": "HCOOH",
                    "CO": "CH₃OH",
                    "C": "CH₄",
                    "[H+]": "H",
                    "O=C[O-]": "HCOO",
                    "[OH-]": "HO⁻",
                    "[CH3-]": "CH₃",
                    "C[O-]": "CH₃O",
                    "[C-4]": "C",
                    "[O-2]": "O",
                    "[C-]#[O+]": "CO",
                    "O=[C-]O": "COOH",
                    "[O-]CO": "H₂COOH",
                    "O=O": "O₂",
                    "[C-2]=[O+][O-]": "COO",
                    "[CH-3]": "CH",
                    "[C-2]=[OH+]": "COH",
                    "[O-]O": "HO₂",
                    "[CH-]=[O+][O-]": "CH-OO",
                    "[C-2]=[O+]O": "C-OOH",
                    "[CH-]=O": "CHO",
                    "[O-]C[O-]": "H₂COO",
                    "[CH2-2]": "CH₂",
                    "[CH-]=[OH+]": "CHOH",
                    "C=O": "CH₂O",
                    "OO": "H₂O₂",
                    "[CH-]=[O+]O": "CH-OO",
                    "[CH2-]O[O-]": "CH₂-OO",
                    "[CH2-]O": "CH₂OH",
                    "[CH2-]OO": "CH₂-OOH",
                    "CO[O-]": "CH₃-OO",
                    }

    def __init__(self, index=0, nature="(reactants)", chem_form="CO2", smiles="O=C=O", reduct_degree=None,
                 struc_form=None):
        self.nature = nature
        self.index = index
        self.chem_form = chem_form
        self.smiles = smiles
        self.reduct_degree = reduct_degree
        self.struc_form = struc_form
        if self.reduct_degree is None:
            self.set_reduct_degree()
        if self.struc_form is None and self.smiles in self.smiles_struc:
            self.struc_form = self.smiles_struc[self.smiles]
        else:
            self.struc_form = self.chem_form

    def set_reduct_degree(self):

        mol = Chem.MolFromSmiles(self.smiles)  # RDkit读入SMILE表达式初始化mol对象
        mol1 = Chem.AddHs(mol)  # 补充氢原子
        atoms = mol1.GetAtoms()  # 获取mol1中原子对象列表
        # print("---------")
        # num_carbon = 0
        reduct_degree = 0
        for atom in atoms:
            # print(f"原子列表：{atom.GetSymbol()}")
            if atom.GetSymbol() == "C":
                num_carbon = +1
                # self.reduct_degree = atom.GetTotalValence()
            if atom.GetSymbol() == "O":
                reduct_degree += 2
            if atom.GetSymbol() == "H":
                reduct_degree += -1
        self.reduct_degree = reduct_degree

    @classmethod
    def initialize(cls, list):
        # species_info, reactant, product = from_stem(file)
        return cls(list[0], list[1], list[2], list[3])  ## index,type,chemical formular, smiles


def from_stem(file):
    species_info, reactant, product = [], [], []
    stem = open(file, "r")
    mode = 1
    for i, line in enumerate(stem):
        if mode == 1:
            species_info_particles = line.split()
            # print(species_info_particles)
            if species_info_particles[0] == "--":
                if "reactions" in species_info_particles:
                    mode = 2
                continue
            species_info_particles[0] = int(species_info_particles[0])
            if len(species_info_particles) == 3:
                species_info_particles.insert(1, "(intermediate)")
            species_info.append(species_info_particles)
        if mode == 2:
            reactant_particles = []
            product_particles = []
            reaction_info = line.split()
            reaction_info.pop(0)
            if "<-->" in reaction_info:
                ind = reaction_info.index("<-->")
                for i in range(1, ind + 1, 2):
                    # print(reaction_info[ind - i])
                    reactant_particles.append(reaction_info[ind - i])
                for i in range(1, len(reaction_info) - ind, 2):
                    if "+" == reaction_info[ind + i] or "<-->" == reaction_info[ind + i]:
                        break
                    product_particles.append(reaction_info[ind + i])
            # print(f"reactant:{reactant_particles}")
            # print(f"product:{product_particles}")
            # print(reaction_info)
            reactant.append(reactant_particles)
            product.append(product_particles)
    return species_info, reactant, product


def get_reduct_degree_dict(species_info):
    species_list = []
    reduct_degree_dict = {}
    for i in range(0, len(species_info)):
        species_list.append(Species.initialize(species_info[i]))
        # print(species_list[i].reduct_degree)
        # {reduct_degree:num}
        if species_list[i].reduct_degree not in reduct_degree_dict:
            reduct_degree_dict[species_list[i].reduct_degree] = 1
        else:
            reduct_degree_dict[species_list[i].reduct_degree] += 1
    return species_list, reduct_degree_dict


def get_rgb():
    # RGB颜色
    rgb = array([(255, 0, 0),
                 (0, 255, 0),
                 (0, 0, 255),
                 (255, 255, 0),
                 (255, 0, 255),
                 (0, 255, 255),
                 (0, 0, 0),
                 (128, 0, 0),
                 (0, 128, 0),
                 (0, 0, 128),
                 (128, 128, 0),
                 (128, 0, 128),
                 (0, 128, 128),
                 (128, 128, 128),
                 (255, 128, 0),
                 (128, 255, 0),
                 (0, 255, 128),
                 (255, 0, 128),
                 (128, 0, 255),
                 (255, 128, 128),
                 (128, 128, 255),
                 (128, 255, 128),
                 (255, 255, 128),
                 (255, 128, 255),
                 (128, 255, 255),
                 (255, 64, 64),
                 (64, 255, 64),
                 (64, 64, 255),
                 (255, 255, 64),
                 (255, 64, 255),
                 (64, 255, 255),
                 (64, 64, 64)])
    return rgb / 255


def set_G(species_list, reduct_degree_dict):
    G = nx.Graph()
    reduct_degree_list = list(reduct_degree_dict.keys())
    reduct_degree_list.sort()
    # set the location of nodes
    for i, species in enumerate(species_list):
        G.add_node(species.index,
                   pos=(reduct_degree_dict[species.reduct_degree],
                        reduct_degree_list.index(species.reduct_degree)),
                   label=species.struc_form)
        reduct_degree_dict[species.reduct_degree] -= 2  # 为什么要减3？
    # set the edges
    combined_list = []
    for i in range(0, len(reactant)):
        # print(f"react:{reactant[i]},product:{product[i]}"
        combined = [[int(x), int(y), {"index": i}] for x in reactant[i] for y in product[i]]
        # print(len(reactant))
        # print(combined)
        combined_list.append(combined)
        G.add_edges_from(combined)
    return G, combined_list


def plot(G, combined_list, caption):
    rgb = get_rgb()
    plt.figure(figsize=(8, 10))
    # 获取节点位置信息
    # pos = nx.get_node_attributes(G, 'pos')
    pos = nx.kamada_kawai_layout(G)
    print(pos)
    # 获取节点标签信息
    # labels = nx.get_node_attributes(G, 'label')
    labels = nx.get_node_attributes(G, 'label')
    # 绘制节点
    ax = plt.gca()
    nx.draw_networkx_nodes(G, pos, node_size=100, node_shape="o", node_color="white")
    # 绘制边
    for i in combined_list:
        # print(i)
        for j in i:
            # print(j)
            x1, x2 = pos[j[0]][0], pos[j[1]][0]
            y1, y2 = pos[j[0]][1], pos[j[1]][1]
            if y1 != y2 and x1 != x2:
                k = ((x2 - x1) / (y2 - y1)) / 10
                if y1 - 3 > 0 and y2 - 3 > 0:
                    k = -3 * ((x2 - x1) / (y2 - y1)) / 9
                if (x2 - x1) ** 2 < 1:
                    k = -2.5 * ((x2 - x1) / (y2 - y1)) / 10
                a = k * (y2 - y1) + (x2 + x1) / 5
            elif y1 == y2:
                k = -0.5
                a = k * (y2 + y1) + (x2 + x1) / 10
            elif x1 == x2:
                a = -1 / (x2 - 0.5)

            ax.annotate("",
                        xy=(x2, y2),
                        xycoords="data",
                        xytext=(x1, y1),
                        textcoords='data',
                        arrowprops=dict(arrowstyle="->", color=rgb[j[2]["index"]],
                                        shrinkA=20, shrinkB=20,
                                        patchA=None, patchB=None,
                                        connectionstyle=f"arc3,rad={a}"
                                        ),

                        )
    # 绘制节点标签
    nx.draw_networkx_labels(G, pos, labels)
    plt.axis('off')  # 关闭坐标轴
    plt.title(caption)
    # 显示图形
    plt.show()


import networkx as nx
from matplotlib import pyplot as plt

if __name__ == '__main__':
    path = 'CRN.txt'
    species_info, reactant, product = from_stem(path)
    species_list, reduct_degree_dict = get_reduct_degree_dict(species_info)
    G, combined_list = set_G(species_list, reduct_degree_dict)
    plot(G, combined_list, 'stem6')
