import numpy as np
from rdkit import Chem
import colorsys
from rdkit.Chem import MolStandardize
# ₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ ₊ ₋ ₌ ₍ ₎ ₐ ₑ ₒ ₓ ₔ ₕ ₖ ₗ ₘ ₙ ₚ ₛ |⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ ⁼ ˂ ˃ ⁽ ⁾ ˙ * ′ ˙ ⁿ º
class Species(object):
    smiles_struc = {"O=C=O": "CO2",
            "[H][H]": "H2",
            "O=CO": "HCOOH",
            "CO": "CH3OH",
            "C": "CH4",
            "[H+]": "H",
            "O=C[O-]": "HCOO",
            "[OH-]": "OH",
            "[CH3-]": "CH3",
            "C[O-]": "CH3O",
            "[C-4]": "C",
            "[O-2]": "O",
            "[C-]#[O+]": "CO",
            "O=[C-]O": "COOH",
            "[O-]CO": "H2COOH",
            "O=O": "O2",
            "[C-2]=[O+][O-]": "COO",
            "[CH-3]": "CH",
            "[C-2]=[OH+]": "COH",
            "[O-]O": "HO2",
            "[CH-]=[O+][O-]": "CH=OO",
            "[C-2]=[O+]O": "C=OOH",
            "[CH-]=O": "CHO",
            "[O-]C[O-]": "H2COO",
            "[CH2-2]": "CH2",
            "[CH-]=[OH+]": "CHOH",
            "C=O": "CH2O",
            "OO": "H2O2",
            "[CH-]=[O+]O": "CHOOH",
            "[CH2-]O[O-]": "CH2OO",
            "[CH2-]O": "CH2OH",
            "[CH2-]OO": "CH2OOH",
            "CO[O-]": "CH3OO",
            "N": "NH3",
            "[N-3]" : "N",
            "[NH2-]": "NH2",
            "[NH-2]" : "NH",
            "N#N" : "N2",
            "O-2" : "O",
            "[O-]O" :" OOH",
            "[N-]=N":"N2H",
            "[N-]=[NH2+]": "NNH2",
            "[N-]=O": "NO",
            "N=N": "N2H2",
            "[NH-]N": "NHNH2",
            "N=O": "HNO",
            "NN" : "N2H4",
            "N[O-]": "H2NO",
            "[N-]=[OH+]":'NOH',
            "OO" : "H2O2",
            "[NH-]O": "NHOH",
            "NO": "NH2OH"
        }

    def __init__(self, index=0, nature="(reactants)", chem_form="CO2", smiles="O=C=O", reduct_degree=None,
                 struc_form=None, smole=False, atom_num=0):
        self.atom_num = atom_num
        self.smole = smole
        self.nature = nature
        self.index = index
        self.chem_form = chem_form
        self.smiles = self.normalized_smiles(smiles)
        self.reduct_degree = reduct_degree
        self.struc_form = struc_form
        if self.reduct_degree is None:
            self.set_reduct_degree()
        if self.struc_form is None and self.smiles in self.smiles_struc:
            self.struc_form = self.smiles_struc[self.smiles]
        else:
            self.struc_form = self.chem_form
        self.decide_smol()

    def set_reduct_degree(self):
        mol = Chem.MolFromSmiles(self.smiles)  # RDkit读入SMILE表达式初始化mol对象
        mol1 = Chem.AddHs(mol)  # 补充氢原子
        atoms = mol1.GetAtoms()  # 获取mol1中原子对象列表
        reduct_degree = 0
        for atom in atoms:
            if atom.GetSymbol() == "O":
                reduct_degree += 2
            if atom.GetSymbol() == "H":
                reduct_degree += -1
        self.reduct_degree = reduct_degree


    @staticmethod
    def normalized_smiles(smi):
        """
        归一化分子式，选择最大的分子，去除负电荷
        """
        normizer = MolStandardize.normalize.Normalizer()
        # lfc = MolStandardize.fragment.LargestFragmentChooser() # 选择最大的分子
        # uc = MolStandardize.charge.Uncharger() # 去除负电荷
        
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = normizer.normalize(mol)
            # mol = lfc.choose(mol)
            # mol = uc.uncharge(mol)
            smi_normalized = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
            print(f"{smi} -> {smi_normalized}")
            return smi_normalized
        else:
            return None
        
    def decide_smol(self):
        decide = True
        mol = Chem.MolFromSmiles(self.smiles)  # RDkit读入SMILE表达式初始化mol对象
        mol1 = Chem.AddHs(mol)  # 补充氢原子
        atoms = mol1.GetAtoms()  # 获取mol1中原子对象列表
        for atom in atoms:
            if atom is not None:
                self.atom_num += 1
            if atom.GetSymbol() == "C" or atom.GetSymbol() == "N":
                decide = False
        if self.atom_num <= 2 and decide == True:
            self.smole = True

    @classmethod
    def initialize(cls, list):
        # species_info, reactant, product = from_stem(file)
        return cls(list[0], list[1], list[2], list[3])  ## index,type,chemical formular, smiles


def from_stem(episode:list):
    #读取CRN文件信息
    species_info, reactant, product = [], [], []
    stem = episode
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
        species_list.append(Species.initialize(species_info[i]))    #生成species的对象
        # print(species_list[i].reduct_degree)
        # {reduct_degree:num}
        if species_list[i].reduct_degree not in reduct_degree_dict:
            reduct_degree_dict[species_list[i].reduct_degree] = 1
        else:
            reduct_degree_dict[species_list[i].reduct_degree] += 1
    return species_list, reduct_degree_dict


def get_rgb(n, total_colors=100):
    """ Convert an integer to an RGB color with strong contrast within the given range. """
    # Normalize the input integer to the range [0, 1]
    normalized = n / total_colors
    
    # Convert the normalized value to HSV color space for better color distribution
    hue = normalized
    saturation = 0.8 + 0.2 * (n % 2)  # Alternate between slightly different saturation levels
    value = 0.8 + 0.2 * (n % 3)       # Alternate between slightly different brightness levels
    
    # Convert HSV to RGB
    rgb = colorsys.hsv_to_rgb(hue, saturation, value)
    
    # Scale RGB values to 0-255
    rgb = np.array(tuple(int(min(c * 255, 255)) for c in rgb))# Ensure values are within 0-255
    
    return rgb/255


def set_G(reactant, product, species_list, reduct_degree_dict):
    G = nx.Graph()
    smole_list = []
    smole_label = {}
    reduct_degree_list = list(reduct_degree_dict.keys())
    reduct_degree_list.sort()
    # set the location of nodes
    for i, species in enumerate(species_list):
        node_color = "#3fc1c9"
        # print(species.nature)
        if species.smole:
            smole_list.append(species.index)
            smole_label[species.index] = species.struc_form
        else:
            if species.nature == "(Reactant)":
                node_color = "#fc5185"
            elif species.nature == "(Product)":
                node_color = "#fc5185"
            G.add_node(species.index,
                       pos=(reduct_degree_dict[species.reduct_degree],
                            reduct_degree_list.index(species.reduct_degree)),
                       label=species.struc_form,
                       color=node_color)
            reduct_degree_dict[species.reduct_degree] -= 2
    # set the edges
    combined = None
    combined_list = []
    for i in range(0, len(reactant)):
        reaction_type = ""
        current_reactant = reactant[i][:]
        currenet_product = product[i][:]
        for j in range(0, len(reactant[i])):
            if int(reactant[i][j]) in smole_list:
                reaction_type = int(reactant[i][j])
                current_reactant.remove(reactant[i][j])
            else:
                if reaction_type == "":
                    reaction_type = int(reactant[i][j])
        #print(f"react:{reactant[i]},product:{product[i]}")
        # print(len(reactant))
        # print(combined)
        for k in range(0, len(product[i])):
            if int(product[i][k]) in smole_list:
                currenet_product.remove(product[i][k])
        # if current_reactant == []:
        #     current_reactant = list(set(reactant[i][:]))
        if currenet_product == [] or current_reactant == []:
                current_reactant==[]
                currenet_product==[]
        # if currenet_product == []:
        #     currenet_product = list(set(product[i][:]))
        combined = [[int(x), int(y), {"type": reaction_type}] for x in current_reactant for y in currenet_product]
        combined_list.append(combined)
        G.add_edges_from(combined)
    return G, combined_list, smole_label


def set_equation(species_list, reatant, product):
    #产生反应式
    equation = []
    species_dict = {}
    for i in species_list:
        # print(i.index)
        species_dict[str(i.index)] = i
        # print(species_dict.keys())
    for i,j in zip(reatant,product):
        if len(i) == 1 and len(j) == 1:
            if species_dict[i[0]].smole == True and species_dict[j[0]].smole == True:
                equation.append(species_dict[i[0]].struc_form + " = " + species_dict[j[0]].struc_form)
        if len(i) == 2 and len(j) == 1:
            if species_dict[i[0]].smole == True and species_dict[i[1]].smole == True and species_dict[j[0]].smole == True:
                equation.append(species_dict[i[0]].struc_form + "+" + species_dict[i[1]].struc_form + " = " +
                                species_dict[j[0]].struc_form)
        if len(i) == 1 and len(j) == 2:
            if species_dict[i[0]].smole == True and species_dict[j[0]].smole == True and species_dict[j[1]].smole == True:
                equation.append(species_dict[i[0]].struc_form + " = " + species_dict[j[0]].struc_form + "+" +
                                species_dict[j[1]].struc_form)
        if len(i) == 2 and len(j) == 2:
            if species_dict[i[0]].smole == True and species_dict[i[1]].smole == True and species_dict[j[0]].smole == True and species_dict[j[1]].smole == True:
                equation.append(species_dict[i[0]].struc_form + "+" + species_dict[i[1]].struc_form + " = " +
                                species_dict[j[0]].struc_form + "+" + species_dict[j[1]].struc_form)
    return equation


def plot(G, combined_list, caption, smole_label, equation, show_equation=True):
    plt.figure(figsize=(10, 12),dpi=300)
    # 获取节点位置信息
    # pos = nx.get_node_attributes(G, 'pos')
    pos = nx.kamada_kawai_layout(G)   #设置节点布局
    # print("pos:",pos)
    # pos = nx.spectral_layout(G)
    # pos = nx.spring_layout(G)
    # pos = nx.circular_layout(G)
    def shift_y_values(pos, delta_x,delta_y):
        label_pos={}
        for key, value in pos.items():
            label_x = value[0] + delta_x
            label_y = value[1] + delta_y
            if label_x <= -0.05:
                label_x = value[0] - delta_x
            if label_y >= 0.05:
                label_y = value[1] - delta_y
            label_pos[key] = np.array([label_x,label_y])
        return label_pos

    delta_x = -0.00
    delta_y = 0.00
    label_pos = shift_y_values(pos,delta_x, delta_y)
    color = nx.get_node_attributes(G, 'color')
    for i in smole_label.keys():
        if len(color) < len(pos):
            if i in pos and i not in color:
                G.add_node(i, label=smole_label[i], color="#3399FF")
                color[i] = "#3399FF"
    # 获取节点标签信息
    labels = nx.get_node_attributes(G, 'label')
    # print(labels)
    # print(smole_label)
    # 绘制节点
    ax = plt.gca()
    nx.draw_networkx_nodes(G, pos, node_size=1000, node_shape="o", node_color=[x for x in color.values()])
    # 绘制边
    color_line = None
    x, y = 0.8, 0.55
    for i in combined_list:
        # print(i)
        for j in i:
            # print(j)
            x1, x2 = pos[j[0]][0], pos[j[1]][0]
            y1, y2 = pos[j[0]][1], pos[j[1]][1]
            color_line =get_rgb(j[2]["type"], len(combined_list))
            ax.annotate("",
                        xy=(x2, y2),
                        xycoords="data",
                        xytext=(x1, y1),
                        textcoords='data',
                        arrowprops=dict(arrowstyle="-", color=color_line,
                                        shrinkA=20, shrinkB=20,
                                        patchA=None, patchB=None,
                                        connectionstyle=f"arc3,rad=0",
                                        linewidth=4
                                        )

                        )
            if j[2]["type"] in smole_label.keys():
                line1, = ax.plot([], [], label=f'{smole_label[j[2]["type"]]}', linestyle='-',
                                 color=color_line)  # Invisible line with alpha=0
                smole_label.pop(j[2]["type"])
    if show_equation:
        bbox = {'facecolor': 'white', 'alpha': 0.3, 'pad': 5}
        plt.text(x, y, '\n'.join(equation), fontsize=12, ha='center', va='center',color='red',
        bbox =bbox)


    # 绘制节点标签
    #nx.draw_networkx_labels(G, pos, labels,font_size=20, font_weight='bold', bbox={'boxstyle': 'round4,pad=0.1', 'facecolor': 'white', 'edgecolor': 'black', 'linewidth':3})
    nx.draw_networkx_labels(G, label_pos, labels,font_size=12,font_color="black",font_weight='bold')
    plt.axis('off')  # 关闭坐标轴
    # plt.title(caption)
    plt.legend(bbox_to_anchor=(0.9, 0.5), fontsize=12)
    plt.tight_layout()
    # plt.subplots(constrained_layout=True)

def from_output(file):
    ### 提取stem
    fragment = []
    with open(file, 'r', encoding='GB18030') as output:  # Specify the encoding as 'utf-8'
        out_info = list(output)
    mode = "pass"  # ["pass","read","end"]
    for i, line in enumerate(out_info):

        if "involved species" in line:
            mode = "read"
            temp = []
        if mode == "read" and line == "\n":
            mode = "end"
        if mode == "pass":
            continue
        if mode == "read":
            temp.append(line)
            # print('temp add')
        if mode == "end":
            mode = "pass"
            fragment.append(temp)
    return fragment

import networkx as nx
from matplotlib import pyplot as plt
import os
def draw_net():
    filename = 'CRNGenerator_log.txt'
    path = os.path.join(os.getcwd(), filename)
    fragment = from_output(path)
    for i,episode in enumerate(fragment):
        species_info, reactant, product = from_stem(episode)
        species_list, reduct_degree_dict = get_reduct_degree_dict(species_info)
        G, combined_list, smol_label = set_G(reactant, product, species_list, reduct_degree_dict)
        equation = set_equation(species_list, reactant, product)
        plot(G, combined_list, f"stem{i}", smol_label, equation, show_equation=True)
        # fold_path = "AutoCat\\utils\\crn\\img_subcrn\\" # dir
        file_name = f"stem{i}.png"
        # save_path = fold_path + file_name
        plt.savefig(file_name)

if __name__ == '__main__':
    draw_net()