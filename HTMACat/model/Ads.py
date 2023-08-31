"""Created on Sat Mar 18 09:00:12 2023.

@author: YuxiaoLan
"""

from HTMACat.catkit.gen.adsorption import Builder
from HTMACat.catkit.gratoms import *
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import math
from HTMACat.Extract_info import *
from HTMACat.model.Substrate import Slab
from HTMACat.catkit.gen.adsorption import AdsorptionSites
from HTMACat.model.Structure import Structure
import networkx.algorithms.isomorphism as iso
from ase import Atoms
from HTMACat.model.Species import init_from_ads, ABS_Species


class Adsorption(Structure):
    def __init__(self, species: list, sites: list, spec_ads_stable=None, substrate=Slab()):
        if spec_ads_stable is None:
            spec_ads_stable = {
                "NH3": [1],
                "NH2": [2],
                "NH": [2, 4],
                "N": [2, 4],
                "O": [2, 4],
                "OH": [2, 4],
                "NO": [2, 4],
                "H2O": [1],
                "H": [2, 4],
            }
        assert isinstance(species, list), "species should be a list of Species class"
        assert isinstance(species[0], ABS_Species), "species should be a list of Species class"
        assert isinstance(sites, list), "sites should be a list"
        #assert sites[0] in ['1', '2'], 'Supports only "1" "2" adsorption sites type for ads!'
        self.species = species
        self.sites = sites
        self.spec_ads_stable = spec_ads_stable
        self.substrate = substrate

    def set_species(self, species):
        self.species = species

    def add_species(self, species):
        self.species.append(species)

    def add_sites(self, sites):
        self.sites.append(sites)

    def get_sites(self):
        return self.sites # ' '.join(self.sites)

    def out_file_name(self):
        ads = []
        ads_str = self.species[0].out_file_name()
        substrate_str = self.substrate.out_file_name()
        vasp_file_str: str = "_".join([substrate_str, ads_str])
        return vasp_file_str

    def out_print(self):
        species_str = []
        for i in range(len(self.species)):
            species_str.append(self.species[i].out_print())

        species_str = " and ".join(species_str)
        substrate_str = self.substrate.out_print()
        print_str = f"{species_str} adsorption on {substrate_str}"
        return print_str

    def construct(self):
        if self.get_sites()[0] == '1':
            if len(self.get_sites()) == 1:
                slabs_ads = self.Construct_single_adsorption()
            else:
                ele = ''.join(self.get_sites()[1:]) ### wzj 20230518
                slabs_ads = self.Construct_single_adsorption(ele=ele)
        elif self.get_sites()[0] == '2':
            slabs_ads = self.Construct_double_adsorption()
        else:
            raise ValueError('Supports only "1" "2" adsorption sites for ads!')

        if self.substrate.is_dope():
            slabs_ads = self.remove_same(slabs_ads)

        return slabs_ads

    def remove_same(self, slabs_ads):
        # To choose the config whose neighbor list includes the doped atoms
        p1 = self.substrate.property["p1"]
        p1_symb = self.substrate.property["p1_symb"]
        slabs_ads_near = []
        for slb in slabs_ads:
            (
                bind_adatoms,
                bind_adatoms_symb,
                bind_type_symb,
                adspecie,
                bind_surfatoms,
                bind_surfatoms_symb,
            ) = get_binding_adatom(slb)
            if self.get_sites() == "1" and (set(p1_symb) & set(bind_surfatoms_symb[0])):
                slabs_ads_near += [slb]
            elif (
                self.get_sites() == "2"
                and (set(p1_symb) & set(bind_surfatoms_symb[0]))
                or (set(p1_symb) & set(bind_surfatoms_symb[0]))
            ):
                slabs_ads_near += [slb]
        return slabs_ads_near
    
    def vec_to_nearest_neighbor_site(self, slab, site_coords):
        # site_coords: a list containing 1 (ads) or 2 (coads) coord(s) of adsorption site(s)
        # return: the vector from the adsorption site of specie A to its nearest neighboring adsorption site
        # [The v_nn of specie A is always the v_nn of specie B too !!!]
        coord_images = site_coords
        num_origin_sites = len(coord_images)
        for d,dim in enumerate(slab._pbc):
            if dim:
                for i in range(num_origin_sites):
                    coord_images.append(coord_images[i]+slab._cellobj[d]) # coord_images[0]: the coord of specie A
                    coord_images.append(coord_images[i]-slab._cellobj[d])
        imagesites_distances = [np.sqrt(np.sum(np.square(v-coord_images[0]))) for v in coord_images]
        idx_tmp = np.argmin(imagesites_distances[1:]) + 1
        v_nn = coord_images[idx_tmp] - coord_images[0]
        #print('imagesites_distances--:', imagesites_distances, idx_tmp)
        #print('vec_to_neigh_imgsite--:', v_nn)
        return v_nn
    
    def dist_of_nearest_diff_neigh_site(self, slab, site_coords):
        # site_coords: a list of 2 site coords (different species)
        coord_images = site_coords
        num_origin_sites = len(coord_images)
        for d,dim in enumerate(slab._pbc):
            if dim:
                coord_images.append(coord_images[1]+slab._cellobj[d]) # coord_images[0]: the coord of specie A
                coord_images.append(coord_images[1]-slab._cellobj[d])
        imagesites_distances = [np.sqrt(np.sum(np.square(v[:2]-coord_images[0][:2]))) for v in coord_images]
        d = np.min(imagesites_distances[1:])
        return d

    def Construct_single_adsorption(self, ele=None):
        # generate surface adsorption configuration
        slab_ad = []
        slabs = self.substrate.construct()
        for i, slab in enumerate(slabs):
            site = AdsorptionSites(slab)
            coordinates = site.get_coordinates()
            builder = Builder(slab)
            ads_use, ads_use_charges = self.species[0].get_molecule()
            if not ele is None:
                if ele == '+':
                    bond_atom_ids = np.where(np.array(ads_use_charges)>0)[0]
                elif ele == '-':
                    bond_atom_ids = np.where(np.array(ads_use_charges)<0)[0]
                else:
                    chemical_symbols = np.array(ads_use.get_chemical_symbols())
                    bond_atom_ids = np.where(chemical_symbols==ele)[0]
                for j, coord in enumerate(coordinates):
                    vec_to_neigh_imgsite = self.vec_to_nearest_neighbor_site(slab=slab, site_coords=[coord])
                    for bond_id in bond_atom_ids:
                        slab_ad += [builder._single_adsorption(ads_use, bond=bond_id, site_index=j,
                                                               rotation_mode ='vertical to vec_to_neigh_imgsite',
                                                               rotation_args ={'vec_to_neigh_imgsite':vec_to_neigh_imgsite},
                                                               direction_mode='decision_boundary')]
                        #if len(bond_atom_ids) > 1:
                        #    slab_ad += [builder._single_adsorption(ads_use, bond=bond_id, site_index=j, direction_mode='decision_boundary', direction_args=bond_atom_ids)]
            else:
                for j, coord in enumerate(coordinates):
                    vec_to_neigh_imgsite = self.vec_to_nearest_neighbor_site(slab=slab, site_coords=[coord])
                    slab_ad += [builder._single_adsorption(ads_use, bond=0, site_index=j,
                                                           rotation_mode ='vertical to vec_to_neigh_imgsite',
                                                           rotation_args ={'vec_to_neigh_imgsite':vec_to_neigh_imgsite})]
        return slab_ad

    def Construct_double_adsorption(self):
        slab_ad = []
        slabs = self.substrate.construct()
        for i, slab in enumerate(slabs):
            site = AdsorptionSites(slab)
            builder = Builder(slab)
            ads_use, ads_use_charges = self.species[0].get_molecule()
            edges = site.get_adsorption_edges()
            for j, edge01 in enumerate(edges):
                slab_ad += [builder._double_adsorption(ads_use, bonds=[0, 1], edge_index=j)]
        return slab_ad

    @classmethod
    def from_input(cls, init_list, substrates, species_dict=None):
        ads = []
        for i in init_list:
            spec1 = init_from_ads(i[0], species_dict)
            sites1 = str(i[1])
            for j in substrates:
                ads.append(cls([spec1], list(sites1), substrate=j))
        return ads


class Coadsorption(Adsorption):
    def __init__(self, species: list, sites: list, spec_ads_stable=None, substrate=Slab()):
        super().__init__(species, sites, spec_ads_stable, substrate)
        assert len(species) == 2, "Coads need Two adsorption Species, but %d was given" % len(
            species
        )

    def construct(self):
        ### Change by RxChen, 2023 7 26: ZhaojieWang 20230828:
        if ([self.get_sites()[0][0], self.get_sites()[1][0]] == ['1','1']):
            ele = [''.join(self.get_sites()[0][1:]),''.join(self.get_sites()[1][1:])]
            slabs_ads = self.Construct_coadsorption_11(ele=ele)
        elif ([self.get_sites()[0][0], self.get_sites()[1][0]] == ['1','2']):
            slabs_ads = self.Construct_coadsorption_12()
        elif ([self.get_sites()[0][0], self.get_sites()[1][0]] == ['2','2']):
            slabs_ads = self.Construct_coadsorption_22()
        else:
            raise ValueError("Supports only '1' or '2' adsorption sites for coads!")### end
        if self.substrate.is_dope():
            slabs_ads = self.remove_same(slabs_ads)
        return slabs_ads

    def out_file_name(self):
        ads = []
        for every_species in self.species:
            ads.append(every_species.out_file_name())
        ads_str = "_".join(ads)
        substrate_str = self.substrate.out_file_name()
        vasp_file_str = "_".join([substrate_str, ads_str])
        return vasp_file_str

    def remove_same(self, slabs_ads):
        # To choose the config whose neighbor list includes the doped atoms
        p1 = self.substrate.property["p1"]
        p1_symb = self.substrate.property["p1_symb"]
        slabs_ads_near = []
        for slb in slabs_ads:
            (
                bind_adatoms,
                bind_adatoms_symb,
                bind_type_symb,
                adspecie,
                bind_surfatoms,
                bind_surfatoms_symb,
            ) = get_binding_adatom(slb)
            bind_surfatoms_symb_all = sum(bind_surfatoms_symb, [])
            if set(p1_symb) & set(bind_surfatoms_symb_all):
                slabs_ads_near += [slb]
        return slabs_ads_near

    def Construct_coadsorption_11(self, dist_min_xoy='default', ele=['','']):
        slab_ad = []
        ads_type = self.spec_ads_stable
        dis_inter = self.substrate.get_dis_inter()
        slabs = self.substrate.construct()
        if 'default' == dist_min_xoy:
            la = np.sqrt(np.dot(slabs[0]._cellobj[0], slabs[0]._cellobj[0]))
            lb = np.sqrt(np.dot(slabs[0]._cellobj[1], slabs[0]._cellobj[1]))
            dist_min_xoy = min(0.4*min(la,lb), 4.0)
        for i, slab in enumerate(slabs):
            site01 = AdsorptionSites(slab)
            coordinates01 = site01.get_coordinates() ### add
            builder01 = Builder(slab)
            # generate surface adsorption configuration
            ads1_use, ads1_use_charges = self.species[0].get_molecule()
            ads2_use, ads2_use_charges = self.species[1].get_molecule()
            bond_atom_ids_list = [[],[]]
            if not ele[0] in ['', '+', '-']:
                chemical_symbols01 = np.array(ads1_use.get_chemical_symbols())
                bond_atom_ids_list[0] = np.where(chemical_symbols01==ele[0])[0]
            elif '+' == ele[0]:
                bond_atom_ids_list[0] = np.where(np.array(ads1_use_charges)>0)[0]
            elif '-' == ele[0]:
                bond_atom_ids_list[0] = np.where(np.array(ads1_use_charges)<0)[0]
            else:
                bond_atom_ids_list[0] = np.array([0])
            if not ele[1] in ['', '+', '-']:
                chemical_symbols02 = np.array(ads2_use.get_chemical_symbols())
                bond_atom_ids_list[1] = np.where(chemical_symbols02==ele[1])[0]
            elif '+' == ele[1]:
                bond_atom_ids_list[1] = np.where(np.array(ads2_use_charges)>0)[0]
            elif '-' == ele[1]:
                bond_atom_ids_list[1] = np.where(np.array(ads2_use_charges)<0)[0]
            else:
                bond_atom_ids_list[1] = np.array([0])
            for k, coord01 in enumerate(coordinates01):
                slab_tmp = builder01._single_adsorption(ads1_use, bond=0, site_index=k, auto_construct=True,
                                                        enable_rotate_xoy=False,
                                                        direction_mode='decision_boundary') # 此处bond参数随意，有吸附即可
                site02 = AdsorptionSites(slab_tmp)
                coordinates02 = site02.get_coordinates()
                ###print('len(coordinates01), len(coordinates02) =', len(coordinates01), len(coordinates02))
                # 第一个物种放到表面上之后打破了对称性才使得第二个物种的可吸附位点情况变多
                for j, coord02 in enumerate(coordinates02):
                    vec_to_neigh_imgsite = self.vec_to_nearest_neighbor_site(slab=slab, site_coords=[coord01,coord02])
                    dis = np.linalg.norm(coord01 - coord02)
                    dis_xoy = self.dist_of_nearest_diff_neigh_site(slab_tmp, [coord01,coord02])
                    if dis < float(dis_inter[0]):
                        continue
                    elif dis > float(dis_inter[1]):
                        continue
                    elif dis_xoy < dist_min_xoy: # 控制两吸附位点的距离（在xoy平面内）不能过小
                        continue
                    else:
                        # 确定最近邻位点连线后正式构建前一个物种的吸附表面
                        for bond_id01 in bond_atom_ids_list[0]:
                            for bond_id02 in bond_atom_ids_list[1]:
                                slab = builder01._single_adsorption(ads1_use, bond=bond_id01, site_index=k, auto_construct=True,
                                                                    enable_rotate_xoy=True,
                                                                    rotation_mode ='vertical to vec_to_neigh_imgsite',
                                                                    rotation_args ={'vec_to_neigh_imgsite':vec_to_neigh_imgsite},
                                                                    direction_mode='decision_boundary')
                                builder02 = Builder(slab)
                                slab_ = builder02._single_adsorption(
                                            ads2_use, bond=bond_id02, site_index=j, auto_construct=True,
                                            enable_rotate_xoy=True,
                                            rotation_mode ='vertical to vec_to_neigh_imgsite',
                                            rotation_args ={'vec_to_neigh_imgsite':vec_to_neigh_imgsite},
                                            direction_mode='decision_boundary'
                                        )
                                slab_ad += [slab_]

        typ = {None: 0, "top": 1, "bri": 2, "fcc": 3, "hcp": 3, "4-fold": 4}
        # view(slab_ad)
        slab_ad_final = []
        for j, adslab in enumerate(slab_ad):
            (
                bind_adatoms,
                bind_adatoms_symb,
                adspecie,
                bind_type_symb,
                bind_surfatoms,
                bind_surfatoms_symb,
            ) = get_binding_adatom(adslab)
            adspecie_tmp, bind_type_symb_tmp = [], []
            for k, spe in enumerate(adspecie):
                if spe in ads_type.keys():
                    adspecie_tmp += [spe]
                    bind_type_symb_tmp += [bind_type_symb[k]]
            if len(adspecie_tmp) < 2:
                # print('Can not identify the config!')
                slab_ad_final += [adslab]
            elif typ.get(bind_type_symb_tmp[0]) in ads_type.get(adspecie_tmp[0]) and typ.get(
                bind_type_symb_tmp[1]
            ) in ads_type.get(adspecie_tmp[1]):
                slab_ad_final += [adslab]
        return slab_ad_final

    def Construct_coadsorption_12(self):
        slab_ad = []
        slabs = self.substrate.construct()
        dis_inter = self.substrate.get_dis_inter()
        for i, slab in enumerate(slabs):
            site01 = AdsorptionSites(slab)
            # print(site01.get_symmetric_sites())
            builder01 = Builder(slab)
            coordinate01 = site01.get_coordinates()
            # generate surface adsorption configuration
            ads1_use, ads1_use_charges = self.species[0].get_molecule()
            ads2_use, ads2_use_charges = self.species[1].get_molecule()
            for k, coord01 in enumerate(coordinate01):
                slab = builder01._single_adsorption(ads1_use, bond=0, site_index=k)
                # after adsorbing an atoms
                # site analysis and surface builder
                site02 = AdsorptionSites(slab)
                builder02 = Builder(slab)
                # obtain the edges and ID and coordinates of all sites
                edge02 = site02.get_adsorption_edges()
                site_coord02 = site02.get_coordinates(unique=False)
                site_type02 = site02.get_periodic_sites(screen=True)
                dic_site = {}
                keys = [site_type02[i] for i in range(len(site_type02))]
                values = [site_coord02[i] for i in range(len(site_type02))]
                dic_site = dict(zip(keys, values))
                # emiliate the too near site
                for j, edge in enumerate(edge02):
                    coord10 = dic_site.get(edge[0])
                    coord11 = dic_site.get(edge[1])
                    if coord10 is not None:
                        dis1 = np.linalg.norm(coord01 - coord10)
                    else:
                        dis1 = dis_inter[0] + 0.1
                    if coord11 is not None:
                        dis2 = np.linalg.norm(coord01 - coord11)
                    else:
                        dis2 = dis_inter[1] + 0.1
                    dis = min(dis1, dis2)

                    if dis < dis_inter[0]:
                        continue
                    elif dis > dis_inter[1]:
                        continue
                    else:
                        slab_ad += [
                            builder02._double_adsorption(ads2_use, bonds=[0, 1], edge_index=j)
                        ]
        return slab_ad

    def Construct_coadsorption_22(self):
        slab_ad = []
        dis_inter = self.substrate.get_dis_inter()
        slabs = self.substrate.construct()
        for i, slab in enumerate(slabs):
            site01 = AdsorptionSites(slab)
            # print(site01.get_symmetric_sites())
            builder01 = Builder(slab)
            # coordinate01 = site01.get_coordinates()
            edge01 = site01.get_adsorption_edges()
            site_coord01 = site01.get_coordinates(unique=False)
            site_type01 = site01.get_periodic_sites(screen=True)
            dic_site01 = {}
            keys = [site_type01[i] for i in range(len(site_type01))]
            values = [site_coord01[i] for i in range(len(site_type01))]
            dic_site01 = dict(zip(keys, values))
            # generate surface adsorption configuration
            ads1_use, ads1_use_charges = self.species[0].get_molecule()
            ads2_use, ads2_use_charges = self.species[1].get_molecule()
            for k, edge01 in enumerate(edge01):
                slab = builder01._double_adsorption(ads1_use, bonds=[0, 1], edge_index=k)
                coord00 = dic_site01.get(edge01[0])
                coord01 = dic_site01.get(edge01[1])
                ## after adsorbing an atoms
                # site analysis and surface builder
                site02 = AdsorptionSites(slab)
                builder02 = Builder(slab)
                # obtain the edges and ID and coordinates of all sites
                edge02 = site02.get_adsorption_edges()
                site_coord02 = site02.get_coordinates(unique=False)
                site_type02 = site02.get_periodic_sites(screen=True)
                dic_site02 = {}
                keys = [site_type02[i] for i in range(len(site_type02))]
                values = [site_coord02[i] for i in range(len(site_type02))]
                dic_site02 = dict(zip(keys, values))
                ## emiliate the too near site
                for j, edge in enumerate(edge02):
                    coord10 = dic_site02.get(edge[0])
                    coord11 = dic_site02.get(edge[1])
                    if not coord10 is None:
                        dis1 = np.linalg.norm(coord01 - coord10)
                    else:
                        dis1 = dis_inter[0] + 0.1
                    if not coord11 is None:
                        dis2 = math.hypot(
                            coord00[0] - coord11[0],
                            coord00[1] - coord11[1],
                            coord00[2] - coord11[2],
                        )
                    else:
                        dis2 = dis_inter[1] + 0.1
                    dis = min(dis1, dis2)

                    if dis < dis_inter[0]:
                        continue
                    elif dis > dis_inter[1]:
                        continue
                    else:
                        slab_ad += [
                            builder02._double_adsorption(ads2_use, bonds=[0, 1], edge_index=j)
                        ]
        return slab_ad

    @classmethod
    def from_input(cls, init_list, substrates, species_dict=None):
        ads = []
        for i in init_list:
            spec1 = init_from_ads(i[0], species_dict)
            spec2 = init_from_ads(i[1], species_dict)
            sites1 = str(i[2])
            sites2 = str(i[3])
            for j in substrates:
                ads.append(cls([spec1, spec2], [sites1, sites2], substrate=j))
        return ads


def ads_from_input(ads_model, substrate, species_dict=None):
    ads = []
    new_ads = []
    ads_type_list = ["ads", "coads"]
    for key, value in ads_model.items():
        if key == "coads":
            new_ads = Coadsorption.from_input(value, substrate, species_dict)
        elif key == "ads":
            new_ads = Adsorption.from_input(value, substrate, species_dict)
        else:
            msg = ",".join(ads_type_list)
            warn_msg = (
                "Only support species type: %s, Your input %s part in Species will be dismiss"
                % (msg, key)
            )
            raise Warning(warn_msg)
        ads = ads + new_ads
    return ads
