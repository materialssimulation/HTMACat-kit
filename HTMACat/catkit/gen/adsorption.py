from . import defaults
from . import utils
from . import symmetry
from ase.build import molecule, fcc111
from ase.io import write
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations

from copy import deepcopy
import HTMACat.catkit as catkit
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import numpy as np
import scipy, copy

radii = defaults.get("radii")


class AdsorptionSites:
    """Adsorption site object."""

    def __init__(self, slab, surface_atoms=None, tol=1e-5):
        """Create an extended unit cell of the surface sites for use in identifying other sites.

        Parameters
        ----------
        slab : Gatoms object
            The slab associated with the adsorption site network to be
            attached.
        tol : float
            Absolute tolerance for floating point errors.
        """
        index, coords, offsets = utils.expand_cell(slab, cutoff=5.0)
        if surface_atoms is None:
            surface_atoms = slab.get_surface_atoms()
        if surface_atoms is None:
            raise ValueError("Slab must contain surface atoms")

        extended_top = np.where(np.in1d(index, surface_atoms))[0]

        self.tol = tol
        self.coordinates = coords[extended_top].tolist()
        self.connectivity = np.ones(extended_top.shape[0]).tolist()
        self.r1_topology = [[i] for i in np.arange(len(extended_top))]
        self.index = index[extended_top]

        sites = self._get_higher_coordination_sites(coords[extended_top])
        self.r2_topology = sites["top"][2]

        # Put data into array format
        selection = ["bridge", "hollow", "4fold"]
        for i, k in enumerate(selection):
            coordinates, r1top, r2top = sites[k]

            if k in ["hollow", "4fold"]:
                r2top = [[] for _ in coordinates]

            self.connectivity += (np.ones(len(coordinates)) * (i + 2)).tolist()
            self.coordinates += coordinates
            self.r1_topology += r1top
            self.r2_topology += r2top

        self.coordinates = np.array(self.coordinates)
        self.connectivity = np.array(self.connectivity, dtype=int)
        self.r1_topology = np.array(self.r1_topology, dtype=object)
        self.r2_topology = np.array(self.r2_topology, dtype=object)
        self.frac_coords = np.dot(self.coordinates, np.linalg.pinv(slab.cell))
        self.slab = slab

        screen = (
            (self.frac_coords[:, 0] > 0 - self.tol)
            & (self.frac_coords[:, 0] < 1 - self.tol)
            & (self.frac_coords[:, 1] > 0 - self.tol)
            & (self.frac_coords[:, 1] < 1 - self.tol)
        )

        self.screen = screen
        self._symmetric_sites = None

    def get_connectivity(self, unique=True):
        """Return the number of connections associated with each site."""
        if unique:
            sel = self.get_symmetric_sites()
        else:
            sel = self.get_periodic_sites()

        return self.connectivity[sel]

    def get_coordinates(self, unique=True):
        """Return the 3D coordinates associated with each site."""
        if unique:
            sel = self.get_symmetric_sites()
        else:
            sel = self.get_periodic_sites()

        return self.coordinates[sel]

    def get_topology(self, unique=True):
        """Return the indices of adjacent surface atoms."""
        topology = [self.index[top] for top in self.r1_topology]
        topology = np.array(topology, dtype=object)
        if unique:
            sel = self.get_symmetric_sites()
        else:
            sel = self.get_periodic_sites()

        return topology[sel]

    def _get_higher_coordination_sites(self, top_coordinates, allow_obtuse=True):
        """Find all bridge and hollow sites (3-fold and 4-fold) given an input slab based Delaunay
        triangulation of surface atoms of a super-cell.

        TODO: Determine if this can be made more efficient by
        removing the 'sites' dictionary.

        Parameters
        ----------
        top_coordinates : ndarray (n, 3)
            Cartesian coordinates for the top atoms of the unit cell.

        Returns
        -------
        sites : dict of 3 lists
            Dictionary sites containing positions, points, and neighbor lists.
        """
        sites = {
            "top": [top_coordinates, [], [[] for _ in top_coordinates]],
            "bridge": [[], [], []],
            "hollow": [[], [], []],
            "4fold": [[], [], []],
        }

        dt = scipy.spatial.Delaunay(sites["top"][0][:, :2])
        neighbors = dt.neighbors
        simplices = dt.simplices

        for i, corners in enumerate(simplices):
            cir = scipy.linalg.circulant(corners)
            edges = cir[:, 1:]

            # Inner angle of each triangle corner
            vec = sites["top"][0][edges.T] - sites["top"][0][corners]
            uvec = vec.T / np.linalg.norm(vec, axis=2).T
            angles = np.sum(uvec.T[0] * uvec.T[1], axis=1)

            # Angle types
            right = np.isclose(angles, 0)
            obtuse = angles < -self.tol

            rh_corner = corners[right]
            edge_neighbors = neighbors[i]

            if obtuse.any() and not allow_obtuse:
                # Assumption: All simplices with obtuse angles
                # are irrelevant boundaries.
                continue

            bridge = np.sum(sites["top"][0][edges], axis=1) / 2.0

            # Looping through corners allows for elimination of
            # redundant points, identification of 4-fold hollows,
            # and collection of bridge neighbors.
            for j, c in enumerate(corners):
                edge = sorted(edges[j])

                if edge in sites["bridge"][1]:
                    continue

                # Get the bridge neighbors (for adsorption vector)
                neighbor_simplex = simplices[edge_neighbors[j]]
                oc = list(set(neighbor_simplex) - set(edge))[0]

                # Right angles potentially indicate 4-fold hollow
                potential_hollow = edge + sorted([c, oc])
                if c in rh_corner:
                    if potential_hollow in sites["4fold"][1]:
                        continue

                    # Assumption: If not 4-fold, this suggests
                    # no hollow OR bridge site is present.
                    ovec = sites["top"][0][edge] - sites["top"][0][oc]
                    ouvec = ovec / np.linalg.norm(ovec)
                    oangle = np.dot(*ouvec)
                    oright = np.isclose(oangle, 0)
                    if oright:
                        sites["4fold"][0] += [bridge[j]]
                        sites["4fold"][1] += [potential_hollow]
                        sites["top"][2][c] += [oc]
                else:
                    sites["bridge"][0] += [bridge[j]]
                    sites["bridge"][1] += [edge]
                    sites["bridge"][2] += [[c, oc]]

                sites["top"][2][edge[0]] += [edge[1]]
                sites["top"][2][edge[1]] += [edge[0]]

            if not right.any() and not obtuse.any():
                hollow = np.average(sites["top"][0][corners], axis=0)
                sites["hollow"][0] += [hollow]
                sites["hollow"][1] += [corners.tolist()]

        # For collecting missed bridge neighbors
        for s in sites["4fold"][1]:
            edges = itertools.product(s[:2], s[2:])
            for edge in edges:
                edge = sorted(edge)
                i = sites["bridge"][1].index(edge)
                n, m = sites["bridge"][1][i], sites["bridge"][2][i]
                nn = list(set(s) - set(n + m))

                if len(nn) == 0:
                    continue
                sites["bridge"][2][i] += [nn[0]]

        return sites

    def get_periodic_sites(self, screen=True):
        """Return an index of the coordinates which are unique by periodic boundary conditions.

        Parameters
        ----------
        screen : bool
            Return only sites inside the unit cell.

        Returns
        -------
        periodic_match : ndarray (n,)
            Indices of the coordinates which are identical by
            periodic boundary conditions.
        """
        periodic_match = np.arange(self.frac_coords.shape[0])

        if screen:
            return periodic_match[self.screen]

        coords = self.frac_coords
        periodic = periodic_match.copy()[self.screen]

        for p in periodic:
            matched = utils.matching_sites(self.frac_coords[p], coords)
            periodic_match[matched] = p

        return periodic_match

    def get_symmetric_sites(self, unique=True, screen=True):
        """Determine the symmetrically unique adsorption sites from a list of fractional
        coordinates.

        Parameters
        ----------
        unique : bool
            Return only the unique symmetrically reduced sites.
        screen : bool
            Return only sites inside the unit cell.

        Returns
        -------
        sites : dict of lists
            Dictionary of sites containing index of site
        """
        symmetry_match = self._symmetric_sites

        if symmetry_match is None:
            sym = symmetry.Symmetry(self.slab, tol=self.tol)

            rotations, translations = sym.get_symmetry_operations(affine=False)
            rotations = np.swapaxes(rotations, 1, 2)

            affine = np.append(rotations, translations[:, None], axis=1)
            points = self.frac_coords
            true_index = self.get_periodic_sites(False)

            affine_points = np.insert(points, 3, 1, axis=1)
            operations = np.dot(affine_points, affine)
            symmetry_match = np.arange(points.shape[0])

            for i, j in enumerate(symmetry_match):
                if i != j:
                    continue

                d = operations[i, :, None] - points
                d -= np.round(d)
                dind = np.where((np.abs(d) < self.tol).all(axis=2))[-1]
                symmetry_match[np.unique(dind)] = true_index[i]

            self._symmetric_sites = symmetry_match

        if screen:
            periodic = self.get_periodic_sites()
            symmetry_match = symmetry_match[periodic]
        if unique:
            return np.unique(symmetry_match)

        else:
            return symmetry_match

    def get_adsorption_vectors(self, unique=True, screen=True):
        """Returns the vectors representing the furthest distance from the neighboring atoms.

        Returns
        -------
        vectors : ndarray (n, 3)
            Adsorption vectors for surface sites.
        """
        top_coords = self.coordinates[self.connectivity == 1]
        if unique:
            sel = self.get_symmetric_sites(screen=screen)
        else:
            sel = self.get_periodic_sites(screen=screen)
        coords = self.coordinates[sel]
        r1top = self.r1_topology[sel]
        r2top = self.r2_topology[sel]

        vectors = np.empty((coords.shape[0], 3))
        for i, s in enumerate(coords):
            plane_points = np.array(list(r1top[i]) + list(r2top[i]), dtype=int)
            vectors[i] = utils.plane_normal(top_coords[plane_points])

        return vectors

    def get_adsorption_edges(self, symmetric=True, periodic=True):
        """Return the edges of adsorption sties defined as all regions with adjacent vertices.

        Parameters
        ----------
        symmetric : bool
            Return only the symmetrically reduced edges.
        periodic : bool
            Return edges which are unique via periodicity.

        Returns
        -------
        edges : ndarray (n, 2)
            All edges crossing ridge or vertices indexed by the expanded
            unit slab.
        """
        vt = scipy.spatial.Voronoi(self.coordinates[:, :2], qhull_options=f"Qbb Qc Qz C{1e-2}")

        select, lens = [], []
        for i, p in enumerate(vt.point_region):
            select += [vt.regions[p]]
            lens += [len(vt.regions[p])]

        dmax = max(lens)
        regions = np.zeros((len(select), dmax), int)
        mask = np.arange(dmax) < np.array(lens)[:, None]
        regions[mask] = np.concatenate(select)

        site_id = self.get_symmetric_sites(unique=False, screen=False)
        site_id = site_id + self.connectivity / 10
        per = self.get_periodic_sites(screen=False)

        uper = self.get_periodic_sites()
        edges, symmetry, uniques = [], [], []
        for i, p in enumerate(uper):
            poi = vt.point_region[p]
            voi = vt.regions[poi]

            for v in voi:
                nr = np.where(regions == v)[0]

                for n in nr:
                    edge = sorted((p, n))

                    if n in uper[: i + 1] or edge in edges:
                        continue

                    if (np.in1d(per[edge], per[uper[:i]]).any()) and periodic:
                        continue

                    sym = sorted(site_id[edge])
                    if sym in symmetry:
                        uniques += [False]
                    else:
                        uniques += [True]
                        symmetry += [sym]

                    edges += [edge]

        edges = np.array(edges)
        if symmetric:
            edges = edges[uniques]

        return edges

    def ex_sites(self, index, select="inner", cutoff=0):
        """Get site indices inside or outside of a cutoff radii from a provided periodic site
        index.

        If two sites are provided, an option to return the mutually inclusive points is also
        available.
        """
        per = self.get_periodic_sites(False)
        sym = self.get_symmetric_sites()
        edges = self.get_adsorption_edges(symmetric=False, periodic=False)
        coords = self.coordinates[:, :2]

        if isinstance(index, int):
            index = [index]

        if not cutoff:
            for i in per[index]:
                sd = np.where(edges == i)[0]
                select_coords = coords[edges[sd]]
                d = np.linalg.norm(np.diff(select_coords, axis=1), axis=2)
                cutoff = max(d.max(), cutoff)
        cutoff += self.tol

        diff = coords[:, None] - coords[sym]
        norm = np.linalg.norm(diff, axis=2)
        neighbors = np.array(np.where(norm < cutoff))

        neighbors = []
        for i in index:
            diff = coords[:, None] - coords[per[i]]
            norm = np.linalg.norm(diff, axis=2)
            if select == "mutual" and len(index) == 2:
                neighbors += [np.where(norm < cutoff)[0].tolist()]
            else:
                neighbors += np.where(norm < cutoff)[0].tolist()

        if select == "inner":
            return per[neighbors]
        elif select == "outer":
            return np.setdiff1d(per, per[neighbors])
        elif select == "mutual":
            return np.intersect1d(per[neighbors[0]], per[neighbors[1]])

    def plot(self, savefile=None):
        """Create a visualization of the sites."""
        top = self.connectivity == 1
        other = self.connectivity != 1
        dt = scipy.spatial.Delaunay(self.coordinates[:, :2][top])

        fig = plt.figure(figsize=(6, 4), frameon=False)
        ax = fig.add_axes([0, 0, 1, 1])

        ax.triplot(dt.points[:, 0], dt.points[:, 1], dt.simplices.copy())
        ax.plot(dt.points[:, 0], dt.points[:, 1], "o")
        ax.plot(self.coordinates[:, 0][other], self.coordinates[:, 1][other], "o")
        ax.axis("off")

        if savefile:
            plt.savefig(savefile, transparent=True)
        else:
            plt.show()


class Builder(AdsorptionSites):
    """Initial module for creation of 3D slab structures with attached adsorbates."""

    def __repr__(self):
        formula = self.slab.get_chemical_formula()
        string = f"Adsorption builder for {formula} slab.\n"
        sym = len(self.get_symmetric_sites())
        string += f"unique adsorption sites: {sym}\n"
        con = self.get_connectivity()
        string += f"site connectivity: {con}\n"
        edges = self.get_adsorption_edges()
        string += f"unique adsorption edges: {len(edges)}"

        return string

    def add_adsorbates(self, adsorbates, indices):
        """"""
        slab = self.slab.copy()
        for i, ads in enumerate(adsorbates):
            bond = np.where(ads.get_tags() == -1)[0][0]

            slab = self._single_adsorption(
                ads, slab=slab, bond=bond, site_index=indices[i], symmetric=False
            )

        return slab

    def add_adsorbate(self, adsorbate, bonds=None, index=0, auto_construct=True, **kwargs):
        """Add and adsorbate to a slab. If the auto_constructor flag is False, the atoms object
        provided will be attached at the active site.

        Parameters
        ----------
        adsorbate : gratoms object
            Molecule to connect to the surface.
        bonds : int or list of 2 int
            Index of adsorbate atoms to be bonded.
        index : int
            Index of the site or edge to use as the adsorption position. A
            value of -1 will return all possible structures.
        auto_construct : bool
            Whether to automatically estimate the position of atoms in larger
            molecules or use the provided structure.

        Returns
        -------
        slabs : gratoms object
            Slab(s) with adsorbate attached.
        """
        if bonds is None:
            # Molecules with tag -1 are designated to bond
            bonds = np.where(adsorbate.get_tags() == -1)[0]

        if len(bonds) == 0:
            raise ValueError("Specify the index of atom to bond.")

        elif len(bonds) == 1:
            if index == -1:
                slab = []
                for i, _ in enumerate(self.get_symmetric_sites()):
                    slab += [
                        self._single_adsorption(
                            adsorbate,
                            bond=bonds[0],
                            site_index=i,
                            auto_construct=auto_construct,
                            **kwargs,
                        )
                    ]
            elif isinstance(index, (list, np.ndarray)):
                slab = []
                for i in index:
                    slab += [
                        self._single_adsorption(
                            adsorbate,
                            bond=bonds[0],
                            site_index=i,
                            auto_construct=auto_construct,
                            **kwargs,
                        )
                    ]
            else:
                slab = self._single_adsorption(
                    adsorbate,
                    bond=bonds[0],
                    site_index=index,
                    auto_construct=auto_construct,
                    **kwargs,
                )

        elif len(bonds) == 2:
            if index == -1:
                slab = []
                edges = self.get_adsorption_edges()
                for i, _ in enumerate(edges):
                    slab += [
                        self._double_adsorption(adsorbate, bonds=bonds, edge_index=i, **kwargs)
                    ]
            else:
                slab = self._double_adsorption(adsorbate, bonds=bonds, edge_index=index, **kwargs)

        else:
            raise ValueError("Only mono- and bidentate adsorption supported.")

        return slab
    #xzq_test
    @staticmethod
    def inertia_tensor(positions, masses):
        """计算分子的惯性矩张量"""
        centroid = np.average(positions, axis=0, weights=masses)  # 计算质心
        positions -= centroid  # 将坐标平移到质心
        tensor = np.zeros((3, 3))
        for i in range(len(positions)):
            r = positions[i]
            m = masses[i]
            tensor += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
        return tensor
    @staticmethod
    def get_principal_axes(inertia_tensor):
        """对惯性矩张量进行对角化，返回主惯性轴"""
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)
        sorted_indices = np.argsort(eigenvalues)[::-1]  # 按惯性矩大小排序（从大到小）
        principal_axes = eigenvectors[:, sorted_indices]
        return principal_axes
    
    def _single_adsorption(
        self,
        adsorbate,
        bond,
        slab=None,
        site_index=0,
        auto_construct=True,
        enable_rotate_xoy=True,
        rotation_mode='vnn',
        rotation_args={},
        direction_mode='default',
        direction_args={},  # 用于指定惯性矩模式
        site_coord=None,
        z_bias=0,
        symmetric=True,
        verbose=False
    ):
        """Bond and adsorbate by a single atom."""
        if slab is None:
            slab = self.slab.copy()
        atoms = adsorbate.copy()
        atoms.set_cell(slab.cell)

        if symmetric:
            ind = self.get_symmetric_sites()[site_index]
            vector = self.get_adsorption_vectors()[site_index]
        else:
            ind = self.get_periodic_sites()[site_index]
            vector = self.get_adsorption_vectors(unique=False)[site_index]

        # Improved position estimate for site.
        u = self.r1_topology[ind]
        r = radii[slab[self.index[u]].numbers]
        top_sites = self.coordinates[self.connectivity == 1]

        numbers = atoms.numbers[bond]
        R = radii[numbers]
        base_position = utils.trilaterate(top_sites[u], r + R, vector)

        # Zhaojie Wang   20230910(precise adsorption coord)
        if not site_coord is None:
            base_position = site_coord

        if auto_construct:
            if direction_mode == 'asphericity':
                # 根据分子的形状调整朝向，使其“平躺”在表面上
                masses = atoms.get_masses()
                positions = atoms.get_positions()
                #xzq
                # 计算惯性矩张量
                inertia_tensor_value = self.inertia_tensor(positions, masses)

                # 获取主惯性轴
                principal_axes = self.get_principal_axes(inertia_tensor_value)

                # 自动处理三种惯性矩方向
                slabs_list = []
                for inertia_mode in range(1, 4):  # 遍历第一、第二、第三惯性矩
                    atoms_copy = atoms.copy()  # ✅ 每次都从原始结构复制
                    base_position = [0.0, 0.0, 0.0]
                    adsorption_vector = principal_axes[:, inertia_mode - 1]                    
                    atoms_copy.rotate(adsorption_vector, [0, 0, 1])

                    if enable_rotate_xoy and rotation_mode == 'vnn' and rotation_args != {}:
                        principle_axe = utils.solve_principle_axe_pca(atoms_copy.get_positions())
                        if abs(rotation_args['vec_to_neigh_imgsite'][0]) < 1e-8:
                            target_vec = [1, 0, 0]
                        elif abs(rotation_args['vec_to_neigh_imgsite'][1]) < 1e-8:
                            target_vec = [0, 1, 0]
                        else:
                            target_vec = [-1/rotation_args['vec_to_neigh_imgsite'][0], 
                                          1/rotation_args['vec_to_neigh_imgsite'][1], 0]
                        atoms_copy.rotate([principle_axe[0], principle_axe[1], 0], target_vec)
                # =============== ✅ 关键部分：设置 Z 高度 ===============
                    z_coordinates = atoms_copy.get_positions()[:, 2]
                    min_z = np.min(z_coordinates)

                    final_positions = slab.get_positions()
                    max_z = np.max(final_positions[:, 2])

                    
                    if abs(site_coord[2]) < 1e-6:  # 如果为 0 或非常接近 0
                        effective_distance = 2.5
                    else:
                        effective_distance = site_coord[2]

                    base_position[2] = round(max_z - min_z + effective_distance, 2)

                   # xzq如果输入的 xy 为 [0, 0]，则将分子放在 slab 的 xy 中心
                    if site_coord[0] == 0 and site_coord[1] == 0:
                       center_x, center_y = utils.center_slab(final_positions)  # 计算 slab 的 xy 中心
                       base_position[0] = round(center_x, 1)  # 将分子放置在 xy 中心
                       base_position[1] = round(center_y, 1)
                    else:
                       # 否则，使用输入的 site_coord 作为 xy 坐标
                       base_position[0] = round(site_coord[0], 1)
                       base_position[1] = round(site_coord[1], 1)

                    atoms_copy.translate(base_position)  # 确保所有构型都以正确的 z 轴间距平移

                    # 将分子平移到指定的 site_coord 位置
                    vec_tmp = site_coord - atoms_copy.get_positions()[bond]  # 计算平移量
                    atoms_copy.translate([vec_tmp[0], vec_tmp[1], 0])  # 只平移 x-y，不改变 z


                    slab_with_adsorbate = slab.copy()
                    slab_with_adsorbate += atoms_copy
                    slabs_list.append(slab_with_adsorbate)

                return slabs_list  # 返回三种吸附构型的列表

            elif direction_mode == 'bond_atom':
                # 根据参与吸附的原子确定位向，将物种“扶正”
                adsorption_vector, flag = utils.solve_normal_vector_linearsvc(atoms.get_positions(), bond)
                atoms.rotate(adsorption_vector, [0, 0, 1])
            elif direction_mode == 'asphericity':
                # 根据分子的形状调整朝向，使其“平躺”在表面上
                masses = atoms.get_masses()
                positions = atoms.get_positions()
                #xzq
                # 计算惯性矩张量
                inertia_tensor_value = self.inertia_tensor(positions, masses)

                # 获取主惯性轴
                principal_axes = self.get_principal_axes(inertia_tensor_value)

                # 自动处理三种惯性矩方向
                slabs_list = []
                for inertia_mode in range(1, 4):  # 遍历第一、第二、第三惯性矩
                    atoms_copy = atoms.copy()  # ✅ 每次都从原始结构复制
                    base_position = [0.0, 0.0, 0.0]
                    adsorption_vector = principal_axes[:, inertia_mode - 1]                    
                    atoms_copy.rotate(adsorption_vector, [0, 0, 1])

                    if enable_rotate_xoy and rotation_mode == 'vnn' and rotation_args != {}:
                        principle_axe = utils.solve_principle_axe_pca(atoms_copy.get_positions())
                        if abs(rotation_args['vec_to_neigh_imgsite'][0]) < 1e-8:
                            target_vec = [1, 0, 0]
                        elif abs(rotation_args['vec_to_neigh_imgsite'][1]) < 1e-8:
                            target_vec = [0, 1, 0]
                        else:
                            target_vec = [-1/rotation_args['vec_to_neigh_imgsite'][0], 
                                          1/rotation_args['vec_to_neigh_imgsite'][1], 0]
                        atoms_copy.rotate([principle_axe[0], principle_axe[1], 0], target_vec)
                    # 设置 吸附 高度
                    z_coordinates = atoms_copy.get_positions()[:, 2]
                    min_z = np.min(z_coordinates)

                    final_positions = slab.get_positions()
                    max_z = np.max(final_positions[:, 2])

                    
                    if abs(site_coord[2]) < 1e-6:  # 如果为 0 或非常接近 0
                        effective_distance = 2
                    else:
                        effective_distance = site_coord[2]

                    base_position[2] = round(max_z - min_z + effective_distance, 2)

                    #使用输入的 site_coord 作为 xy 坐标
                    base_position[0] = round(site_coord[0], 1)
                    base_position[1] = round(site_coord[1], 1)

                    atoms_copy.translate(base_position)  # 确保所有构型都以正确的 z 轴间距平移

                    # 将分子平移到指定的 site_coord 位置
                    vec_tmp = site_coord - atoms_copy.get_positions()[bond]  # 计算平移量
                    atoms_copy.translate([vec_tmp[0], vec_tmp[1], 0])  # 只平移 x-y，不改变 z


                    slab_with_adsorbate = slab.copy()
                    slab_with_adsorbate += atoms_copy
                    slabs_list.append(slab_with_adsorbate)

                return slabs_list  # 返回三种吸附构型的列表)
            elif direction_mode == 'hetero': # zjwang 20240815
                # 适用于有明确“官能团”的偏链状的分子
                # 求出从分子杂原子中心到所有原子中心的吸附矢量vec_p，将vec_p旋转到[0,0,1]竖直向上
                # 加旋转：将此时的分子绕z轴分别旋转：
                # 加偏斜：将此时的[0,0,1]分别旋转至：
                # 共3*5=15个atoms对象
                # hetero模式下无条件以最靠近杂原子形心的原子id作为bond
                print('===== hetero group mode =====')
                atoms_list = [] # list of <atoms>
                vec_p, centroid, centroid_hetero = utils.solve_normal_vector_hetero(atoms.get_positions(), atoms.get_chemical_symbols()) # hetero模式最好不设置bond原子，无意义重复（等于平移）
                # 先确定bond原子id
                d_min = 100
                for idx_a,pos in enumerate(atoms.get_positions()):
                    d = np.sqrt( np.sum( np.dot(pos-centroid_hetero,pos-centroid_hetero) ) )
                    if d < d_min:
                        d_min = d
                        bond = idx_a
                print(' bond, symbol, d_min_to_centroid_hetero =', bond, atoms.get_chemical_symbols()[bond], d_min)
                # 再对atoms坐标进行操作
                atoms.rotate(vec_p, [0, 0, 1])
                for deg in [0, 72, 144, 216, 288]: # 先旋转再偏斜，与表面接触处的变化更多样
                    for vec_bias in [[0,0,1], [0.2,0,1], [0,0.2,1], [-0.2,0,1], [0,-0.2,1]]:
                        atoms_list.append(copy.deepcopy(atoms))
                        atoms_list[-1].rotate([0,0,1], vec_bias) # 吸附矢量偏斜
                        atoms_list[-1].rotate(deg, [0,0,1])# 绕z轴旋转
                print('*** len(atoms_list) =', len(atoms_list))
            else: # direction_mode == 'default':
                atoms.rotate([0, 0, 1], vector)
            # 再在xoy平面中旋转物种以避免重叠（思路：投影“长轴”与rotation_args['vec_to_neigh_imgsite']垂直）
            # “长轴”：物种原子在xoy平面上的投影坐标进行PCA降维（2to1），找降维后的主元，其方向即“长轴”方向
            if enable_rotate_xoy and rotation_mode == 'vertical to vec_to_neigh_imgsite' and rotation_args != {}:
                principle_axe = utils.solve_principle_axe_pca(atoms.get_positions())
                ### print('principle_axe:\n', principle_axe)
                if abs(rotation_args['vec_to_neigh_imgsite'][0]) < 1e-8:
                    target_vec = [1, 0, 0]
                elif abs(rotation_args['vec_to_neigh_imgsite'][1]) < 1e-8:
                    target_vec = [0, 1, 0]
                else:
                    target_vec = [-1/rotation_args['vec_to_neigh_imgsite'][0], 1/rotation_args['vec_to_neigh_imgsite'][1], 0]
                ### print('target_vec:\n', target_vec)
                atoms.rotate([principle_axe[0], principle_axe[1], 0], target_vec)
        
        
        # zjwang 20240815
        if direction_mode == 'hetero':
            max_z = np.max(slab.get_positions()[:,2]) #获取slabz轴最大值
            n = len(slab)
            slabs_list = []
            score_configurations = [] # 各个吸附构型(slab)的“分数”
            for ia, a_orig in enumerate(atoms_list):
                slabs_list.append(copy.deepcopy(slab))
                a = copy.deepcopy(a_orig)  # 确保每次处理原始未改动的 adsorbate
                dealt_positions = a.get_positions()
                min_z = np.min(dealt_positions[:,2]) #得到分子z轴最小值
                base_position[2] = max_z - min_z + 2 # 吸附物种最低原子处于slab以上2.2A，随后再尝试降低
                a.translate(base_position)

                # 将分子平移到指定的 site_coord 位置
                xy_translation = site_coord[:2] - a.get_positions()[bond][:2]
                a.translate([xy_translation[0], xy_translation[1], 0])  # 只移动 x-y，不改变 z
                print('Final base_position:', base_position)
    
                # 组合 slab 和吸附物 
                slabs_list[-1] += a
                # Add graph connections
                for metal_index in self.index[u]:
                    slabs_list[-1].graph.add_edge(metal_index, bond + n)
                '''# 评估当前构型并不断尝试将吸附物降低以尽可能贴近表面，直到score不再提高
                score_tmp = utils.score_configuration_hetero(coords=slabs_list[-1].get_positions(),
                                                             symbols=slabs_list[-1].get_chemical_symbols(),
                                                             z_surf=max_z)
                score_before = -100.0
                while score_tmp > score_before:
                    slabs_list[-1] = copy.deepcopy(slab)
                    a.translate([0, 0, -0.02])
                    slabs_list[-1] += a
                    for metal_index in self.index[u]:
                        slabs_list[-1].graph.add_edge(metal_index, bond + n)
                    score_before = score_tmp
                    score_tmp = utils.score_configuration_hetero(coords=slabs_list[-1].get_positions(),
                                                                 symbols=slabs_list[-1].get_chemical_symbols(),
                                                                 z_surf=max_z)
                # 复位：撤销最后一次下降操作
                slabs_list[-1] = copy.deepcopy(slab)
                a.translate([0, 0, -0.02])
                slabs_list[-1] += a
                for metal_index in self.index[u]:
                    slabs_list[-1].graph.add_edge(metal_index, bond + n)
                score_configurations.append(score_before)
                # show & write log (mainly for score_configurations)
                str_log = str(ia) + ' | score = ' + str(np.round(score_configurations[-1],3)).ljust(8) + '(x,y,z): ' + str(base_position)
                ### print(str_log)
                with open('score_log.txt', 'a') as f:
                    f.write(str_log+'\n')
                    '''
            with open('score_log.txt', 'a') as f:
                f.write('\n')
                f.write('Ranking configurations by their scores:\n')
                for i,idx in enumerate(np.argsort(score_configurations)[::-1]):
                    f.write(str(i).ljust(4)+':    '+str(idx).ljust(8)+str(score_configurations[idx])+'\n')
            
            return slabs_list
        
        #lbx
        if base_position[2] == 0.0:
            #z_coordinates = rotated_positions[:, 2]
            z_coordinates = atoms.get_positions()[:, 2]
            min_z = np.min(z_coordinates) #得到分子z轴最小值
            #print(rotated_positions)
            #print(min_z)
            final_positions = slab.get_positions() #slab坐标  
            z_coordinates = final_positions[:, 2]
            max_z = np.max(z_coordinates) #获取slabz轴最大值
            base_position[2] = round(0 - min_z + 2.5 + max_z,1)
            print(f"min_z (adsorbate): {min_z}, max_z (slab): {max_z}, new base_position[2]: {base_position[2]}")
            #计算slab中心坐标
            center_x, center_y = utils.center_slab(final_positions)
            #print("(x, y):", center_x,center_y,base_position[2])
            base_position[0] = round(center_x ,1)
            base_position[1] = round(center_y ,1)
            print("(x, y,z):", base_position[0],base_position[1],base_position[2])
            atoms.translate(base_position)
            #atoms.translate(base_position)
            n = len(slab)
            slab += atoms
            #lbx:slab为增加分子后的坐标,base_position[2]为config输入的z
        else:
            # print('base_position:', len(base_position), '\n', base_position)
            base_position[2] = base_position[2] + z_bias
            atoms.translate(base_position)
            n = len(slab)
            slab += atoms


        # Add graph connections
        for metal_index in self.index[u]:
            slab.graph.add_edge(metal_index, bond + n)

        return slab


    def _double_adsorption(self, adsorbate, bonds=None, edge_index=0):
        """Bond and adsorbate by two adjacent atoms."""
        slab = self.slab.copy()
        atoms = adsorbate.copy()
        atoms.set_cell(slab.cell)

        numbers = atoms.numbers[bonds]
        R = radii[numbers] * 0.95
        edges = self.get_adsorption_edges()
        coords = self.coordinates[edges[edge_index]]

        U = self.r1_topology[edges[edge_index]]
        for i, u in enumerate(U):
            r = radii[slab[self.index[u]].numbers] * 0.95
            top_sites = self.coordinates[self.connectivity == 1]
            coords[i] = utils.trilaterate(top_sites[u], R[i] + r)

        vec = coords[1] - coords[0]
        n = np.linalg.norm(vec)
        uvec0 = vec / n
        d = np.sum(radii[numbers]) * 0.95
        dn = (d - n) / 2

        base_position0 = coords[0] - uvec0 * dn
        base_position1 = coords[1] + uvec0 * dn

        # Position the base atoms
        atoms[bonds[0]].position = base_position0
        atoms[bonds[1]].position = base_position1

        # Temporarily break adsorbate bond
        atoms.graph.remove_edge(*bonds)

        vectors = self.get_adsorption_vectors(screen=False, unique=False)
        uvec1 = vectors[edges[edge_index]]
        uvec2 = np.cross(uvec1, uvec0)
        uvec2 /= -np.linalg.norm(uvec2, axis=1)[:, None]
        uvec1 = np.cross(uvec2, uvec0)

        """ ### zjwang 20230426
        branches0 = list(nx.bfs_successors(atoms.graph, bonds[0]))
        if len(branches0[0][1]) != 0:
            uvec = [-uvec0, uvec1[0], uvec2[0]]
            self._branch_bidentate(atoms, uvec, branches0[0])
            for branch in branches0[1:]:
                catkit.gen.molecules._branch_molecule(
                    atoms, branch, adsorption=True)

        branches1 = list(nx.bfs_successors(atoms.graph, bonds[1]))
        if len(branches1[0][1]) != 0:
            uvec = [uvec0, uvec1[0], uvec2[0]]
            self._branch_bidentate(atoms, uvec, branches1[0])
            for branch in branches1[1:]:
                catkit.gen.molecules._branch_molecule(
                    atoms, branch, adsorption=True)
        """

        n = len(slab)
        slab += atoms
        # Add graph connections
        slab.graph.add_edge(*np.array(bonds) + n)
        for i, u in enumerate(U):
            for metal_index in self.index[u]:
                slab.graph.add_edge(metal_index, bonds[i] + n)

        return slab

    def _branch_bidentate(self, atoms, uvec, branch):
        """Return extended positions for additional adsorbates based on provided unit vectors."""
        r, nodes = branch
        num = atoms.numbers[[r] + nodes]
        d = radii[num[1:]] + radii[num[0]]
        c = atoms[r].position

        # Single additional atom
        if len(nodes) == 1:
            coord0 = (
                c
                + d[0] * uvec[0] * np.cos(1 / 3.0 * np.pi)
                + d[0] * uvec[1] * np.sin(1 / 3.0 * np.pi)
            )
            atoms[nodes[0]].position = coord0

        # Two branch system
        elif len(nodes) == 2:
            coord0 = (
                c
                + d[0] * uvec[1] * np.cos(1 / 3.0 * np.pi)
                + 0.866 * d[0] * uvec[0] * np.cos(1 / 3.0 * np.pi)
                + 0.866 * d[0] * uvec[2] * np.sin(1 / 3.0 * np.pi)
            )
            atoms[nodes[0]].position = coord0

            coord1 = (
                c
                + d[1] * uvec[1] * np.cos(1 / 3.0 * np.pi)
                + 0.866 * d[1] * uvec[0] * np.cos(1 / 3.0 * np.pi)
                + 0.866 * d[1] * -uvec[2] * np.sin(1 / 3.0 * np.pi)
            )
            atoms[nodes[1]].position = coord1

        else:
            raise ValueError("Too many bonded atoms to position correctly.")


def get_adsorption_sites(slab, surface_atoms=None, symmetry_reduced=True, tol=1e-5):
    """Get the adsorption sites of a slab as defined by surface symmetries of the surface atoms.

    Parameters
    ----------
    slab : Atoms object
        The slab to find adsorption sites for. Must have surface
        atoms defined.
    surface_atoms : array_like (N,)
        List of the surface atom indices.
    symmetry_reduced : bool
        Return the symmetrically unique sites only.
    adsorption_vectors : bool
        Return the adsorption vectors.

    Returns
    -------
    coordinates : ndarray (N, 3)
        Cartesian coordinates of activate sites.
    connectivity : ndarray (N,)
        Number of bonds formed for a given adsorbate.
    symmetry_index : ndarray (N,)
        Arbitrary indices of symmetric sites.
    """
    if surface_atoms is not None:
        slab.set_surface_atoms(surface_atoms)

    sites = AdsorptionSites(slab)

    if symmetry_reduced:
        s = sites.get_symmetric_sites()
    else:
        s = sites.get_periodic_sites()

    coordinates = sites.coordinates[s]
    connectivity = sites.connectivity[s]

    if symmetry_reduced:
        return coordinates, connectivity

    symmetries = sites.get_symmetric_sites(unique=False)
    unique, counts = np.unique(symmetries, return_counts=True)
    symmetry_index = np.arange(len(unique)).repeat(counts)

    return coordinates, connectivity, symmetry_index


def _get_adsorption_sites(surface_positions, tol=1e-5):
    """Return the positions and topology of adsorption sites based on the provided 3D positions of
    surface atoms for a structure.

    This function is intended for internal use only.

    Parameters
    ----------
    surface_positions : ndarray (N, 3)
        Cartesian coordinates for the surface atoms of a structure.
    tol : float
        Float point precision tolerance.

    Returns
    -------
    positions : ndarray (M, 3)
        Cartesian coordinates for the identifed adsorption sites. This is
        a superset which contrains the surface atoms positions.
    r1topology : list of lists (M, X)
        First nearest-neighbors of each of the adsorption sites associated
        with the positions.
    r2topology : list of lists (M, Y)
        Second nearest-neighbors of each of the adsorption sites associated
        with the positions.
    """
    positions = [surface_positions, [], [], []]
    r1topology = [[[i] for i, _ in enumerate(surface_positions)], [], [], []]
    r2topology = [[[] for _ in surface_positions], [], [], []]

    dt = scipy.spatial.Delaunay(positions[0][:, :2])
    neighbors = dt.neighbors
    simplices = dt.simplices

    for i, corners in enumerate(simplices):
        cir = scipy.linalg.circulant(corners)
        edges = cir[:, 1:]

        # Inner angle of each triangle corner
        vec = positions[0][edges.T] - positions[0][corners]
        uvec = vec.T / np.linalg.norm(vec, axis=2).T
        angles = np.sum(uvec.T[0] * uvec.T[1], axis=1)

        # Angle types
        right = np.isclose(angles, 0)
        obtuse = angles < -tol

        rh_corner = corners[right]
        edge_neighbors = neighbors[i]

        bridge = np.sum(positions[0][edges], axis=1) / 2

        # Looping through corners allows for elimination of
        # redundant points, identification of 4-fold hollows,
        # and collection of bridge neighbors.
        for j, c in enumerate(corners):
            edge = sorted(edges[j])

            if edge in r1topology[1]:
                continue

            # Get the bridge neighbors (for adsorption vector)
            neighbor_simplex = simplices[edge_neighbors[j]]
            oc = list(set(neighbor_simplex) - set(edge))[0]

            # Right angles potentially indicate 4-fold hollow
            potential_hollow = edge + sorted([c, oc])
            if c in rh_corner:
                if potential_hollow in r1topology[3]:
                    continue

                # Assumption: If not 4-fold, this suggests
                # no hollow OR bridge site is present.
                ovec = positions[0][edge] - positions[0][oc]
                ouvec = ovec / np.linalg.norm(ovec)
                oangle = np.dot(*ouvec)
                oright = np.isclose(oangle, 0)
                if oright:
                    positions[3] += [bridge[j]]
                    r1topology[3] += [potential_hollow]
                    r2topology[3] += []
                    r2topology[0][c] += [oc]
            else:
                positions[1] += [bridge[j]]
                r1topology[1] += [edge]
                r2topology[1] += [[c, oc]]

            r2topology[0][edge[0]] += [edge[1]]
            r2topology[0][edge[1]] += [edge[0]]

        if not right.any() and not obtuse.any():
            hollow = np.average(positions[0][corners], axis=0)
            positions[2] += [hollow]
            r1topology[2] += [corners.tolist()]
            r2topology[2] += []

    # For collecting missed bridge neighbors
    for s in r1topology[3]:
        edges = itertools.product(s[:2], s[2:])
        for edge in edges:
            edge = sorted(edge)
            i = r1topology[1].index(edge)
            n, m = r1topology[1][i], r2topology[1][i]
            nn = list(set(s) - set(n + m))

            if len(nn) == 0:
                continue
            r2topology[1][i] += [nn[0]]

    positions = np.array([j for i in positions for j in i])
    r1topology = np.array([j for i in r1topology for j in i])
    r2topology = np.array([j for i in r2topology for j in i])

    return positions, r1topology, r2topology


def symmetry_equivalent_points(fractional_coordinates, atoms, tol=1e-5):
    """Return the symmetrically equivalent points from a list of provided fractional coordinates.

    Parameters
    ----------
    fractional_coordinates : ndarray (N ,3)
        Fractional coordinates to find symmetrical equivalence
        between.
    atoms : Atoms object
        Atoms object to use the unit cell, positions, and pbc of.
    tol : float
        Float point precision tolerance.

    Returns
    -------
    symmetry_match : ndarray (N,)
        Indices of fractional coordinates which are unique or
        matching.
    """
    sym = symmetry.Symmetry(atoms, tol=1e-2)
    affine_matrices = sym.get_symmetry_operations()
    affine_matrices = np.swapaxes(affine_matrices, 1, 2)

    pbc = ~np.array(atoms.pbc.tolist() + [True])
    nt = np.isclose(affine_matrices[:, pbc, 3], 0).any(axis=1)
    nt &= np.isclose(affine_matrices[:, pbc, 2], 1).any(axis=1)
    affine_matrices = affine_matrices[nt]

    affine_points = np.insert(fractional_coordinates, 3, 1, axis=1)
    operations = np.dot(affine_points, affine_matrices)[:, :, :3]

    periodic_match = np.arange(fractional_coordinates.shape[0])
    symmetry_match = periodic_match.copy()
    for i, j in enumerate(symmetry_match):
        if i != j:
            continue

        d = operations[i, :, None] - fractional_coordinates
        d -= np.round(d)
        dind = np.where((np.abs(d) < tol).all(axis=2))[-1]
        symmetry_match[np.unique(dind)] = periodic_match[i]

    return symmetry_match
