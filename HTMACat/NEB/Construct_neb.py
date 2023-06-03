from ase import io
from ase.geometry import get_distances, Cell
from ase.neb import NEB
from ase.neb import interpolate
from ase.constraints import FixAtoms
from ase.visualize import view
import numpy as np


def Construct_neb(struct, method_inter="linear", nimages=6, dcut=0.25, d_scaled=0.5):
    """Construct a Nudged Elastic Band (NEB) for transitioning between two structures.

    Parameters
    ----------
    struct : list
        List of two strings representing the initial and final states in VASP file format.
    method_inter : str, optional
        Method for interpolating the positions of the middle images. Either 'linear' or 'idpp'. Default is 'linear'.
    nimages : int, optional
        Number of images in the band, excluding the initial and final states. Default is 6.
    dcut : float, optional
        Distance cutoff (in Angstroms) used to determine which atoms need to be fixed during the NEB calculation. Default is 0.25.
    d_scaled : float, optional
        Distance (scaled by the cell lengths) used to modify the positions of atoms that are too far apart between the initial and final structures. Default is 0.5.

    Returns
    -------
    images : list
        List of images representing the NEB calculation. The first and last images are the initial and final states, respectively, and the middle images are interpolated.

    Examples
    --------
    >>> struct = ['initial.vasp', 'final.vasp']
    >>> images = Construct_neb(struct, method_inter='idpp', nimages=10, dcut=0.2, d_scaled=0.3)s
    """
    ###Read initial and final states:
    initial = io.read(struct[0], format="vasp")
    final = io.read(struct[1], format="vasp")

    ###fixed atoms distance and indices,modify the positions
    p1_scaled = initial.get_scaled_positions()
    p2_scaled = final.get_scaled_positions()
    p1 = initial.get_positions()
    p2 = final.get_positions()

    # print(p2[-1],final.get_scaled_positions()[-1])
    length = initial.get_cell_lengths_and_angles()[0:3]
    Mask = []
    Natom = []
    if len(p1) == len(p2):
        for i in range(len(p1)):
            # get the distance of different axises
            # check the distance of every axis is or not larger than 0.5*cell length
            # If larger, add the length to modify the position
            for j in range(3):
                if (p2_scaled[i][j] - p1_scaled[i][j]) > d_scaled:
                    if abs(p2_scaled[i][j] - p1_scaled[i][j]) > abs(
                        p2_scaled[i][j] - (p1_scaled[i][j] + 1)
                    ):
                        p1_scaled[i][j] = p1_scaled[i][j] + 1
                        initial.set_scaled_positions(p1_scaled)
                    else:
                        pass
                elif (p2_scaled[i][j] - p1_scaled[i][j]) < -d_scaled:
                    if abs(p2_scaled[i][j] - p1_scaled[i][j]) > abs(
                        p2_scaled[i][j] + 1 - p1_scaled[i][j]
                    ):
                        p2_scaled[i][j] = p2_scaled[i][j] + 1
                        final.set_scaled_positions(p2_scaled)
                    else:
                        pass
                else:
                    continue

            # get the Euclidean Distance
            # vector1,vector2=p1[i],p2[i];dis = np.sqrt(np.sum(np.square(vector1-vector2)))
            dis = get_distances(p1[i], p2[i])[1][0][0]
            # p1_t=np.array([p1_scaled[i][0]*length[0],p1_scaled[i][1]*length[1],p1_scaled[i][2]*length[2]])
            # p2_t=np.array([p2_scaled[i][0]*length[0],p2_scaled[i][1]*length[1],p2_scaled[i][2]*length[2]])
            # dis_t= np.sqrt(np.sum(np.square(p1_t-p2_t)))
            # print(dis,dis_t)
            # according to the dcut(unit:Angstrom),choose the atoms needed to be fixed
            if dis > dcut:
                Mask += [False]
            else:
                Mask += [True]
                Natom += [i]

    else:
        print("The numbers of atoms between two structure are not equal !!!")
    # view(final)
    ###Make a band consisting of 6 images:
    images = [initial]
    images += [initial.copy() for i in range(nimages)]
    images += [final]
    neb = NEB(images)
    ###Interpolate idpp or linear the potisions of the six middle images:
    if method_inter == "linear":
        neb.interpolate(method="linear")
    elif method_inter == "idpp":
        neb.interpolate(method="idpp")
    else:
        print("the interpolate method no exist!")
        # ima=neb.iterimages()

    ###fixed atoms
    ##method 1:constraint = FixAtoms(indices=Natom)
    ##method 2:constraint = FixAtoms(mask=Mask)
    constraint = FixAtoms(mask=Mask)
    for image in images:
        # image.calc = Vasp(xc='PBE')
        image.set_constraint(constraint)
    return images


if __name__ == "__main__":
    import os
    from ase.io.vasp import write_vasp, read_vasp

    os.system("rm -rf 0*")
    ima = Construct_neb(["POSstart", "POSend"], method_inter="idpp", d_scaled=0.75, dcut=0.25)
    for i in range(len(ima)):
        os.mkdir(f"0{i}", 0o777)
        write_vasp("0%s/POSCAR" % (i), ima[i], direct=True, sort=[""], vasp5=True)
    print(f"{int(i)-1} images are constructed")
