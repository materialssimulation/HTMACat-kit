from ase.io import read
from catkit.gen.utils import to_gratoms
from collections import Counter
import os


def mig_atoms(Ffile, dis=[0.50, 0.50, 0.00]):
    """Migrate atoms in VASP POSCAR files.

    Parameters
    ----------
    Ffile : str
        A string that specifies the filename of the VASP POSCAR file.
    dis : list of float, optional
        A list of three floats that specifies the displacement of each atom
        along the x-, y-, and z-axis, respectively. Default is [0.50, 0.50, 0.00].

    Returns
    -------
    None

    Notes
    -----
    This function reads in a VASP POSCAR file, applies a displacement to each
    atom in the file, and writes out a new VASP POSCAR file with the migrated
    atoms. The displacement is specified by the `dis` parameter.

    Examples
    --------
    >>> mig_atoms('POSCAR', dis=[0.20, 0.20, 0.20])
    """
    # print(Ffile)
    for i, poscar in enumerate(list(Ffile)):
        poscar_bk = "".join([poscar, "_bk"])
        os.system(f"mv {poscar} {poscar_bk}")

        struc = read(poscar_bk, format="vasp")
        struc2 = to_gratoms(struc)
        symbols = struc2.get_chemical_symbols()
        dic_sym = dict(Counter(symbols))
        Num = sum(list(dic_sym.values()))
        constraint = struc2.constraints
        if constraint:
            Nline = 9
        else:
            Nline = 8

        pos = open(poscar_bk)
        # posfile = ''.join([poscar,'-new'])
        pos_new = open(poscar, "w")
        for j, line in enumerate(pos):
            if j < Nline:
                pos_new.write(line)
            elif j >= (int(Nline) + int(Num)):
                pos.close()
                pos_new.close()
                break
            else:
                coord = line.split().copy()
                Xtmp = float(coord[0]) + float(dis[0])
                Ytmp = float(coord[1]) + float(dis[1])
                Ztmp = float(coord[2]) + float(dis[2])
                if Xtmp > 1:
                    Xtmp = Xtmp - 1
                if Ytmp > 1:
                    Ytmp = Ytmp - 1
                if Ztmp > 1:
                    Ztmp = Ztmp - 1
                coord[0] = str(Xtmp)
                coord[1] = str(Ytmp)
                coord[2] = str(Ztmp)
                coordn = " ".join(coord)
                pos_new.write("%s\n" % (coordn))
        pos.close()
        pos_new.close()
    # os.system('mv POSstart POSstart_bk')
    # os.system('mv POSstart-new POSstart')
    # os.system('mv POSend POSend_bk')
    # os.system('mv POSend-new POSend')
    # return 0


# mig_atoms(['optmk/CONTCAR'])

# if __name__ == '__main__':
