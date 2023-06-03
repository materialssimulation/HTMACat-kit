#!/usr/bin/python
# calculate the atom moving distance between initial and final structure for NEB calculation.
#     1.Usage : Dismove.py dcut
#     2.dcut : when atom moving distance > dcut, the atom id should be exported into Nmoveatom used to NEB calculation.
#     3.the distance in fraction coord  is exported into POStmp file
#     4.the distance in real xoord is exported into POStmp2 file
#     5.the diatance > dcut is exported into POSfinal file
#     6.the POStmp POStmp2 POSfinal are put into dismove dir
#     7.the POSstart & POSend are used in vasp coord file (tansforming cif into vasp in vesta)
# writen by JqYang in 20160512
#####################################################################################################################
import os
import sys
import linecache
import math


def Dismove(dcut):
    """calculate the atom moving distance between initial and final structure for NEB calculation.

    Parameters
    ----------
    dcut:float
        when atom moving distance > dcut, the atom id should be exported into Nmoveatom used to NEB calculation.

    Returns
    -------
    None

    Notes
    -----
    This function needs to read two files named POSstart and POSend, and write the processing results to three files:
    POStmp, POStmp2, and POSfinal.At the same time, this function will also create a file called Nmovetatom to record
    the sequence numbers of atoms whose movement distance is greater than the threshold. It is necessary to ensure that
    there is no folder named remove in the current working directory before calling this function, otherwise the folder
    will be deleted.
    """
    # dcut = float(sys.argv[1])
    POSstart = open("POSstart")
    POSend = open("POSend")
    POStmp = open("POStmp", "w")
    POStmp2 = open("POStmp2", "w")
    for i in range(2):
        head = POSstart.readline()
        head2 = POSend.readline()
        POStmp.write(head)
        POStmp2.write(head)

    # to x lattice parameter
    head = POSstart.readline()
    head2 = POSend.readline()
    xcoord = head.split()
    xcoord = float(xcoord[0])
    POStmp.write(head)
    POStmp2.write(head)

    # to y lattice parameter
    head = POSstart.readline()
    head2 = POSend.readline()
    ycoord = head.split()
    ycoord = float(ycoord[1])
    POStmp.write(head)
    POStmp2.write(head)

    # to z lattice parameter
    head = POSstart.readline()
    head2 = POSend.readline()
    zcoord = head.split()
    zcoord = float(zcoord[2])
    POStmp.write(head)
    POStmp2.write(head)

    head = POSstart.readline()
    head2 = POSend.readline()
    POStmp.write(head)
    POStmp2.write(head)

    head = POSstart.readline()
    head2 = POSend.readline()
    POStmp.write(head)
    POStmp2.write(head)

    # to the number of atoms
    Nelement = head.split()
    Tatom = 0
    for i in range(len(Nelement)):
        Tatom = Tatom + int(Nelement[i])

    head = POSstart.readline()
    head2 = POSend.readline()
    POStmp.write(head)
    POStmp2.write(head)

    POSftmp = open("POSfinal", "w")
    Nmovefile = open("Nmoveatom", "w")
    POStmp.write("Natom" + "\t\t\t\t " + "x" + "\t\t\t\t " + "y" + "\t\t\t\t " + "z" + "\n")
    POStmp2.write(
        "Natom"
        + "\t\t\t\t "
        + "x"
        + "\t\t\t\t "
        + "y"
        + "\t\t\t\t "
        + "z"
        + "\t\t\t\t"
        + "Dmove"
        + "\n"
    )
    POSftmp.write("dcut=" + str(dcut) + "\n")
    POSftmp.write("Natom" + "\t\t\t\t" + "Dmove-final" + "\n")
    # to atom move distance
    for i in range(int(Tatom)):
        head = POSstart.readline()
        head2 = POSend.readline()
        #     print head2
        head = head.split()
        head2 = head2.split()
        x = float(head2[0]) - float(head[0])
        if x > 0.5:
            x = x - 1
        if x < -0.5:
            x = x + 1
        y = float(head2[1]) - float(head[1])
        if y > 0.5:
            y = y - 1
        if y < -0.5:
            y = y + 1
        z = float(head2[2]) - float(head[2])
        xvalue = float(x) * float(xcoord)
        yvalue = float(y) * float(ycoord)
        zvalue = float(z) * float(zcoord)
        Dmove = math.sqrt(
            float(xvalue) * float(xvalue)
            + float(yvalue) * float(yvalue)
            + float(zvalue) * float(zvalue)
        )
        i = i + 1
        if Dmove > dcut:
            POSftmp.write("\t\t" + str(i) + "\t\t " + str(Dmove) + "\n")
            Nmovefile.write(str(i) + "\t")

        #     print  xvalue
        #     print x
        #     print y
        #     print z
        if x != 0 and y != 0 and z != 0:
            POStmp.write(
                "\t\t" + str(i) + "\t\t " + str(x) + "\t\t " + str(y) + "\t\t " + str(z) + " \n"
            )
            POStmp2.write(
                "\t\t"
                + str(i)
                + "\t\t "
                + str(xvalue)
                + "\t\t "
                + str(yvalue)
                + "\t\t "
                + str(zvalue)
                + "\t\t"
                + str(Dmove)
                + " \n"
            )
    POSstart.close()
    POSend.close()
    Nmovefile.close()
    POStmp.close()
    POStmp2.close()
    POSftmp.close()
    os.system("rm -rf dismove")
    os.system("mkdir dismove")
    # os.system('pwd')
    os.system("cp POStmp* POSf* dismove")
    os.system("rm POStmp* POSf*")
