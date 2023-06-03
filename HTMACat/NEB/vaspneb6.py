#!/usr/bin/python

###############################################################################
### This python scripts is used to generate the images for neb calcualtions ###
###    1.it will generate six images                                        ###
###    2.you should prepare three file: POSstart POSend and Nmoveatom       ###
###    3.please write the number of moving atoms in Nmoveatom file          ###
###        1) like 1 2 3 5                                                  ###
###        2) you can write 0(-1) in Nmoveatom make all atoms T(F)          ###
###        3) if -2 is in Nmoveatom, the selected atom will be fixed        ###
###                         Usage: ./vaspneb6.py                            ###
###              Written by xiaoliu (xiaoliu@hust.edu.cn)                   ###
###                      Last revised on Aug 17 2013                        ###
###############################################################################

###########neb-first################

import os


def vaspneb(num, dcut):
    """Linear interpolate n images with fixing the atoms within the dcut of distance change."""


dirnametmp = ["00"]
dirname = "".join(dirnametmp)
# print dirname
# os.chdir(dirname)
os.system("ls > dirnamefile")

dir = open("dirnamefile")
firstline = dir.readline()
# print firstline
if dirname in firstline:
    # os.system('echo 1')
    os.system("rm -rf 0* con")
    os.system("mkdir 00 01 02 03 04 05 06 07")
    os.system("echo it is not a first neb calculation")
    os.system("echo ")
else:
    # os.system('echo 2')
    os.system("mkdir 00 01 02 03 04 05 06 07")
    os.system("echo it is a first neb calculation")
    os.system("echo ")

os.system("rm -f dirnamefile")

###########neb-boundray################

# import os
POSstart = open("POSstart")
POSend = open("POSend")
POSstarttmp = open("POSstarttmp", "w")
POSendtmp = open("POSendtmp", "w")
# pos = open('pos','w')

# handle the head of POSstart and POSend
# handle POSstart
tmp = POSstart.readline()
POSstarttmp.write("System = NEB image \n")

for i in range(5):
    head = POSstart.readline()
    POSstarttmp.write(head)

# tmp = POSstart.readline()

head = POSstart.readline()
POSstarttmp.write(head)
Nelement = head.split()
Tatom = 0
for i in range(len(Nelement)):
    Tatom = Tatom + int(Nelement[i])

head = POSstart.readline()
POSstarttmp.write(head)

# handle POSend
tmp = POSend.readline()
POSendtmp.write("System = NEB image \n")

for i in range(5):
    head = POSend.readline()
    POSendtmp.write(head)

# tmp = POSend.readline()

head = POSend.readline()
POSendtmp.write(head)
Nelement = head.split()
Tatom = 0
for i in range(len(Nelement)):
    Tatom = Tatom + int(Nelement[i])

head = POSend.readline()
POSendtmp.write(head)

# handle the coordinate of POSstart and POSend
for j in range(Tatom):
    Coor_start = POSstart.readline().split()
    Coor_end = POSend.readline().split()
    for i in range(3):
        # print Coor_start[i]
        # print Coor_end[i]
        # print float(Coor_start[i])
        if abs(float(Coor_start[i]) - float(Coor_end[i])) > 0.95:
            #        print j+1
            POSstarttmp.write(Coor_start[i])
            POSstarttmp.write("  ")
            POSendtmp.write(Coor_start[i])
            POSendtmp.write("  ")
        else:
            POSstarttmp.write(Coor_start[i])
            POSstarttmp.write("  ")
            POSendtmp.write(Coor_end[i])
            POSendtmp.write("  ")
    POSstarttmp.write("\n")
    POSendtmp.write("\n")

POSstart.close()
POSend.close()
POSstarttmp.close()
POSendtmp.close()
# os.system('cat POSstarttmp POSendtmp > pos')
# pos.close()
# os.system('rm -f POSstarttmp POSendtmp')
os.system("echo all atoms have been moved into an unit cell")

###########neb-interpolate###############

# suppose we have 8 structures
# import os
# os.system('mkdir 00 01 02 03 04 05 06 07')
rootdir = os.getcwd()

for i in range(8):
    itmp = i
    dirnametmp = ["0", str(i)]
    dirname = "".join(dirnametmp)
    # print dirname
    os.chdir(dirname)
    os.system("cp ../POSstarttmp .")
    os.system("cp ../POSendtmp .")
    POSstart = open("POSstarttmp")
    POSend = open("POSendtmp")
    POSCAR = open("POSCAR", "w")

    # handle the head of POSCAR
    for i in range(8):
        head = POSend.readline()

    for i in range(6):
        head = POSstart.readline()
        POSCAR.write(head)

    head = POSstart.readline()
    POSCAR.write(head)
    Nelement = head.split()
    Tatom = 0
    for i in range(len(Nelement)):
        Tatom = Tatom + int(Nelement[i])

    head = POSstart.readline()
    POSCAR.write(head)

    # handle the coordinate of POSCAR
    for i in range(Tatom):
        Coordinatestart = POSstart.readline().split()
        Coordinateend = POSend.readline().split()
        Coordinatetmp = Coordinatestart
        # print Coordinatetmp[0]
        # print Coordinatestart[0]
        for i in range(3):
            # print Coordinatestart[i]
            # print Coordinateend[i]
            # print float(Coor_start[i])
            Coordinatetmp[i] = (
                float(Coordinatetmp[i])
                + float(itmp) * (float(Coordinateend[i]) - float(Coordinatestart[i])) / 7
            )
            # print Coordinatetmp[i]
            POSCAR.write(str(round(Coordinatetmp[i], 6)).ljust(10))
            POSCAR.write("  ")
        POSCAR.write("\n")

    POSCAR.close()
    POSstart.close()
    POSend.close()
    os.system("rm POSstarttmp POSendtmp")
    os.chdir(rootdir)

os.system("rm -f POSstarttmp POSendtmp")
os.system("echo ")
os.system("echo six images have been generated by linear interpolation")
os.system("echo ")

###########neb-selective################

# please write the number of moving atoms in Nmoveatom file
# like 1 2 3 5
# you can write 0(-1) in Nmoveatom make all atoms T(F)
# suppose we have 10 structures
# import os
rootdir = os.getcwd()
os.system("mkdir con")
# Nmoveatom = open('Nmoveatom', 'r')

for i in range(8):
    dirnametmp = ["0", str(i)]
    dirname = "".join(dirnametmp)
    # print dirname
    os.chdir(dirname)
    os.system("cp ../Nmoveatom .")
    POSCAR = open("POSCAR")
    POSfinal = open("POStmp", "w")
    Nmoveatom = open("Nmoveatom")
    # POSfinal.write('debug')
    # os.chdir(rootdir)

    # handle the head of POSCAR
    for i in range(6):
        head = POSCAR.readline()
        POSfinal.write(head)

    head = POSCAR.readline()
    POSfinal.write(head)
    Nelement = head.split()
    Tatom = 0
    for i in range(len(Nelement)):
        Tatom = Tatom + int(Nelement[i])

    POSfinal.write("Selective Dynamic \n")
    head = POSCAR.readline()
    POSfinal.write(head)
    # read the number of moving atoms
    Nmatom = Nmoveatom.readline().split()
    Num = []
    for i in range(len(Nmatom)):
        Num = Num + [int(Nmatom[i])]

    TMP = 1
    # handle the coordinate of POSCAR
    for i in range(Tatom):
        Coordinate = POSCAR.readline()
        # print  Coordinate.strip()
        if -2 in Num:
            TMP = -2
            if 0 in Num:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                POSfinal.write("   T   T   T \n")
                TMP = 0
            elif -1 in Num:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                POSfinal.write("   F   F   F \n")
                TMP = -1
            elif i + 1 in Num:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                # POSfinal.write('   T   T   T \n')
                POSfinal.write("   F   F   F \n")
            else:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                # POSfinal.write('   F   F   F \n')
                POSfinal.write("   T   T   T \n")
        else:
            if 0 in Num:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                POSfinal.write("   T   T   T \n")
                TMP = 0
            elif -1 in Num:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                POSfinal.write("   F   F   F \n")
                TMP = -1
            elif i + 1 in Num:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                POSfinal.write("   T   T   T \n")
            else:
                POSfinal.write("   ")
                POSfinal.write(Coordinate.strip().ljust(34))
                POSfinal.write("   F   F   F \n")

    POSCAR.close()
    POSfinal.close()
    Nmoveatom.close()
    os.system("rm Nmoveatom")
    os.system("mv POStmp POSCAR")
    os.chdir(rootdir)

os.system("cp 00/POSCAR con/POSCAR0")
os.system("cp 01/POSCAR con/POSCAR1")
os.system("cp 02/POSCAR con/POSCAR2")
os.system("cp 03/POSCAR con/POSCAR3")
os.system("cp 04/POSCAR con/POSCAR4")
os.system("cp 05/POSCAR con/POSCAR5")
os.system("cp 06/POSCAR con/POSCAR6")
os.system("cp 07/POSCAR con/POSCAR7")
# os.system('cp 08/POSCAR con/POSCAR8')
# os.system('cp 09/POSCAR con/POSCAR9')

if TMP == -1:
    os.system("echo all atoms have been set to fixed")
elif TMP == 0:
    os.system("echo all atoms have been set to mobile")
elif TMP == -2:
    os.system("echo only the selective atoms have been set to fixed")
else:
    os.system("echo only the selective atoms have been set to mobile")

os.system("echo ")
os.system("echo All settings have done! Congratulations!")
os.system("echo ")
