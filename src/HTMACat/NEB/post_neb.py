# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import os


def calE_neb(nimage=6):
    ## read data and calculate barrier
    Ener_neb = open("react.log", "r+")
    neb_E = []
    neb_B = []
    Ei, Ef, Er = 0, 0, 0
    for i, line in enumerate(Ener_neb):
        #print(line)
        if i == 0:
            reaction = line.strip()
        elif i == 1:
            tmp = line.strip()
            Ei = tmp.split(',', 2)[0].split('=')[1].split()[0]
            Ef = tmp.split(',', 2)[1].split('=')[1].split()[0]
            Ent = tmp.split(',', 2)[2].split('=')[1].split()[0]
            neb_E += [Ei]
            neb_B += [i - 1, round(float(0), 2)]
        elif i < (2 + nimage):
            neb_E += [line.strip()]
            neb_B += [i - 1, round((float(neb_E[i - 1]) - float(neb_E[0])), 2)]
        else:
            i = i - 1
            break
    neb_E += [Ef]
    neb_B += [i, round(float(Ent), 2)]
    Ener_neb.close()
    return reaction, neb_B, neb_E

    ## plot the figure


def plt_neb(react, neb, show=True, name='neb'):
    X = neb[::2]
    Y = neb[1::2]
    barrier = max(Y)
    enthalpy = Y[-1]
    plt.figure()
    path = plt.plot(X, Y, 'go-')
    for x, y in zip(X, Y):
        plt.text(x, y, '%.2f' % y, ha='center', va='bottom', fontsize=11)
    plt.xlabel("Reaction Coordination")
    plt.ylabel("Reaction Energy(eV)")
    plt.title("Reaction Pathway")
    plt.legend(labels=[f"{react}:Ea={barrier}eV;Ef={enthalpy}eV"], loc=8)
    if show:
        plt.show()
        plt.savefig(f'{name}.png', dpi=300, bbox_inches='tight')
    return barrier, enthalpy


### post processing
if __name__ == '__main__':
    react, neb_B, neb_E = calE_neb()
    Ea, Ef = plt_neb(react, neb_B, show=True)
    print(react)
    print("Ea:", Ea)
    print("Ef:", Ef)
    Ener = open("react.log", "a+")
    Ener.write("-----Reaction Info-----\n")
    Ener.write(f"Reaction:{react}\n")
    Ener.write(f"{neb_E}\n")
    Ener.write(f"{neb_B}\n")
    Ener.write(f"Barrier:{Ea}\n")
    Ener.write(f"Enthalpy:{Ef}\n")
    Ener.write("----------------------\n")
    Ener.close()

    dir = os.getcwd()
    os.chdir("../")
    Nebfile = open("neb.csv", "a+")
    Nebfile.write(f"Reaction:{react}\n")
    Nebfile.write(f"Barrier:{Ea};Enthalpy:{Ef}\n")
    Nebfile.close()
    os.chdir(f"{dir}")
