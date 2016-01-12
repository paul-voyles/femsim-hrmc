import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class EAM(object):
    def __init__(self, file):
        if '.eam.fs' in file:
            self.read_eam_fs(file)

    def _f(self, i, r):
        bin = int(r/self.dr)
        return self.f[i][bin]

    def _rho(self, i, j, r):
        bin = int(r/self.dr)
        return self.rho[i][j][bin]

    def _phi(self, i, j, r):
        bin = int(r/self.dr)
        return self.phi[i][j][bin]

    def read_eam_fs(self, file):
        content = open(file).readlines()
        n = 3
        nelements = int(content[n].split()[0])
        atom_syms = []
        for i in range(nelements):
            atom_syms.append(content[n].split()[i+1])
        n = 4
        nrho, drho, nr, dr, eam_max_r = (float(x) for x in content[n].split())
        nrho = int(nrho)
        nr = int(nr)
        self.nrho = nrho
        self.drho = drho
        self.nr = nr
        self.dr = dr
        self.eam_max_r = eam_max_r
        self.nelements = nelements

        self.znum = [0 for i in range(nelements)]
        self.mass = [0 for i in range(nelements)]
        self.latt_const = [0 for i in range(nelements)]
        self.f = [[0. for i in range(nrho)] for j in range(nelements)]
        self.rho = [[[0. for i in range(nr)] for j in range(nelements)] for k in range(nelements)]
        self.phi = [[[0. for i in range(nr)] for j in range(nelements)] for k in range(nelements)]

        n = 5
        # For each element
        for i in range(nelements):
            self.znum[i] = content[n].split()[0]
            self.mass[i] = content[n].split()[1]
            self.latt_const[i] = content[n].split()[2]
            n += 1
            # Read nrho values into f [0-10k], [30k-40k]
            for k in range(nrho):
                self.f[i][k] = float(content[n])
                n += 1
            # For each element-pair (*not* unique)
            for j in range(nelements):
                # Read nr values into rho [10k-20k], [20k-30k], [40k-50k], [50k-60k]
                for k in range(nr):
                    self.rho[i][j][k] = float(content[n])
                    n += 1

        # For each element
        for i in range(nelements):
            # For each *unique* element-pair
            for j in range(nelements):
                if i >= j:
                    # Read nr values into phi [50k-60k], [60k-70k], [70k-80k]
                    for k in range(nr):
                        x = float(content[n])
                        self.phi[i][j][k] = x#/((k+1)*dr)
                        self.phi[j][i][k] = x#/((k+1)*dr)
                        n += 1
        print(n)


    def read_eam_alloy(self, file):
        pass

def main():
    potentialfile = sys.argv[1]
    eam = EAM(potentialfile)
    print(eam._phi(1, 1, 4.0))
    print(eam._phi(1, 1, 5.0))

    xx1 = []
    yy1 = []
    xx2 = []
    yy2 = []
    xx3 = []
    yy3 = []
    with open('temp.txt', 'w') as f:
        xmin, ymin = 0, 0
        for i in range(eam.nr):
            r = i*eam.dr
            if r >= eam.eam_max_r - eam.dr:
                break
            y = eam._phi(0, 0, r)
            #y = eam._rho(0, 0, r)
            f.write("{0}  {1}\n".format(r, y))
            if y <= 2.5:
                xx1.append(r)
                yy1.append(y)
            if y < ymin:
                ymin = y
                xmin = r
            y = eam._phi(0, 1, r)
            #y = eam._rho(0, 1, r)
            if y <= 2.5:
                xx2.append(r)
                yy2.append(y)
            y = eam._phi(1, 1, r)
            #y = eam._rho(1, 1, r)
            if y <= 2.5:
                xx3.append(r)
                yy3.append(y)
        print(xmin, ymin)
    print(xx1[0], xx2[0], xx3[0])
    #assert eam.rho[1][1] == eam.rho[0][1]
    plt.plot(xx1, yy1)
    plt.plot(xx2, yy2)
    plt.plot(xx3, yy3)
    zero = [0.0 for i in range(len(xx1))]
    plt.plot(xx1, zero)
    plt.show()


if __name__ == '__main__':
    main()
