import math

def calculate_rotations(nphi, npsi, ntheta, verbose=True):
    dphi = 2*math.pi / nphi
    dpsi = math.pi / npsi

    rot = [[0.0 for j in range(3)] for i in range(nphi*npsi*ntheta)]

    num_rot = 1  # Start at 1 so the initial entires are set to 0,0,0
    for i in range(0, nphi):
        phi = i * dphi

        for j in range(0, npsi):
            psi = j * dpsi

            _ntheta = int(math.sin(psi)*ntheta) - 1
            try:
                dtheta = math.pi/_ntheta
            except ZeroDivisionError:
                pass

            for k in range(1, _ntheta):
                theta = k * dtheta
                rot[num_rot][0] = phi
                rot[num_rot][1] = psi
                rot[num_rot][2] = theta
                num_rot = num_rot + 1

    if verbose:
        for phi, psi, theta in rot[0:num_rot]:
            print(phi, psi, theta, sep='  ')

    return rot[0:num_rot]


def main():
    nphi = 1  # 1
    npsi = 20  # 20 -- Currently, you need to double this number for input in the paramfile
    ntheta = 20 # 20

    rots = calculate_rotations(nphi, npsi, ntheta, verbose=True)
    print("Number of rotations:", len(rots))



if __name__ == '__main__':
    main()
