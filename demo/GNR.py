from numpy import *


def GNR(N, M):
    Energy = 0.0
    Ene = 0.0
    t = 2.7
    p = 0.0
    d = 0.0
    orbitals = [1, 0]
    hamiltonian = []
    zeroes = [0, 0, 0, 0]
    ene = [Energy, 0, 0, Ene]
    coupling1 = [t, 0, 0, p]
    coupling2 = [t * 1.12, 0, 0, p]
    orbitals = orbitals + zeroes
    hamiltonian.append(orbitals)

    for j in range(M):
        for i in range(N):
            n = i + 1 + j * N
            p = [n, n]
            p = p + ene
            hamiltonian.append(p)

    for j in range(1, M - 1, +4):
        for i in range(1, N):
            n = i + 1 + j * N
            m = i + (j + 1) * N
            p = [n, m]
            p = p + coupling1
            hamiltonian.append(p)
    #      hamiltonian.append([m, n, t, p, d])

    for j in range(3, M - 1, +4):
        for i in range(0, N - 1):
            n = i + 1 + j * N
            m = i + 2 + (j + 1) * N
            p = [n, m]
            p = p + coupling1
            hamiltonian.append(p)
    #     hamiltonian.append([m, n, t, p, d])

    # nell'if ripristinare il fattore t*1.12
    for j in range(0, M - 1, +4):
        for i in range(N):
            n = i + 1 + j * N
            m = i + 1 + (j + 1) * N
            if i == 0:
                p = [n, m]
                p = p + coupling2
                hamiltonian.append(p)
            #          hamiltonian.append([m, n, t*1.12, p, d])
            else:
                p = [n, m]
                p = p + coupling1
                hamiltonian.append(p)
        #         hamiltonian.append([m, n, t, p, d])

    for j in range(1, M - 1, +4):
        for i in range(N):
            n = i + 1 + j * N
            m = i + 1 + (j + 1) * N
            p = [n, m]
            p = p + coupling1
            hamiltonian.append(p)
    #            hamiltonian.append([m, n, t, p, d])

    # nell'if ripristinare il fattore t*1.12
    for j in range(2, M - 1, +4):
        for i in range(N):
            n = i + 1 + j * N
            m = i + 1 + (j + 1) * N
            if i == (N - 1):
                p = [n, m]
                p = p + coupling2
                hamiltonian.append(p)
            #                 hamiltonian.append([m, n, t*1.12, p, d])
            else:
                p = [n, m]
                p = p + coupling1
                hamiltonian.append(p)
    #                hamiltonian.append([m, n, t, p, d])

    for j in range(3, M - 1, +4):
        for i in range(N):
            n = i + 1 + j * N
            m = i + 1 + (j + 1) * N
            p = [n, m]
            p = p + coupling1
            hamiltonian.append(p)
    #         hamiltonian.append([m, n, t, p, d])

    print hamiltonian

    return array(hamiltonian, dtype=complex)
