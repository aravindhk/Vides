#!/usr/bin/python

def CNT(N,M):    
    Energy = 0.0
    t = 2.7
    # The following two lines has to be restored for cell-hamiltonians
#    M=(M*4);
#    N=(int)(N/4);
    hamiltonian = []
    for j in range(M):
	for i in range(N):
	    n = i + 1 + j*N
	    hamiltonian.append([n, n, Energy])
    
    for j in range(1, M-1, +4):
	for i in range(1, N):
 	    n = i + 1 + j*N
 	    m = i + (j+1)*N 
	    hamiltonian.append([n, m, t])
            #hamiltonian.append([m, n, t])

    for j in range(3, M-1, +4):
	for i in range(0, N-1):
 	    n = i + 1 + j*N
 	    m = i + 2 + (j+1)*N 
	    hamiltonian.append([n, m, t])
            #hamiltonian.append([m, n, t])

    for j in range(M-1):
	for i in range(N):
 	    n = i + 1 + j*N
 	    m = i + 1 + (j+1)*N 
	    hamiltonian.append([n, m, t])
            #hamiltonian.append([m, n, t])

    for j in range(1, M-1, 4):
        n = j*N + 1
        m = (j+1)*N + N  
	hamiltonian.append([n, m, t])
        #hamiltonian.append([m, n, t])

    for j in range(3, M-1, 4):
        n = (j+1)*N
        m = (j+1)*N + 1
	hamiltonian.append([n, m, t])
        #hamiltonian.append([m, n, t])
	
    return array(hamiltonian)

def GNR(N,M):
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
    coupling2 = [t*1.12, 0, 0, p]
    orbitals = orbitals + zeroes
    hamiltonian.append(orbitals)

    for j in range(M):
	for i in range(N):
	    n = i + 1 + j*N
            p = [n,n]
            p =  p + ene
	    hamiltonian.append(p)
    
    for j in range(1, M-1, +4):
	for i in range(1, N):
 	    n = i + 1 + j*N
 	    m = i + (j+1)*N 
            p = [n,m]
            p =  p + coupling1
	    hamiltonian.append(p)
      #      hamiltonian.append([m, n, t, p, d])

    for j in range(3, M-1, +4):
	for i in range(0, N-1):
 	    n = i + 1 + j*N
 	    m = i + 2 + (j+1)*N 
            p = [n,m]
            p =  p + coupling1
	    hamiltonian.append(p)
       #     hamiltonian.append([m, n, t, p, d])


# nell'if ripristinare il fattore t*1.12  
    for j in range(0, M-1, +4):
	for i in range(N):
 	    n = i + 1 + j*N
 	    m = i + 1 + (j+1)*N
	    if i == 0:
                  p = [n,m]
                  p =  p + coupling2
	          hamiltonian.append(p)
        #          hamiltonian.append([m, n, t*1.12, p, d])
            else :
                  p = [n,m]
                  p =  p + coupling1
        	  hamiltonian.append(p)
         #         hamiltonian.append([m, n, t, p, d])

    for j in range(1, M-1, +4):
	for i in range(N):
 	    n = i + 1 + j*N
 	    m = i + 1 + (j+1)*N 
            p = [n,m]
            p =  p + coupling1
	    hamiltonian.append(p)
#            hamiltonian.append([m, n, t, p, d])

# nell'if ripristinare il fattore t*1.12 
    for j in range(2, M-1, +4):
	for i in range(N):
 	    n = i + 1 + j*N
 	    m = i + 1 + (j+1)*N
	    if i == (N-1):
                  p = [n,m]
                  p =  p + coupling2
	          hamiltonian.append(p)
 #                 hamiltonian.append([m, n, t*1.12, p, d])
            else :
                  p = [n,m]
                  p =  p + coupling1
        	  hamiltonian.append(p)
  #                hamiltonian.append([m, n, t, p, d])

    for j in range(3, M-1, +4):
	for i in range(N):
 	    n = i + 1 + j*N
 	    m = i + 1 + (j+1)*N 
            p = [n,m]
            p =  p + coupling1
	    hamiltonian.append(p)
   #         hamiltonian.append([m, n, t, p, d])

    print hamiltonian

    return array(hamiltonian,dtype=complex)


from NanoTCAD_ViDES import *
from pylab import *
# atoms per slice (or cell)
atoms = 2
# slices (or cells)
slices = 8

h = GNR(atoms, slices)
h[0][0]=1.0
H = Hamiltonian(atoms, slices)
H.Eupper = 3
H.Elower = 0
H.H = h[:,0:3]
H.mu1=1.8
H.mu2=1.8
H.dE=1e-2;
H.ctest=1+1j;
H.Ei=-H.Phi
H.charge_T()


GNR2=nanoribbon(2,1)
GNR2.Eupper=3;
GNR2.Elower=0;
GNR2.mu1=1.8
GNR2.mu2=1.8;
GNR2.dE=1e-2;
GNR2.charge_T();

#plot(H.E,H.T)
#hold
#plot(GNR2.E,GNR2.T,'o')
#show()

plot(H.charge);
hold
plot(GNR2.charge,'o');
show()
