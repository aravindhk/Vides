from NanoTCAD_ViDES import *
from GNR import *

h = GNR(5, 8)
H = Hamiltonian(5, 8)

H.H = h

H.Elower = -3
H.Eupper = 3
H.dE = 0.05
H.eta = 1e-5

H.charge_T()

file_out = [H.E, H.T]
savetxt("transmission.dat", transpose(file_out))

Nend = int((H.Eupper - H.Elower)/H.dE)
plot(H.E[:Nend], H.T[:Nend])
show()
