from NanoTCAD_ViDES import *

GNR = nanoribbon(3, 1)
GNR.Elower = -3
GNR.Eupper = 3
GNR.dE = 0.05
GNR.Phi = -0.2 * GNR.z
GNR.mu2 = -0.3;
GNR.charge_T()
plot(GNR.E, GNR.T);
show()
print GNR.current();
