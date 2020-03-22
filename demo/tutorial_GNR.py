from NanoTCAD_ViDES import *
GNR=nanoribbon(6,3)
GNR.Phi=-0.2*GNR.z
GNR.mu2=-0.3;
#GNR.Eupper=2;
#GNR.Elower=-2;
GNR.charge_T()
plot(GNR.E,GNR.T);
show()
print "Current\n";
print GNR.current();
