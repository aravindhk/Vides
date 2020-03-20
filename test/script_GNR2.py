from NanoTCAD_ViDES import *
GNR=nanoribbon_fast_ohmic(6,1.2);
GNR.mu1=1.8
GNR.mu2=1.8
GNR.dE=1e-2
GNR.Eupper=2;
GNR.Elower=1;
GNR.charge_T()


GNR2=nanoribbon(6,1)
GNR2.mu1=1.8
GNR2.mu2=1.8
GNR2.dE=1e-2
GNR2.Eupper=2;
GNR2.Elower=1;
GNR2.charge_T()

plot(GNR.E,GNR.T)
hold

ind=nonzero(abs(GNR2.T)>10);
GNR2.T[ind]=0
plot(GNR2.E,GNR2.T)
show()


plot(GNR.charge)
hold
plot(GNR2.charge)
show()
