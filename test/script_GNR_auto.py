#!/usr/bin/python
from numpy import *
from NanoTCAD_ViDES import *

#from pylab import *
# atoms per slice (or cell)
atoms = 2
# slices (or cells)
slices = 8

writeout("ciao");

H=nanoribbon_fast_ohmic(2,2)
H.Eupper = 3
H.Elower = 0
H.mu1=1.8
H.mu2=1.8
H.dE=1e-2;
H.ctest=1+1j;
H.charge_T()


GNR2=nanoribbon(2,1)
GNR2.Eupper=3;
GNR2.eta=1e-5;
GNR2.Elower=0;
GNR2.mu1=1.8
GNR2.mu2=1.8;
GNR2.dE=1e-2;
GNR2.charge_T();

plot(H.E,H.T)
hold
plot(GNR2.E,GNR2.T,'o')
show()

Np=size(GNR2.charge);
if (max(abs(GNR2.charge-H.charge[:Np]))<1e-3):
    string="PASSED \n"
else:
    string="NOT PASSED \n"


print(size(GNR2.E))
print(size(H.E))

if (max(abs(GNR2.E-H.E[:size(GNR2.E)]))<1e-3):
    string2="PASSED \n"
else:
    string2="NOT PASSED \n"


fp=open("GNR_test","w");
string3="TEST on GNR \n Test on T: %s Test on charge: %s" %(string,string2);
fp.write(string3);
fp.close()


plot(H.charge);
hold
plot(GNR2.charge,'o');
show()
