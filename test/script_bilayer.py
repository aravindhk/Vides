#!/usr/bin/python
from NanoTCAD_ViDES import *
from pylab import *
# atoms per slice (or cell)
atoms=2;
slices=8;
thop=-2.7;
tp=-0.35
#4*pi/(3*sqrt(3)*0.144);
delta=0.144*sqrt(3);
k=4
V=0.6

h=zeros((34,3),dtype=complex);
h[0][0]=1;
for i in range(1,17):
    h[i][0]=i
    h[i][1]=i
#    if ((i%2)==0):
#        h[i][2]=V*0.5;
#    else:
#        h[i][2]=-V*0.5;

h[17][0]=1;
h[17][1]=3;
h[17][2]=thop;

h[18][0]=3;
h[18][1]=5;
h[18][2]=thop+thop*exp(k*delta*1j);

h[19][0]=3;
h[19][1]=6;
h[19][2]=tp;

h[20][0]=5;
h[20][1]=7;
h[20][2]=thop;

h[21][0]=7;
h[21][1]=9;
h[21][2]=thop+thop*exp(-k*delta*1j);

h[22][0]=7;
h[22][1]=10;
h[22][2]=tp;

h[23][0]=9;
h[23][1]=11;
h[23][2]=thop;

h[24][0]=11;
h[24][1]=13;
h[24][2]=thop+thop*exp(k*delta*1j);

h[25][0]=11;
h[25][1]=14;
h[25][2]=tp;

h[26][0]=13;
h[26][1]=15;
h[26][2]=thop;

h[27][0]=2;
h[27][1]=4;
h[27][2]=thop;

h[28][0]=4;
h[28][1]=6;
h[28][2]=thop+thop*exp(-k*delta*1j);

h[29][0]=6;
h[29][1]=8;
h[29][2]=thop;

h[30][0]=8;
h[30][1]=10;
h[30][2]=thop+thop*exp(k*delta*1j);

h[31][0]=10;
h[31][1]=12;
h[31][2]=thop;

h[32][0]=12;
h[32][1]=14;
h[32][2]=thop+thop*exp(-k*delta*1j);

h[33][0]=14;
h[33][1]=16;
h[33][2]=thop;

H = Hamiltonian(atoms, slices)
H.Elower = 0
H.Eupper = 9
H.H = h
H.mu1=H.mu2=0.3
H.dE=1e-2;
H.eta=1e-5
for i in range(0,16):
    if ((i%2)==0):
        H.Phi[i]=V*0.5;
    else:
        H.Phi[i]=-V*0.5;
H.charge_T()

savetxt("H",h)

plot(H.E,H.T,'-o')
show()
