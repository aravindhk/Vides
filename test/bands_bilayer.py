#!/usr/bin/python
from NanoTCAD_ViDES import *
from pylab import *

# atoms per slice (or cell)
atoms=2;
slices=8;
thop=-2.7;
tp=-0.35
V=0.6;

#4*pi/(3*sqrt(3)*0.144);
delta=0.144*sqrt(3);

h=zeros((34,3),dtype=complex);
h[0][0]=1;
for i in range(1,17):
    h[i][0]=i
    h[i][1]=i
    if ((i%2)==0):
        h[i][2]=V*0.5;
    else:
        h[i][2]=-V*0.5;


#kmin=-pi/delta;
kmin=0;
k=kmin;
Np=30;
kmax=pi/delta
#kmax=10
dk=(kmax-kmin)/Np;
Emax=10
Emin=0;
kvect=linspace(kmin,kmax,Np)
E=arange(Emin,Emax,0.9999e-2);
X,Y=meshgrid(E,kvect);
Ne=size(E)
Z=zeros((Ne,Np));
i=0;
#k=4;

h[17][0]=1;
h[17][1]=3;
h[17][2]=thop;

h[18][0]=3;
h[18][1]=5;


h[19][0]=3;
h[19][1]=6;
h[19][2]=tp;

h[20][0]=5;
h[20][1]=7;
h[20][2]=thop;

h[21][0]=7;
h[21][1]=9;

    
h[22][0]=7;
h[22][1]=10;
h[22][2]=tp;

h[23][0]=9;
h[23][1]=11;
h[23][2]=thop;

h[24][0]=11;
h[24][1]=13;

    
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

    
h[29][0]=6;
h[29][1]=8;
h[29][2]=thop;

h[30][0]=8;
h[30][1]=10;

h[31][0]=10;
h[31][1]=12;
h[31][2]=thop;

h[32][0]=12;
h[32][1]=14;


h[33][0]=14;
h[33][1]=16;
h[33][2]=thop;

hold=h;
#k=0;
while (k<(kmax*0.99)):
#while (k==0):
    print "*********************"
    print "Size H, type H \n"
    print "*********************\n\n\n\n\n\n\n"

    h[18][2]=thop+thop*exp(k*delta*1j);
    h[21][2]=thop+thop*exp(-k*delta*1j);
    h[24][2]=thop+thop*exp(k*delta*1j);
    h[28][2]=thop+thop*exp(-k*delta*1j);
    h[30][2]=thop+thop*exp(k*delta*1j);
    h[32][2]=thop+thop*exp(-k*delta*1j);

    print h-hold


#H = Hamiltonian(atoms, slices)
#H.Elower = 0
#H.Eupper = 9
#H.H = h
#H.mu1=H.mu2=0.3
#H.dE=1e-2;
#H.eta=1e-10
#H.charge_T()


    H = Hamiltonian(atoms, slices)
    H.Eupper = Emax
    H.Elower = Emin
    H.H = h
    H.mu1=H.mu2=0.3
    H.dE=1e-2;

    hold=h

    H.charge_T()
    print H.E

    Z[:,i]=H.T[0:Ne];
    print i
    i=i+1;
    print k
    k=k+dk

#    plot(H.E,H.T)
#    show()
    del H
    
plt.imshow(Z, interpolation='bilinear', cmap=cm.gray,
           origin='lower', extent=[Emin,Emax,kmin,kmax])
show()


#plot(H.charge);
#hold
#plot(GNR2.charge,'o');
#show()
