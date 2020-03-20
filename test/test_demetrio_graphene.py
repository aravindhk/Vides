#!/usr/bin/python
from NanoTCAD_ViDES import *
from pylab import *
# atoms per slice (or cell)
atoms=1;
#N=46;
N=5;
slices=4*N;
thop=-2.7;

#4*pi/(3*sqrt(3)*0.144);
delta=0.144*sqrt(3);

#kmin=-pi/delta;
kmin=8;
k=kmin;
Np=30;
#kmax=pi/delta
kmax=9.0
#kmax=10
dk=(kmax-kmin)/Np;
GNR=nanoribbon(1,N*2.5*0.144+(N-1)*0.144*0.5+2*0.144)
zzz=GNR.z[0:-2];
del GNR;

Efield=0;
Emax=0.2;
Emin=-0.2;
dE=1e-3;
kvect=linspace(kmin,kmax,Np)
Ne=int(floor((Emax-Emin)/dE));
Ne=NEmax
E=linspace(Emin,Emax,Ne);
X,Y=meshgrid(E,kvect);
Z=zeros((Ne,Np));

h=zeros((2*slices,3),dtype=complex);
h[0][0]=1;
for i in range(1,slices+1):
    h[i][0]=i
    h[i][1]=i

kk=1;
for ii in range(slices+1,2*slices):
    
    if ((ii%2)==1):

        h[ii][0]=kk;
        h[ii][1]=kk+1;
        h[ii][2]=thop;
    kk=kk+1;

i=0;
while (k<(kmax*0.99)):
#while (k==4):

    flaggo=0;
    kk=1;
    for ii in range(slices+1,2*slices):
        if ((ii%2)==0):
            h[ii][0]=kk;
            h[ii][1]=kk+1;
            if ((flaggo%2)==0):
                h[ii][2]=thop+thop*exp(k*delta*1j);
            else:
                h[ii][2]=thop+thop*exp(-k*delta*1j);
            flaggo=flaggo+1;
        kk=kk+1;


    H = Hamiltonian(atoms, slices)

    H.z=zzz;
    H.Eupper = Emax
    H.Elower = Emin
    H.H = h
    H.mu1=H.mu2=0.3
    H.dE=dE;



    H.Phi=Efield/H.z[-1]*H.z;

    H.charge_T()
    
    Z[:,i]=H.T;
    print i
    i=i+1;
    k=k+dk

#    plot(H.E,H.T)
#    show()
    if (k<(kmax*0.99)): 
        del H
    
#plt.imshow(Z, interpolation='bilinear', cmap=cm.gray,
#           origin='lower', extent=[kmin,kmax,Emin,Emax])
#show()


T=sum(Z,1)*2*dk/(2*pi)*1e3;
plot(H.E,T);
G0=loadtxt("g0.dat");
plot(G0[:,0],G0[:,1]);
show()


# atoms per slice (or cell)
atoms=1;
N=47;
#N=5;
slices=4*N;
thop=-2.7;

#4*pi/(3*sqrt(3)*0.144);
delta=0.144*sqrt(3);

#kmin=-pi/delta;
kmin=8;
k=kmin;
Np=30;
#kmax=pi/delta
kmax=9.0
#kmax=10
dk=(kmax-kmin)/Np;
GNR=nanoribbon(1,N*2.5*0.144+(N-1)*0.144*0.5+2*0.144)
zzz=GNR.z[0:-2];
del GNR;

Efield=-0.005*zzz[-1];
Emax=0.2;
Emin=-0.2;
dE=1e-3;
kvect=linspace(kmin,kmax,Np)
Ne=int(floor((Emax-Emin)/dE));
Ne=NEmax
E=linspace(Emin,Emax,Ne);
X,Y=meshgrid(E,kvect);
Z=zeros((Ne,Np));

h=zeros((2*slices,3),dtype=complex);
h[0][0]=1;
for i in range(1,slices+1):
    h[i][0]=i
    h[i][1]=i

kk=1;
for ii in range(slices+1,2*slices):
    
    if ((ii%2)==1):

        h[ii][0]=kk;
        h[ii][1]=kk+1;
        h[ii][2]=thop;
    kk=kk+1;

i=0;
while (k<(kmax*0.99)):
#while (k==4):

    flaggo=0;
    kk=1;
    for ii in range(slices+1,2*slices):
        if ((ii%2)==0):
            h[ii][0]=kk;
            h[ii][1]=kk+1;
            if ((flaggo%2)==0):
                h[ii][2]=thop+thop*exp(k*delta*1j);
            else:
                h[ii][2]=thop+thop*exp(-k*delta*1j);
            flaggo=flaggo+1;
        kk=kk+1;


    H = Hamiltonian(atoms, slices)

    H.z=zzz;
    H.Eupper = Emax
    H.Elower = Emin
    H.H = h
    H.mu1=H.mu2=0.3
    H.dE=dE;



    H.Phi=Efield/H.z[-1]*H.z;

    H.charge_T()

    Z[:,i]=H.T;
    print i
    i=i+1;
    k=k+dk

#    plot(H.E,H.T)
#    show()
    if (k<(kmax*0.99)): 
        del H
    
#plt.imshow(Z, interpolation='bilinear', cmap=cm.gray,
#           origin='lower', extent=[kmin,kmax,Emin,Emax])
#show()

print T
T=sum(Z,1)*2*dk/(2*pi)*1e3;
plot(H.E,T);
G=loadtxt("g.dat");
plot(G[:,0],G[:,1]);
show()
