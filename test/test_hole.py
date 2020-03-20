from NanoTCAD_ViDES import *

#IDEAL RIBBON
#[x,y,z]=get_xyz_from_file("ribbon3.xyz");
#H_GNR=create_H_from_xyz(x,y,z,1,-2.7,1.45);
#H=Hamiltonian(3,12);
#H.H=H_GNR;
#H.Eupper=1.5;
#H.Elower=-1.5;
#H.dE=0.1;
#H.charge_T()

#RIBBON WITH A HOLE IN CORRESPONDENCE OF ATOM 28 (NEW METHOD)
[x,y,z]=get_xyz_from_file("ribbon5_hole.xyz");
[H_GNR_hole,n,Nc]=create_H_from_xyz(x,y,z,1,-2.7,1.45,3);
H_hole=Hamiltonian(n,Nc);
H_hole.H=H_GNR_hole;
H_hole.Eupper=1.5;
H_hole.Elower=-1.5;
H_hole.dE=0.1;
H_hole.charge_T()

#RIBBON WITH A HOLE IN CORRESPONDENCE OF ATOM 28 (OLD METHOD)
[x1,y1,z1]=get_xyz_from_file("ribbon5.xyz");
[H_GNR_old,n,Nc]=create_H_from_xyz(x1,y1,z1,1,-2.7,1.45,3);
ind=nonzero((real(H_GNR_old[:,0])==28)&(real(H_GNR_old[:,1])==32))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==28)&(real(H_GNR_old[:,1])==33))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==23)&(real(H_GNR_old[:,1])==28))[0];
H_GNR_old[ind,2]=0;
H_old=Hamiltonian(n,Nc);
H_old.H=H_GNR_old;
H_old.Eupper=1.5;
H_old.Elower=-1.5;
H_old.dE=0.1;
H_old.charge_T()

print x,y,z

#plot(H.E,H.T)
#hold
plot(H_hole.E,H_hole.T,'o')
plot(H_old.E,H_old.T)
show()
