from NanoTCAD_ViDES import *

#I create the ribbon
#GNR=nanoribbon(5,3);
#save_format_xyz("ribbon5_long.xyz",GNR.x[:-10],GNR.y[:-10],GNR.z[:-10],"C")
#exit(0);

#RIBBON WITH A HOLE IN CORRESPONDENCE OF ATOM 28 (NEW METHOD)
[x,y,z]=get_xyz_from_file("ribbon5_long_edge.xyz");
[H_GNR_hole,n,Nc]=create_H_from_xyz(x,y,z,1,-2.7,1.45,3);
H_hole=Hamiltonian(n,Nc);
H_hole.H=H_GNR_hole;
H_hole.Eupper=1.5;
H_hole.Elower=-1.5;
H_hole.dE=0.1;
H_hole.charge_T()

#RIBBON WITH A HOLE IN CORRESPONDENCE OF ATOM 28 (OLD METHOD)
[x1,y1,z1]=get_xyz_from_file("ribbon5_long.xyz");
[H_GNR_old,n,Nc]=create_H_from_xyz(x1,y1,z1,1,-2.7,1.45,3);
#ind=nonzero((real(H_GNR_old[:,0])==16)&(real(H_GNR_old[:,1])==21))[0];
#H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==56)&(real(H_GNR_old[:,1])==61))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==61)&(real(H_GNR_old[:,1])==66))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==66)&(real(H_GNR_old[:,1])==71))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==51)&(real(H_GNR_old[:,1])==56))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==56)&(real(H_GNR_old[:,1])==62))[0];
H_GNR_old[ind,2]=0;

ind=nonzero((real(H_GNR_old[:,0])==57)&(real(H_GNR_old[:,1])==62))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==62)&(real(H_GNR_old[:,1])==67))[0];
H_GNR_old[ind,2]=0;

ind=nonzero((real(H_GNR_old[:,0])==50)&(real(H_GNR_old[:,1])==55))[0];
H_GNR_old[ind,2]=0;
ind=nonzero((real(H_GNR_old[:,0])==55)&(real(H_GNR_old[:,1])==60))[0];
H_GNR_old[ind,2]=0;

ind=nonzero((real(H_GNR_old[:,0])==60)&(real(H_GNR_old[:,1])==65))[0];
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
