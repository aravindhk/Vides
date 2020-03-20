from NanoTCAD_ViDES import *
#a=array([5.4,0,0,10,6])
#rank = MPI.COMM_WORLD.Get_rank()
rank=0;
a=array([5.43,0,1,3.85,10.5])
[x,y,z]=atoms_coordinates_nanowire(a);
save_format_xyz("Z2.xyz",x/10.0,y/10.0,z/10.0,"Si");
[HH,n,Nc]=create_H_from_xyz(x,y,z,10,thop_Si,3,4);
ind=argsort(HH[:,1]);
savetxt("H",(real(HH[:,0:10])),fmt='%10.2f');
savetxt("H.sort",(real(HH[ind,0:10])),fmt='%10.2f');
H=Hamiltonian(n,Nc);
H.bands=1;
H.H=HH;
H.Eupper=10;
H.Elower=0;
H.eta=1e-12;
H.dE=0.1
#MPIze(H)
H.charge_T()
#H.charge_T()
a=[H.E,H.T]
if (rank==0): savetxt("T.SNW",transpose(a));
#MPI.Finalize()

string="transmission-test-Si.dat" 
a=loadtxt(string);
    
plot(H.E,H.T);
hold
plot(a[:,0],a[:,1],'o');
show();
