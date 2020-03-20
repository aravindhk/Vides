from NanoTCAD_ViDES import *

# I create the grid
xg=nonuniformgrid(array([-1,0.2,0,0.05,1,0.2]))

flake=graphene(30);
Nk=30;
flake.kmin=7;
flake.kmax=10;
flake.dk=double(flake.kmax-flake.kmin)/Nk

grid=grid2D(xg,flake.y,flake.x,flake.y);

# Now I define the gate regions
top_gate=gate("hex",grid.xmax,grid.xmax,10,20);
bottom_gate=gate("hex",grid.xmin,grid.xmin,10,20);
top_gate.Ef=-0.6;
bottom_gate.Ef=0;

# I take care of the solid
SiO2=region("hex",grid.xmin,grid.xmax,grid.ymin,grid.ymax)
SiO2.eps=3.9;

p=interface2D(grid,top_gate,bottom_gate,SiO2);

#p.normpoisson=1e-3;
#solve_Poisson(grid,p)
p.Phi=loadtxt("Phi.MPI");
p.MPI_kt="yes"

dope_reservoir(grid,p,flake,5e-3,array([-1,1,grid.ymin,10]));
dope_reservoir(grid,p,flake,5e-3,array([-1,1,20,grid.ymax]));

p.normpoisson=1e-1;
p.normd=5e-2;
solve_self_consistent(grid,p,flake);

savetxt("Phi.out",p.Phi);
savetxt("ncar.out",p.free_charge);
a=[flake.E,flake.T];
savetxt("T.out",transpose(a));
MPI_Finalize()

#temp=p.free_charge/grid.dVe/1e-27;
#section("x",temp,0.1,grid);
#print CNT.current()
