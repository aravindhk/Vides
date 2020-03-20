from NanoTCAD_ViDES import *

# The width of the nanoribbon is 1.37 nm, and it is 15 nm long
GNR=nanoribbon(6,15);

# I create the grid
xg=nonuniformgrid(array([-2,0.3,0,0.2,2,0.3]))
yg=nonuniformgrid(array([-1,0.3,0,0.2,1.37,0.2,2.37,0.3]));
grid=grid3D(xg,yg,GNR.z,GNR.x,GNR.y,GNR.z);

# I define Schottky contacts
GNR.contact="Schottky"

# Now I define the gate regions
top_gate=gate("hex",grid.xmax,grid.xmax,grid.ymin,grid.ymax,grid.zmin,grid.zmax)
bottom_gate=gate("hex",grid.xmin,grid.xmin,grid.ymin,grid.ymax,grid.zmin,grid.zmax)
top_gate.Ef=-0.3;
bottom_gate.Ef=-0.3;
GNR.mu2=-0.2
# I take care of the solid
SiO2=region("hex",-2,2,-2,2,grid.zmin,grid.zmax);
SiO2.eps=3.9;

p=interface3D(grid,top_gate,bottom_gate,SiO2);

p.normpoisson=1e-3;
solve_Poisson(grid,p)

p.normpoisson=1e-1;
p.normd=5e-2;
solve_self_consistent(grid,p,GNR);

plot(GNR.E,GNR.T)
show()

