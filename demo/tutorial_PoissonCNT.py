from NanoTCAD_ViDES import *

#I create the (13,0) CNT, 10 nm long
CNT=nanotube(13,10);

#I create the grid
xg=linspace(-1,1,20);
yg=xg;
grid=grid3D(xg,yg,CNT.z,CNT.x,CNT.y,CNT.z);

#I save the computed grid
savetxt("gridx.out",grid.gridx);
savetxt("gridy.out",grid.gridy);
savetxt("gridz.out",grid.gridz);

# Now I define the gate regions
# The device is a double gate and the oxide thickness is 
# equal to 0.5 nm
# The lateral spacing is 0.5 nm, too
top_gate=gate("hex",1,1,-1,1,grid.gridz[0],grid.gridz[grid.nz-1])
bottom_gate=gate("hex",-1,-1,-1,1,grid.gridz[0],grid.gridz[grid.nz-1])
top_gate.Ef=-0.1;
bottom_gate.Ef=0;


# I take care of the region embedding the CNT, which is SiO2
SiO2=region("hex",-1,1,-1,1,0,grid.gridz[grid.nz-1]);
SiO2.eps=3.9;

# I then define the interface
p=interface3D(grid,top_gate,bottom_gate,SiO2);
p.normpoisson=1e-1;
# Actually the charge term is imposed to zero, so I 
# solve the Laplace equation
solve_Poisson(grid,p);

# I pass the computed potential to the CNT instance.
# Note that the I need to grid.swap array to map points
# in the 3D domain, to points belonging to the CNT domain
CNT.Phi=p.Phi[grid.swap];
CNT.charge_T();

section("z",p.Phi,5,grid);

