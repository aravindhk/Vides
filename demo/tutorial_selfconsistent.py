from NanoTCAD_ViDES import *

# I create the grid
x1=array([-2,-1.70172,-1.43483,-1.19604,-0.982388,-0.791223,-0.620181,0]);
x2=array([0.620181,0.791223,0.982388,1.19604,1.43483,1.70172,2]);
xg=concatenate((x1,x2),1);
y1=array([-2,-1.70172,-1.43483,-1.19604,-0.982388,-0.791223,-0.620181]);
y2=array([0.620181,0.791223,0.982388,1.19604,1.43483,1.70172,2]);
yg=concatenate((y1,y2),1);
zg=array([-0.072,30.024]);

CNT=nanotube(13,30);

grid=grid3D(xg,yg,zg,CNT.x,CNT.y,CNT.z);

# Now I define the gate regions
top_gate=gate("hex",grid.xmax,grid.xmax,grid.ymin,grid.ymax,grid.zmin,grid.zmax)
bottom_gate=gate("hex",grid.xmin,grid.xmin,grid.ymin,grid.ymax,grid.zmin,grid.zmax)
top_gate.Ef=-0.6;
bottom_gate.Ef=0;
CNT.Nmodes=2;
CNT.mu2=-0.5
# I take care of the solid
SiO2=region("hex",-2,2,-2,2,grid.zmin,grid.zmax);
SiO2.eps=3.9;

p=interface3D(grid,top_gate,bottom_gate,SiO2);

p.normpoisson=1e-3;
solve_Poisson(grid,p)
p.modespace="yes"

p.normpoisson=1e-1;
p.normd=5e-2;
solve_self_consistent(grid,p,CNT);
savetxt("Phi.out",p.Phi);

temp=p.free_charge/grid.dVe/1e-27;
section("x",temp,0.1,grid);
print CNT.current()

