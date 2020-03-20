from NanoTCAD_ViDES import *

# I define the nanotube
CNT=nanotube(13,15);

# I create the grid
x=nonuniformgrid(array([-2,0.3,0,0.2,2,0.3]))
y=x;

grid=grid3D(x,y,CNT.z,CNT.x,CNT.y,CNT.z);

# Now I define the gate regions
top_gate=gate("hex",2,2,-2,2,5,10)
bottom_gate=gate("hex",-2,-2,-2,2,5,10)

# I take care of the solid
SiO2=region("hex",-2,2,-2,2,grid.gridz[0],grid.gridz[grid.nz-1]);
SiO2.eps=3.9;

# I create the interface
p=interface3D(grid,top_gate,bottom_gate,SiO2);

# I dope the reservoirs
dope_reservoir(grid,p,CNT,5e-3,
               array([grid.xmin,grid.xmax,grid.ymin,grid.ymax,0,5]));

dope_reservoir(grid,p,CNT,5e-3,
               array([grid.xmin,grid.xmax,grid.ymin,grid.ymax,10,15]));


# I work in the mode space, using 2 modes
p.modespace="yes"
CNT.Nmodes=2;

# Vds = 0.1 V
CNT.mu2=-0.1;

# I start the Vgs sweep. In particular 0<=Vgs<=0.5 V, with 
# with 0.05V as voltage step
Vgmin=0.0;
Vgmax=0.5;
Vgstep=0.05;

#I create the vectors in which I store the data
vg=zeros(20);
current=zeros(20);

counter=0;
Vgs=Vgmin;
while (Vgs<=Vgmax):
    # I set the Fermi level of the top and bottom gate
    top_gate.Ef=-Vgs;
    set_gate(p,top_gate);
    bottom_gate.Ef=-Vgs;
    set_gate(p,bottom_gate);

    #If the first voltage, then I compute the initial solution
    if (Vgs==Vgmin):
        # I compute the initial solution
        p.normpoisson=1e-3;
        print "Computing the initial solution"
        solve_init(grid,p,CNT);

    p.normpoisson=1e-1;
    p.normd=5e-2;
    solve_self_consistent(grid,p,CNT);
    vg[counter]=Vgs;
    current[counter]=CNT.current();
    counter=counter+1;
    Vgs=Vgs+Vgstep;

tempo=[vg,current]
savetxt("transfer.out",transpose(tempo));

