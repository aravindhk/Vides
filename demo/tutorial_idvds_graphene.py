from NanoTCAD_ViDES import *
import sys

rank = MPI.COMM_WORLD.Get_rank()

# I create the grid
xg = nonuniformgrid(array([-10, 1, 0, 0.05, 10, 1]))

FLAKE = graphene(50);

acc = 0.144;
kF = 2 * pi / (3 * sqrt(3) * acc);
kymax = kF + 2;
Nky = 64;
dk = (kymax - kF) / (Nky * 0.5);
FLAKE.kmax = kF + dk * Nky * 0.5;
FLAKE.kmin = kF - dk * Nky * 0.5;
FLAKE.dk = dk;

grid = grid2D(xg, FLAKE.y, FLAKE.x, FLAKE.y);

# I take care of the solid
SiO2 = region("hex", grid.xmin, grid.xmax, grid.ymin, grid.ymax)
SiO2.eps = 3.9;

top_gate = gate("hex", grid.xmax, grid.xmax, 15, 35);
bottom_gate = gate("hex", grid.xmin, grid.xmin, 15, 35);

top_gate.Ef = -0.1;
bottom_gate.Ef = -0.1;

p = interface2D(grid, SiO2, top_gate, bottom_gate);

p.MPI_kt = "yes"

fraction = 5e-3
dope_reservoir(grid, p, FLAKE, fraction, array([-1, 1, grid.ymin, 15]));
dope_reservoir(grid, p, FLAKE, fraction, array([-1, 1, 35, grid.ymax]));

solve_init(grid, p, FLAKE);

Vdsmin = 0.05;
Vdsmax = 0.5;
Vdstep = 0.05;

Np = int(abs(Vdsmin - Vdsmax) / Vdstep) + 1;
vg = zeros(Np);
current = zeros(Np);
p.underel = 0.1;

counter = 0;
Vds = Vdsmin;
while (Vds <= Vdsmax):

    FLAKE.mu2 = -Vds;
    p.normpoisson = 1e-1;
    p.normd = 5e-3;
    solve_self_consistent(grid, p, FLAKE);
    vg[counter] = Vds;
    current[counter] = FLAKE.current();
    # I save the output files
    if (rank == 0):
        string = "./datiout_idvds/Phi%s.out" % Vds;
        savetxt(string, p.Phi);
        string = "./datiout_idvds/ncar%s.out" % Vds;
        savetxt(string, p.free_charge);
        a = [FLAKE.E, FLAKE.T];
        string = "./datiout_idvds/T%s.out" % Vds;
        savetxt(string, transpose(a));
        string = "./datiout_idvds/jayn%s.out" % Vds;
        fp = open(string, "w");
        string2 = "%s" % current[counter];
        fp.write(string2);
        fp.close();
    counter = counter + 1;
    Vds = Vds + Vdstep;

tempo = [vg, current]
savetxt("./datiout_idvds/idvds.out", transpose(tempo));
MPI.Finalize()
