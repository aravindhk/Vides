from NanoTCAD_ViDES import *
import sys
from module_TMD import *

rank = 0

# I create the grid
xg = nonuniformgrid(array([-2.0, 1, 0, 0.05, 2.0, 1]))

FLAKE = TMD(30.0, "n");

acc = FLAKE.acc;
kF = 2 * pi / (3 * sqrt(3) * acc);
kymax = pi / FLAKE.delta;
Nky = 32.0;
dk = kymax / Nky;
FLAKE.kmax = pi / FLAKE.delta;
FLAKE.kmin = 0;
FLAKE.dk = dk;

FLAKE.dE = 0.001
grid = grid2D(xg, FLAKE.y, FLAKE.x, FLAKE.y);
savetxt("gridx.out", grid.gridx)
savetxt("gridy.out", grid.gridy)

# I take care of the solid
Oxide1 = region("hex", grid.xmin, 0, grid.ymin, grid.ymax)
Oxide1.eps = 3.9;

Oxide2 = region("hex", 0, grid.xmax, grid.ymin, grid.ymax)
Oxide2.eps = 3.9;

top_gate = gate("hex", grid.xmax, grid.xmax, 10.0, 20.0);
bottom_gate = gate("hex", grid.xmin, grid.xmin, 10.0, 20.0);

p = interface2D(grid, Oxide1, Oxide2, top_gate, bottom_gate);

fraction_source = 0.01
fraction_drain = 0.01
dope_reservoir(grid, p, FLAKE, fraction_source, array([-1, 1, grid.ymin, 10.0]));
dope_reservoir(grid, p, FLAKE, fraction_drain, array([-1, 1, 20.0, grid.ymax]));

# solve_init(grid,p,FLAKE);

Vgmin = 0.0;
Vgmax = 1.0;
Vgstep = 0.05;

Np = int(abs(Vgmin - Vgmax) / Vgstep) + 1;
vg = zeros(Np);
current = zeros(Np);
p.underel = 0.1;

counter = 0;
Vgs = Vgmin;
FLAKE.mu1 = -0.0
FLAKE.mu2 = -0.1

while (Vgs <= Vgmax):
    bottom_gate.Ef = -Vgs;
    set_gate(p, bottom_gate)
    top_gate.Ef = -Vgs;
    set_gate(p, top_gate)
    p.normpoisson = 1e-1;
    p.normd = 5e-3;
    solve_self_consistent(grid, p, FLAKE);
    vg[counter] = Vgs;
    current[counter] = FLAKE.current();
    # I save the output files
    if (rank == 0):
        string = "./datiout/Phi%s.out" % Vgs;
        savetxt(string, p.Phi);
        string = "./datiout/ncar%s.out" % Vgs;
        savetxt(string, p.free_charge);
        a = [FLAKE.E, FLAKE.T];
        string = "./datiout/T%s.out" % Vgs;
        savetxt(string, transpose(a));
        string = "./datiout/jayn%s.out" % Vgs;
        fp = open(string, "w");
        string2 = "%s" % current[counter];
        fp.write(string2);
        fp.close();
    counter = counter + 1;
    Vgs = Vgs + Vgstep;

tempo = [vg, current]
savetxt("./datiout/idvgs.out", transpose(tempo));

plot(vg[:counter], current[:counter])
show()
