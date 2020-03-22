from NanoTCAD_ViDES import *

xmin = -1.5
xmax = 1.5
ymin = -1.5
ymax = 1.5
zmin = 0
zmax = 15

xg = nonuniformgrid(array([-1.5, 0.2, 0, 0.05, 1.5, 0.2]))
yg = nonuniformgrid(array([-1.5, 0.2, 0, 0.05, 1.5, 0.2]))
zg = linspace(zmin, zmax, 90);

CNT = nanotube(13, 15);

grid = grid3D(xg, yg, CNT.z, CNT.x, CNT.y, CNT.z);
