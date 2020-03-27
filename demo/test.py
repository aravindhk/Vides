from NanoTCAD_ViDES import *
from numpy import genfromtxt

fi = genfromtxt("./datiout_idvds/idvds.out", delimiter = ' ')
plot(fi[:,0],fi[:,1])
show()