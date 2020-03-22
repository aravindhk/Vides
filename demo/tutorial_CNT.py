from NanoTCAD_ViDES import *
CNT=nanotube(13,20);
CNT.Eupper=1;
CNT.Elower=-1;
#CNT.Elower=-1;
CNT.charge_T()
#CNT.Phi=-2*ones(CNT.n*CNT.Nc)
#CNT.Nmodes=4;
#CNT.mode_charge_T()
CNT.charge_T()
plot(CNT.E,CNT.T)
show()
