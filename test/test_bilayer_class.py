from NanoTCAD_ViDES import *

###########################
#FIRST TEST
###########################

BG=bilayer_graphene(10);
BG.dE=1e-3;
BG.eta=1e-5;
BG.kmax=12.5958;
BG.kmin=5.39722;
BG.dk=(BG.kmax-BG.kmin)/32.0;

#-0.2
BG.Eupper = 0.633527 
BG.Elower = -0.433527 
BG.Phi=-0.2*ones(size(BG.Phi));

MPIze_kt(BG);

plot(BG.charge*0.144*sqrt(3)*1e-9);
hold
a=loadtxt("carica.-0.2");
plot(a,'o')
show()

###########################
#SECOND TEST
###########################

BG=bilayer_graphene(10);
BG.dE=1e-3;
BG.eta=1e-5;
BG.kmax=12.5958;
BG.kmin=5.39722;
BG.dk=(BG.kmax-BG.kmin)/32.0;

#0.5
BG.Eupper = 0.258527 
BG.Elower = -0.933527 
BG.Phi=0.5*ones(size(BG.Phi));

MPIze_kt(BG);

plot(BG.charge*0.144*sqrt(3)*1e-9);
hold
a=loadtxt("carica.0.5");
plot(a,'o')
show()


###########################
# THIRD TEST
###########################


BG=bilayer_graphene(10);
BG.dE=1e-3;
BG.eta=1e-5;
BG.kmax=12.5958;
BG.kmin=5.39722;
BG.dk=(BG.kmax-BG.kmin)/32.0;

#lin
BG.Eupper = 0.433527 
BG.Elower = -0.530727 
BG.Phi=BG.y*1e-2;

MPIze_kt(BG);

plot(BG.charge*0.144*sqrt(3)*1e-9);
hold
a=loadtxt("carica.lin");
plot(a,'o')
show()


###########################
# FOURTH TEST
###########################


BG=bilayer_graphene(10);
BG.dE=1e-3;
BG.eta=1e-5;
BG.kmax=12.5958;
BG.kmin=5.39722;
BG.dk=(BG.kmax-BG.kmin)/32.0;

#plusminus
BG.Eupper = 0.533527 
BG.Elower = -0.533527 


i=linspace(0,size(BG.Phi)-1,size(BG.Phi))
i_even=nonzero((i%2)==0);
i_odd=nonzero((i%2)==1);
BG.Phi[i_even]=-0.1;
BG.Phi[i_odd]=0.1;

MPIze_kt(BG);

plot(BG.charge*0.144*sqrt(3)*1e-9);
hold
a=loadtxt("carica.plusminus");
plot(a,'o')
show()
