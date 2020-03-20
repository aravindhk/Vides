#!/usr/bin/python
from NanoTCAD_ViDES import *

############ TEST 1 ################

# I create the Hamiltonian for the SNW
a=array([5.431,0,1,3.85,15])
[x_SNW,y_SNW,z_SNW]=atoms_coordinates_nanowire(a);
save_format_xyz("SNW.xyz",x_SNW/10.0,y_SNW/10.0,z_SNW/10.0,"Si");
[HH,n,Nc]=create_H_from_xyz(x_SNW,y_SNW,z_SNW,10,onsite_Si,thop_Si,3,4);
SNW=Hamiltonian(n,Nc);

SNW.H=HH;

SNW.Elower=0;
SNW.Eupper=10;
SNW.dE=0.1;
SNW.eta=1e-7;
SNW.charge_T();


material="Si"
string="transmission-test-%s.dat" %material;
b=loadtxt(string);

if (max(abs(SNW.T[:size(b[:,1])]-b[:,1]))<5e-5):
    string1="%s PASSED \n" %material;
else:
    string1="%s NOT PASSED \n" %material;

plot(SNW.E,SNW.T);
hold
plot(b[:,0],b[:,1],'o');
#print max(abs(SNW.T[:size(b[:,1])]-b[:,1]))
show()




############ TEST 2 ################

# I create the Hamiltonian for the SNW
a=array([5.431,0,0,7.85,20])
[x_SNW,y_SNW,z_SNW]=atoms_coordinates_nanowire(a);
save_format_xyz("SNW.xyz",x_SNW/10.0,y_SNW/10.0,z_SNW/10.0,"Si");
[HH,n,Nc]=create_H_from_xyz(x_SNW,y_SNW,z_SNW,10,onsite_Si,thop_Si,3,4);
SNW=Hamiltonian(n,Nc);

SNW.H=HH;

SNW.Elower=0;
SNW.Eupper=10;
SNW.dE=0.1;
SNW.eta=1e-7;
SNW.charge_T();


material="Si"
string="transmission-test-%s-2.dat" %material;
b=loadtxt(string);
if (max(abs(SNW.T[:size(b[:,1])]-b[:,1]))<1e-1):
    string2="%s PASSED \n" %material;
else:
    string2="%s NOT PASSED \n" %material;


plot(SNW.E,SNW.T);
hold
plot(b[:,0],b[:,1],'o');
#print max(abs(SNW.T[:size(b[:,1])]-b[:,1]))
#
show()


############ TEST 3 ################

# I create the Hamiltonian for the SNW
a=array([5.431,0,1,7.85,20])
[x_SNW,y_SNW,z_SNW]=atoms_coordinates_nanowire(a);
save_format_xyz("SNW.xyz",x_SNW/10.0,y_SNW/10.0,z_SNW/10.0,"Si");
[HH,n,Nc]=create_H_from_xyz(x_SNW,y_SNW,z_SNW,10,onsite_Si,thop_Si,3,4);
SNW=Hamiltonian(n,Nc);

SNW.H=HH;

SNW.Elower=0;
SNW.Eupper=10;
SNW.dE=0.1;
SNW.eta=1e-7;
SNW.charge_T();


material="Si"
string="transmission-test-%s-1.dat" %material;
b=loadtxt(string);
if (max(abs(SNW.T[:size(b[:,1])]-b[:,1]))<5e-5):
    string3="%s PASSED \n" %material;
else:
    string3="%s NOT PASSED \n" %material;

plot(SNW.E,SNW.T);
hold
plot(b[:,0],b[:,1],'o');
#print max(abs(SNW.T[:size(b[:,1])]-b[:,1]))
show()

print(string1,string2,string3)
fp=open("nanowire_test","w");
fp.write(string1);
fp.write(string2);
fp.write(string3);
fp.close()
