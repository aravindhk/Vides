# =====================================================================================
#  Copyright (c) 2010-2012, G. Fiori, G. Iannaccone, University of Pisa
#  
#  This file is released under the BSD license.
#  See the file "license.txt" for information on usage and
#  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# ===================================================================================== 

from numpy import *
from NanoTCAD_ViDESmod import *
from section import *
import sys
import types

writeout("\n")
writeout("-------------------------------------------------------------------\n")
writeout("                         NanoTCAD ViDES ")
writeout("                      Version 1.5 (rel-1-6)")
writeout("                    Last Modified 19 May 2016")
writeout("                     Copyright (C) 2004-2016      \n")
writeout("-------------------------------------------------------------------\n")
writeout("\n")

NEmax=5e3;
DIGIT_PRECISION=20;
max_number_of_cores_on_a_server=8;

#I check if mpi4py is installed on the machine or not
try: 
    from mpi4py import MPI
    mpi4py_loaded = True 
    sizeMPI = MPI.COMM_WORLD.Get_size()
except ImportError: 
    mpi4py_loaded = False 


#I check if scipy is installed
#This is needed for graphene_TOB
try:
    from scipy.optimize import newton
    newton_loaded = True
except ImportError:
    newton_loaded = False


#I check if pylab is installed on the machine or not
try: 
    if (mpi4py_loaded):
        if (sizeMPI<=max_number_of_cores_on_a_server):
            from pylab import *
            pylab_loaded = True 
    else:
        from pylab import *
        pylab_loaded = True 
#except ImportError: 
except Exception: 
    pylab_loaded = False 
    writeout("pylab not installed on this machine or not set up correctly DISPLAY variable")

#definition of constants
kboltz=1.3807e-23
hbar=1.05459e-34
m0=9.1095e-31
q=1.60219e-19
eps0=8.85e-12
#Slater-Costner parameter for sp3d5s* tight-binding Hamiltonian in Si
thop_Si=array([-1.95933,-4.24135,-1.52230,3.02562,3.15565,-2.28485,-0.80993,4.10364,-1.51801,-1.35554,2.38479,-1.68136,2.58880,-1.81400]);
onsite_Si=array([-2.15168,4.22925,4.22925,4.22925,13.78950,13.78950,13.78950,13.78950,13.78950,19.11650]);

def MPIze(channel):
    if (mpi4py_loaded):
        del channel.E;
        channel.E=zeros(NEmax);
        Eupper_save=channel.Eupper;
        Elower_save=channel.Elower;
        vt=kboltz*channel.Temp/q;
        sizeMPI = MPI.COMM_WORLD.Get_size()
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
            channel.rank=rank;
        # I compute the maximum and the minimum
        # of the energy interval
        if ((channel.Eupper>900)&(channel.Elower<-900)):
            Eupper=max(max(channel.mu1,max(-channel.Phi)),channel.mu2)+0.5*channel.gap()+10*vt;
            Elower=min(min(channel.mu1,min(-channel.Phi)),channel.mu2)-0.5*channel.gap()-10*vt;
        else:
            Eupper=channel.Eupper;
            Elower=channel.Elower;
#        string="Eupper and Elower %s %s " %(Eupper,Elower)
#        if (rank==0): writeout(string)
        E=arange(Elower,Eupper,channel.dE);
        arraydim=size(E)/sizeMPI;

        excess=size(E)-sizeMPI*arraydim
        if (rank<excess):
            channel.Elower=E[rank*(arraydim+1)];
            channel.Eupper=E[(rank+1)*(arraydim+1)-1];
        else:
            channel.Elower=E[(rank-excess)*arraydim+excess*(arraydim+1)];
            if (rank==(sizeMPI-1)):
                channel.Eupper=E[size(E)-1];
            else:
                channel.Eupper=E[(rank-excess+1)*arraydim-1+excess*(arraydim+1)];

#        string="Inizio rank %s %s %s" %(rank,channel.Elower,channel.Eupper)
#        writeout(string)           
        channel.charge_T();
        #writeout("Finito rank "),rank,channel.Elower,channel.Eupper;

        # I send the charge and the transmission coefficient
        if (rank!=0):
            temp=array(channel.charge);
            MPI.COMM_WORLD.Send([temp, MPI.DOUBLE],dest=0,tag=11);
            del temp;
            
            NPE=zeros(1,int);
            NPE[0]=int(ceil((channel.Eupper-channel.Elower)/channel.dE))+1;
            #size(arange(channel.Elower,channel.Eupper,channel.dE));
            #int((channel.Eupper-channel.Elower)/channel.dE);
            #size(nonzero(channel.E));
            temp=array(channel.T[:NPE[0]]);
            temp2=array(channel.E[:NPE[0]]);
#            NPE[0]=size(temp);
            MPI.COMM_WORLD.Send([NPE, MPI.INT],dest=0,tag=10);
            MPI.COMM_WORLD.Send([temp, MPI.DOUBLE],dest=0,tag=12);
            MPI.COMM_WORLD.Send([temp2, MPI.DOUBLE],dest=0,tag=14);
            #writeout("Spedito rank "),rank
            del temp;
            del temp2;
        else:
            channel.charge=array(channel.charge);
            NNEE=int(ceil((channel.Eupper-channel.Elower)/channel.dE))+1;
#size(arange(channel.Elower,channel.Eupper,channel.dE));
#            NNEE=((channel.Eupper-channel.Elower)/channel.dE);
#            size(nonzero(channel.E));
            channel.T=array(channel.T[:NNEE]);
            channel.E=array(channel.E[:NNEE]);
            for i in range(1,sizeMPI):
                temp=empty(size(channel.charge),dtype=double);
                MPI.COMM_WORLD.Recv([temp, MPI.DOUBLE],source=i,tag=11);
                channel.charge=channel.charge+temp;
                del temp;
                NPE=empty(1,int);
                MPI.COMM_WORLD.Recv([NPE, MPI.INT],source=i,tag=10);
                temp=empty(NPE[0],dtype=double);
                MPI.COMM_WORLD.Recv([temp, MPI.DOUBLE],source=i,tag=12);
                temp2=empty(NPE[0],dtype=double);
                MPI.COMM_WORLD.Recv([temp2, MPI.DOUBLE],source=i,tag=14);
                channel.T=concatenate((channel.T,temp));
                channel.E=concatenate((channel.E,temp2));
                del temp;
                del temp2;
                #writeout("Preso rank "),i

        channel.charge = MPI.COMM_WORLD.bcast(channel.charge, root=0)
        channel.T = MPI.COMM_WORLD.bcast(channel.T, root=0)
        channel.E = MPI.COMM_WORLD.bcast(channel.E, root=0)
        channel.Eupper=Eupper_save;
        channel.Elower=Elower_save;
#        MPI.Finalize();
    else:
        writeout("*********************************")
        writeout("MPI not installed on this machine")
        writeout("*********************************")
    return;

def MPIze_kt(channel):
    if (mpi4py_loaded):

        kmin_save=channel.kmin;
        kmax_save=channel.kmax;
        vt=kboltz*channel.Temp/q;
        sizeMPI = MPI.COMM_WORLD.Get_size()
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
            channel.rank=rank;
        # I compute the maximum and the minimum
        # of the wave-vector kt
        kt_max=channel.kmax;
        kt_min=channel.kmin;
        if (rank==0): writeout("kt_max, kt_min"),kt_max,kt_min
        k=arange(kt_min,kt_max,channel.dk);
        arraydim=size(k)/sizeMPI;
        channel.kmin=k[rank*arraydim];
        if (rank==(sizeMPI-1)):
            channel.kmax=k[size(k)-1];
        else:
            channel.kmax=k[(rank+1)*arraydim-1];
        
        channel.charge_T();

        NE=size(channel.E);

        # I send the charge and the transmission coefficient
        if (rank!=0):
            temp=array(channel.charge);
            MPI.COMM_WORLD.Send([temp, MPI.DOUBLE],dest=0,tag=11);
            del temp;

            temp=array(channel.T);
            MPI.COMM_WORLD.Send([temp, MPI.DOUBLE],dest=0,tag=12);
            del temp;
        else:
            channel.charge=array(channel.charge);
            channel.T=array(channel.T);
            for i in range(1,sizeMPI):
                temp=empty(size(channel.charge),dtype=double);
                MPI.COMM_WORLD.Recv([temp, MPI.DOUBLE],source=i,tag=11);
                channel.charge=channel.charge+temp;
                del temp;
                temp=empty(NE,dtype=double);
                MPI.COMM_WORLD.Recv([temp, MPI.DOUBLE],source=i,tag=12);
                channel.T=channel.T+temp;
                del temp;

        channel.charge = MPI.COMM_WORLD.bcast(channel.charge, root=0)
        channel.T = MPI.COMM_WORLD.bcast(channel.T, root=0)
        channel.kmin=kmin_save;
        channel.kmax=kmax_save;
#        MPI.Finalize();

    else:
        writeout("*********************************")
        writeout("MPI not installed on this machine")
        writeout("*********************************")
    return;

def set_gate(interface,gate):
    interface.boundary_conditions[gate.index]=gate.Ef;

def solve_init(grid,interface,channel):
    
    # I get the rank
    if (mpi4py_loaded):
        channel.rank = MPI.COMM_WORLD.Get_rank()
    # I set the rank
    if (mpi4py_loaded):
        rank = MPI.COMM_WORLD.Get_rank()
        interface.rank=rank;
    else:
        interface.rank=0;

    # I first give an estimation of the density of states
    # when computing the flat band potential in the regions
    # where the fixed_charge is not equal to zero, assuming
    # full ionization

    # I save the temperature, mu1, mu2, the potential, n, Nc, Eupper, Elower
#    temp_save=channel.Temp;
    mu1_save=channel.mu1;
    mu2_save=channel.mu2;
    Nc_save=channel.Nc;
    Eupper_save=channel.Eupper;
    Elower_save=channel.Elower;
    boundary_conditions_save=copy(interface.boundary_conditions);
    normpoisson_save=interface.normpoisson;
    
    interface.normpoisson=1e-3;
    # I impose a low-temperature, so to compute the LDOS, instead of the 
    # LDOS multiplied by the Fermi-Dirac
    name=grid.__class__.__name__;
    name_channel=channel.__class__.__name__;
    if (name=="grid3D"):
        if (name_channel=="multilayer_graphene"):
            channel.Nc=8;
            x_save=channel.x
            y_save=channel.y
            z_save=channel.z
            channel.atoms_coordinates();
        else:
            channel.Nc=6;
        channel.Phi=zeros(channel.n*channel.Nc);
        channel.mu1=0;
        channel.mu2=0;
        vt=kboltz*channel.Temp/q;
        channel.Eupper=channel.gap()+10*vt;
        channel.Elower=0;
        # I compute the NEGF
        #    if (interface.modespace=="yes"):
        #        channel.mode_charge_T();
        #    else:
        #        if (interface.MPI=="yes"):
        #            MPIze(channel);
        #        else:
        channel.charge_T();
    
    #    N1D=abs(sum(channel.charge))/(6*channel.Nc)/(3*channel.acc)*1e9;
        Ec=channel.gap()*0.5;
        N1D=sum(abs(channel.charge))/(6*channel.n)/(4*channel.acc)*1e9*exp(Ec/vt);
#    return N1D

        # I compute the mean z: if atoms have a z-coordinate > zmean => I impose the electrochemical potential mu2
        # if atoms have a z-coordinate < zmean => I impose the electrochemical potential mu1
        zmean=(grid.zmin+grid.zmax)*0.5;
        indexS=nonzero((abs(interface.fixed_charge)>1e-20)&(grid.z3D<zmean));
        indexD=nonzero((abs(interface.fixed_charge)>1e-20)&(grid.z3D>=zmean));
        potential=zeros(grid.Np);
        argoS=(abs(interface.fixed_charge[indexS])*grid.surf[indexS,5]/N1D);
        argoD=(abs(interface.fixed_charge[indexD])*grid.surf[indexD,5]/N1D);
        
        potential[indexS]=(vt*(log(exp(argoS)-1))+Ec)*sign(interface.fixed_charge[indexS])+mu1_save;
        potential[indexD]=(vt*(log(exp(argoD)-1))+Ec)*sign(interface.fixed_charge[indexD])+mu2_save;
        
        interface.boundary_conditions[indexS]=potential[indexS];
        interface.boundary_conditions[indexD]=potential[indexD];
        
        
        solve_Poisson(grid,interface);
    elif (name=="grid2D"):
        channel.Nc=8;
        channel.Phi=zeros(channel.n*channel.Nc);
        channel.mu1=0;
        channel.mu2=0;
        vt=kboltz*channel.Temp/q;
        channel.Eupper=channel.gap()+10*vt;
        channel.Elower=0;
        # I compute the NEGF
        #    if (interface.modespace=="yes"):
        #        channel.mode_charge_T();
        #    else:
        #if (interface.MPI_kt=="yes"):
        #    MPIze_kt(channel);
        #else:
        channel.charge_T();
    
        Ec=channel.gap()*0.5;
        N1D=sum(abs(channel.charge))/(8*channel.n)/(8*channel.acc)*1e9*exp(Ec/vt);

        # I compute the mean z: if atoms have a z-coordinate > zmean => I impose the electrochemical potential mu2
        # if atoms have a z-coordinate < zmean => I impose the electrochemical potential mu1
        ymean=(grid.ymin+grid.ymax)*0.5;
        indexS=nonzero((abs(interface.fixed_charge)>1e-20)&(grid.y2D<ymean));
        indexD=nonzero((abs(interface.fixed_charge)>1e-20)&(grid.y2D>=ymean));
        potential=zeros(grid.Np);
        argoS=(abs(interface.fixed_charge[indexS])/N1D);
        argoD=(abs(interface.fixed_charge[indexD])/N1D);
        
        potential[indexS]=(vt*(log(exp(argoS)-1))+Ec)*sign(interface.fixed_charge[indexS])+mu1_save;
        potential[indexD]=(vt*(log(exp(argoD)-1))+Ec)*sign(interface.fixed_charge[indexD])+mu2_save;

#        potential[indexS]=Ec;
#        potential[indexD]=Ec;
        
        interface.boundary_conditions[indexS]=potential[indexS];
        interface.boundary_conditions[indexD]=potential[indexD];

        solve_Poisson(grid,interface);
        
    #going back to the old values
    channel.Nc=Nc_save
    channel.mu2=mu2_save;
    channel.mu1=mu1_save;
    channel.Eupper=Eupper_save;
    channel.Elower=Elower_save;
    interface.boundary_conditions=boundary_conditions_save;
    interface.normpoisson=normpoisson_save;
    if (name_channel=="multilayer_graphene"):
        channel.x=x_save
        channel.y=y_save
        channel.z=z_save
        del x_save,y_save,z_save
    #deleting the save variables
    del mu1_save,mu2_save,Nc_save,Eupper_save,Elower_save,boundary_conditions_save;
    
    return;

    


def solve_self_consistent(grid,interface,channel):
    normad=1e30;
#    Phiold=1.0*interface.Phi;
    interface.Phiold=interface.Phi.copy();
    counter=1;
    if (mpi4py_loaded):
        rank = MPI.COMM_WORLD.Get_rank()
    else:
        rank=0;

    while (normad>interface.normd):        
        # I pass the potential in correspondence of the
        # atoms of the material for which I compute the NEGF
        channel.Phi=interface.Phi[grid.swap]
        # I compute the NEGF

#        channel.Phi=zeros(size(grid.swap));
#        savetxt("Phi.before",interface.Phi[grid.swap]);

        if (interface.modespace=="yes"):
            channel.mode_charge_T();
        else:
            if (interface.MPI=="yes"):
                MPIze(channel);
            elif (interface.MPI_kt=="yes"):
                MPIze_kt(channel);
            else:
                channel.charge_T();

#        savetxt("Phi.temp2",interface.Phi);

#        a=[channel.E,channel.T];
#        savetxt("T.temp",transpose(a));

        if (rank==0):
            writeout("--------------------------------------------")
            string="            CURRENT = %s A/m" %(channel.current());
            writeout(string);
            writeout("--------------------------------------------")
            
        # I pass back the free_charge term to 
        # the 3D domain
        interface.free_charge[grid.swap]=channel.charge
        
        if (rank==0): 
            savetxt("ncar.ini",interface.free_charge);
            savetxt("Phi.ini",interface.Phi);

        # I solve Poisson
        solve_Poisson(grid,interface);
#        normad=sqrt(sum((interface.Phiold-interface.Phi)**2));
#        Phiold=zeros(grid.Np);
        normad=max(abs(interface.Phiold-interface.Phi))
        interface.Phi=interface.Phi+(interface.underel)*(interface.Phiold-interface.Phi)
        del interface.Phiold;
#        del Phiold;
#        Phiold=1.0*interface.Phi;
        interface.Phiold=interface.Phi.copy();
        
        if (rank==0): print()
        string="Iteration # %s; ||Phi-Phiold||2 = %s" %(counter,normad)
        if (rank==0): writeout(string) 
        if (rank==0): print() 
        counter=counter+1;
        if (counter>600):
            return;



def solve_Poisson(grid,interface):
    name=grid.__class__.__name__;
    if (name=="grid3D"):
        solvePoisson(grid,interface);
    elif (name=="grid2D"):
        solvePoisson2D(grid,interface);
    elif (name=="grid1D"):
        solvePoisson1D(grid,interface);
    interface.Phi=array(interface.Phi)
    return;

def nonuniformgrid(argu):
    #This is a wrapper for the nonuniformgridmod function
    #so to convert both the argument and the output to numpy arrays
    #I convert the argument in an array
    argarr=array(argu);
    out=nonuniformgridmod(argarr);
    # I return a pyarray
    outarr=array(out);
    return outarr;

#Fermi-Dirac Function
def Fermi(x):
    return 1/(1+exp(x));

def delete_class(class_obj):
    del_class(class_obj);
    del class_obj;
    return;

# This is the class for the nanotube
class nanotube:
    acc=0.144;
    def __init__(self,n,L):
        self.Nc=int(4*(floor((floor(L/nanotube.acc)-1)/3))+2);
        self.n=n;
        self.Phi=zeros(n*self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.dE=1e-3;
        self.thop=-2.7;
        self.eta=1e-5;
        self.mu1=0;
        self.mu2=0;
        self.Temp=300;
        self.contact="doped";
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(self.n*self.Nc);
        self.Nmodes=n;
        self.x=zeros(n*self.Nc);
        self.y=zeros(n*self.Nc);
        self.z=zeros(n*self.Nc);
        self.L=int(self.Nc/2+((self.Nc-1)-self.Nc*0.5)*0.5)*nanotube.acc;
        self.atoms_coordinates();
        self.rank=0;
    def gap(self):
        return abs(2*self.acc*self.thop*pi/(self.n*sqrt(3)*self.acc));
    def atoms_coordinates(self):
        CNT_atoms_coordinates(self);
        self.x=array(self.x);
        self.y=array(self.y);
        self.z=array(self.z);
        return;
    def charge_T(self):
        CNT_charge_T(self);
        self.E=array(self.E);
        self.T=array(self.T);
        self.charge=array(self.charge);
        return;
    def mode_charge_T(self):
        CNTmode_charge_T(self);
        self.E=array(self.E);
        self.T=array(self.T);
        self.charge=array(self.charge);
        return 
    def current(self):
        vt=kboltz*self.Temp/q;
        E=self.E;
        T=self.T;
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE;
        return sum(arg);


# This is the class for the nanoribbon
class GNRphonon:
    def __init__(self,dimer):
        self.N=1000; # number of points qx (longitudinal direction)
        while (((((self.N)-1)%(dimer/2))!=0) | (((self.N)%2)==0)):
            (self.N)+=1;
        self.dimer=dimer;   # numero dimer lines
        self.rank=0;
        self.phi=0.0; # channel potential (midgap)
        self.numberAC=2; # number of AC modes of different simmetry considered (=2: LA+TA, =1: only LA)
        self.Ecutoff=1.0; # cutoff energy 
        self.delta=2; # integer: it specifies the sampling along the kx direction  
        self.deltak=0;
        self.kyE=zeros(dimer);  # transverse electron wavevector 
        self.qy=zeros(dimer); # transverse phonon wavevector 
        self.kx=zeros(self.N); # longitudinal electron wavevector 
        self.qx=zeros(self.N); # longitudinal phonon wavevector 
        self.qx0=0.0; # fixed value for qx (computation of graphene branches)
        self.qy0=0.0; # fixed value for qy (computation of graphene branches)
        self.kxup=0; # maximum value for kx (computation of rates)
        self.kxdown=0; # minimum value for kx (computation of rates)
        self.dim1=self.N;
        self.dim2=dimer;
        self.dim3=6;
        self.mmin=0;
        self.mmax=dimer-1;
        self.kxmin=0;
        self.kxmax=0;
        self.Phi_r1=39.87*10.0;  # first neighbors
        self.Phi_ti1=17.28*10.0;
        self.Phi_to1=9.89*10.0; 
        self.Phi_r2=7.29*10.0;   # second neighbors
        self.Phi_ti2=-4.61*10.0; 
        self.Phi_to2=-0.82*10.0; 
        self.Phi_r3=-2.64*10.0;  # third neighbors
        self.Phi_ti3=3.31*10.0; 
        self.Phi_to3=0.58*10.0; 
        self.Phi_r4=0.10*10.0;  # fourth neighbors
        self.Phi_ti4=0.79*10.0;  
        self.Phi_to4=-0.52*10.0;
        self.energyE=zeros((self.dim1,(2*self.dim2))) # GNR electron curves
        self.energyP2D=zeros((self.dim1,(self.dim2*self.dim3))) # GNR phonon subbranches
        self.minAC=zeros((self.dim2,3));# minimum of the acoustic subbranches
        self.Egraphene=zeros(self.dim3); # graphene
        self.rateAA=zeros((self.dim1,self.dim2));
        self.rateAE=zeros((self.dim1,self.dim2));
        self.rateOA=zeros((self.dim1,self.dim2));
        self.rateOE=zeros((self.dim1,self.dim2));
        self.Dac=4.5*(1.60219e-19); # deformation potential value (eV)
        self.temp=300; # temperature (K)
        self.thop=2.7; # hopping parameter (eV)
        self.aCC=0.144e-9; # lattice constant (m)
    def electron_GNR(self):
        electron_GNR(self);
        self.kx=array(self.kx);
        self.kyE=array(self.kyE);
        self.energyE=array(self.energyE);
        return;
    def phonon_GNR(self):
        phonon_GNR(self);
        self.qx=array(self.qx);
        self.qy=array(self.qy);
        self.energyP2D=array(self.energyP2D);
        return;
    def phonon_graphene(self):
        phonon_graphene(self);
        self.Egraphene=array(self.Egraphene);
        return;
    def rateACABS(self):
        rateACABS(self);
        self.rateAA=array(self.rateAA);
        return;
    def rateACEM(self):
        rateACEM(self);
        self.rateAE=array(self.rateAE);
        return;
    def rateOPTABS(self):
        rateOPTABS(self);
        self.rateOA=array(self.rateOA);
        return;
    def rateOPTEM(self):
        rateOPTEM(self);
        self.rateOE=array(self.rateOE);
        return;

# This is the class for the nanoribbon
class nanoribbon:
    acc=0.144;
    def __init__(self,n,L):
        self.Nc=int(4*(int((int(L/nanoribbon.acc)-1)/3))+2);
        self.n=n;
        self.Phi=zeros(n*self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.dE=1e-3;
        self.thop=-2.7;
        self.eta=1e-5;
        self.mu1=0;
        self.mu2=0;
        self.Temp=300;
        self.contact="doped";
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(self.n*self.Nc);
        self.defects="no";
        self.roughness="no";
        self.rank=0;
        self.atoms_coordinates();
    def atoms_coordinates(self):
        GNR_atoms_coordinates(self);
        self.x=array(self.x);
        self.y=array(self.y);
        self.z=array(self.z);
        return;
    def gap(self):
        return GNRgap(self);
    def charge_T(self):
        GNR_charge_T(self);
        self.E=array(self.E);
        self.T=array(self.T);
        self.charge=array(self.charge);
        return;
    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);

# This is the class for the graphene
class graphene:
    acc=0.144;
    n=1;
    def __init__(self,L):
        self.Nc=int(4*(floor((floor(L/graphene.acc)-1)/3)));
        self.Phi=zeros(self.Nc);
        self.Ei=zeros(self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.delta=sqrt(3)*graphene.acc;
        self.kmax=pi/self.delta;
        self.kmin=0;
        self.dk=0.1;
        self.dE=1e-3;
        self.thop=-2.7;
        self.eta=1e-8;
        self.mu1=0.0;
        self.mu2=0.0;
        self.Temp=300;
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(self.Nc);
        self.rank=0;
        self.atoms_coordinates();
        self.gap();
        self.T2D="no"
    def atoms_coordinates(self):
        GNR_atoms_coordinates(self);
        self.y=array(self.z);
        self.x=zeros(size(self.y));
        del self.z;
        return;
    def gap(self):
        return 0;
    def charge_T(self):
        # Number of slices and atoms
        slices=self.Nc;
        atoms=1;
        # I define the vector of the k-wave vector
        kvect=arange(self.kmin,self.kmax,self.dk)
        # I start defining the Hamiltonian for the graphene flake
        h=zeros((2*slices,3),dtype=complex);
        h[0][0]=1;
        for i in range(1,slices+1):
            h[i][0]=i
            h[i][1]=i
        
        kk=1;
        for ii in range(slices+1,2*slices):
            if ((ii%2)==1):
                h[ii][0]=kk;
                h[ii][1]=kk+1;
                h[ii][2]=self.thop;
            kk=kk+1;
        
        # I then compute the charge and the T for each energy and k and perform the integral
        i=0;
        k=self.kmin;
        H = Hamiltonian(atoms, slices)
        if (self.T2D=="yes"):
            EE=arange(self.Elower,self.Eupper,self.dE);
            kvect=arange(self.kmin,self.kmax+self.dk,self.dk);
            X,Y=meshgrid(EE,kvect);
            Z=zeros((size(EE),size(kvect)))
        while (k<=(self.kmax+self.dk*0.5)):
            if (self.rank==0): writeout("----------------------------------")
            string="    kx range: [%s,%s] " %(self.kmin,self.kmax);
            if (self.rank==0): writeout(string) 
            string="    iteration %s " %i;
            if (self.rank==0): writeout(string);
            if (self.rank==0): writeout("----------------------------------")
            flaggo=0;
            kk=1;
            # I fill the Hamiltonian for the actual wavevector k in the cycle
            for ii in range(slices+1,2*slices):
                if ((ii%2)==0):
                    h[ii][0]=kk;
                    h[ii][1]=kk+1;
                    if ((flaggo%2)==0):
                        h[ii][2]=self.thop+self.thop*exp(k*self.delta*1j);
                    else:
                        h[ii][2]=self.thop+self.thop*exp(-k*self.delta*1j);
                    flaggo=flaggo+1;
                kk=kk+1;
                
            H.Eupper = self.Eupper;
            H.Elower = self.Elower;
            H.rank=self.rank;
            H.H = h
            H.dE=self.dE;
            H.Phi=self.Phi;
            H.Ei=-self.Phi;
            H.eta=self.eta;
            H.mu1=self.mu1;
            H.mu2=self.mu2;
            H.Egap=self.gap();
            
            
            # I then compute T and the charge for the actual kx
            H.charge_T()
            
            # I sum up all the contribution
            if (i==0):
                self.E=H.E;
                # the factor 2 is because I integrate over kx>0
                self.T=H.T*(2*self.dk/(2*pi));
                self.charge=H.charge*(2*self.dk/(2*pi));
            else:
                # the factor 2 is because I integrate over kx>0
                self.T=self.T+H.T*(2*self.dk/(2*pi));
                self.charge=self.charge+H.charge*(2*self.dk/(2*pi));

            if (self.T2D=="yes"):
                Z[:,i]=H.T[:size(EE)];
            k=k+self.dk
            i=i+1;

        if (self.T2D=="yes"):
            plt.imshow(Z, interpolation='bilinear', cmap=cm.gray,
                       origin='lower', extent=[self.kmin,self.kmax,self.Elower,self.Eupper])
            show()

        del H;
        self.E=array(self.E);
        self.T=array(self.T)*1e9;
        self.charge=array(self.charge)*1e9;
        del kvect,h;
        return;

    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);

# This is the class for the graphene bilayer
class bilayer_graphene:
    acc=0.144;
    acc_p=0.35;
    n=2;
    def __init__(self,L):
        self.Nc=int(4*(floor((floor(L/bilayer_graphene.acc)-1)/3)));
        self.n=2;
        self.Phi=zeros(bilayer_graphene.n*self.Nc);
        self.Ei=zeros(bilayer_graphene.n*self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.delta=sqrt(3)*bilayer_graphene.acc;
        self.kmax=pi/self.delta;
        self.kmin=0;
        self.dk=0.1;
        self.dE=1e-3;
        self.thop=-2.7;
        self.tp=-0.35;
        self.eta=1e-8;
        self.mu1=0.0;
        self.mu2=0.0;
        self.Temp=300;
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(bilayer_graphene.n*self.Nc);
        self.rank=0;
        self.atoms_coordinates();
        self.gap();
        self.T2D="no"
    def atoms_coordinates(self):
        n_save=self.n;
        self.n=1;
        GNR_atoms_coordinates(self);
        ydown=array(self.z);
        yup=ydown-self.acc*0.5;
        NN=size(ydown);
        kkk=0;
        self.y=zeros(2*NN);
        for i in range(0,NN):
            self.y[kkk]=ydown[i];
            self.y[kkk+1]=yup[i];
            kkk=kkk+2;
        self.x=zeros(size(self.y));
        i=linspace(0,size(self.y)-1,size(self.y))
        i_even=nonzero((i%2)==0);
        i_odd=nonzero((i%2)==1);
        self.x[i_even]=0;
        self.x[i_odd]=bilayer_graphene.acc_p;
        del self.z,i,i_even,i_odd;
        self.n=n_save;
        return;
    def gap(self):
        # This is an rough exstimation of 
        # the Energy gap: for sure this is
        # the largest attainable value, within
        # the pz tight-binding model
        return abs(self.tp);
    def charge_T(self):
        # Number of slices and atoms
        slices=self.Nc;
        atoms=self.n;
        # I define the vector of the k-wave vector
        kvect=arange(self.kmin,self.kmax,self.dk)
        # I start defining the Hamiltonian for the bilayer graphene
        h=zeros((4*slices+2*(slices/4)-2,3),dtype=complex);
        h[0][0]=1;
        for i in range(1,2*slices+1):
            h[i][0]=i
            h[i][1]=i
            h[i][2]=0.0;
        
        # I then compute the charge and the T for each energy 
        # and k and perform the integral
        i=0;
        k=self.kmin;
        H = Hamiltonian(atoms, slices)
        if (self.T2D=="yes"):
            EE=arange(self.Elower,self.Eupper,self.dE);
            kvect=arange(self.kmin,self.kmax+self.dk,self.dk);
            X,Y=meshgrid(EE,kvect);
            Z=zeros((size(EE),size(kvect)))
        while (k<=(self.kmax+self.dk*0.5)):
            if (self.rank==0): writeout("----------------------------------")
            string="    kx range: [%s,%s] " %(self.kmin,self.kmax);
            if (self.rank==0): writeout(string);
            string="    k: %s " %k;
            if (self.rank==0): writeout(string);
            if (self.rank==0): writeout("----------------------------------")

            # -------------------------------------------------
            # BEGINNING OF THE HAMILTONIAN DEFINITION
            # FOR THE GRAPHENE BILAYER
            # -------------------------------------------------
            
            # I work on the bottom graphene layer
            kk=1;
            flaggo=0;
            for ii in range(2*slices+1,3*slices):
                if ((ii%2)==1):
                    h[ii][0]=kk;
                    h[ii][1]=kk+2;
                    h[ii][2]=self.thop;
                    kk=kk+2;
                else:
                    h[ii][0]=kk;
                    h[ii][1]=kk+2;
                    if ((flaggo%2)==0):
                        h[ii][2]=self.thop+self.thop*exp(k*self.delta*1j);
                    else:
                        h[ii][2]=self.thop+self.thop*exp(-k*self.delta*1j);
                    kk=kk+2;
                    flaggo=flaggo+1;

            # I work on the top graphene layer
            kk=2;
            flaggo=1;
            for ii in range(3*slices,4*slices-1):
                if ((ii%2)==0):
                    h[ii][0]=kk;
                    h[ii][1]=kk+2;
                    h[ii][2]=self.thop;
                    kk=kk+2;
                else:
                    h[ii][0]=kk;
                    h[ii][1]=kk+2;
                    if ((flaggo%2)==0):
                        h[ii][2]=self.thop+self.thop*exp(k*self.delta*1j);
                    else:
                        h[ii][2]=self.thop+self.thop*exp(-k*self.delta*1j);
                    kk=kk+2;
                    flaggo=flaggo+1;


            # I now work on the perpendicular hopping parameter
            kk=3;
            for ii in range(4*slices-1,4*slices+int(slices/2)-2):
                h[ii][0]=kk;
                h[ii][1]=kk+3;
                h[ii][2]=self.tp;
                kk=kk+4;

            # -------------------------------------------------
            #          END OF THE HAMILTONIAN
            # -------------------------------------------------

            H.Eupper = self.Eupper;
            H.Elower = self.Elower;
            H.H = h
            H.rank=self.rank;
            H.dE=self.dE;
            H.Phi=self.Phi;
            ind_even=arange(0,size(H.Phi),2);
            ind_odd=ind_even+1;
            H.Ei[ind_even]=-(self.Phi[ind_even]+self.Phi[ind_odd])*0.5;
            H.Ei[ind_odd]=-(self.Phi[ind_even]+self.Phi[ind_odd])*0.5;
            H.Ei_flag="no"
            H.eta=self.eta;
            H.mu1=self.mu1;
            H.mu2=self.mu2;
            H.Egap=self.gap();
            
#            return H.H

            # I then compute T and the charge for the actual kx
            H.charge_T()

            # I sum up all the contribution
            if (i==0):
                self.E=H.E;
                # the factor 2 is because I integrate over kx>0
                self.T=H.T*(2*self.dk/(2*pi));
                self.charge=H.charge*(2*self.dk/(2*pi));
#                self.charge=H.charge;
            else:
                # The spin is taken into account in the integral for the current
                # the factor 2 is because I integrate over kx>0
                self.T=self.T+H.T*(2*self.dk/(2*pi));
                # 2 because I take into account
                # that I integrate over kx>0
                self.charge=self.charge+H.charge*(2*self.dk/(2*pi));
            
            if (self.T2D=="yes"):
                Z[:,i]=H.T[:size(EE)];
            k=k+self.dk
            i=i+1;

        if (self.T2D=="yes"):
            plt.imshow(Z, interpolation='bilinear', cmap=cm.gray,
                       origin='lower', extent=[self.kmin,self.kmax,self.Elower,self.Eupper])
            show()

        del H;
        self.E=array(self.E);
        self.T=array(self.T)*1e9;
        self.charge=array(self.charge)*1e9;
#        self.charge=array(self.charge);
        del kvect,h;
        return;

    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);



# This is the class for the general Hamiltonian
class Hamiltonian:
    def __init__(self, n, Nc):
        self.Nc=Nc;
        self.n=n;
        self.x=zeros(n*self.Nc);
        self.y=zeros(n*self.Nc);
        self.z=zeros(n*self.Nc);
        self.Phi=zeros(n*self.Nc);
        self.Ei=zeros(n*self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.dE=0.001;
        self.eta=1e-8;
        self.mu1=0;
        self.mu2=0;
        self.Temp=300;
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(n*self.Nc);
        self.Egap=0;
        self.rank=0;
        # if this flag is set to "yes" then Ei=-Phi
        self.Ei_flag="yes"
# The +1 will be then replaced by the number of orbitals per atoms in the nearest neighbourgh approximation
#	self.H=zeros((((Nc*n)*(Nc*n+1)/2),2+100+10));
    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);
    def charge_T(self):
        if (self.Ei_flag=="yes"):
            self.Ei=-self.Phi;
        H_charge_T(self);
        self.E=array(self.E);
        self.T=array(self.T);
        self.charge=array(self.charge);
    def gap(self):
        return 0.5;
# This is the class for the zincblend structures
# This is the class for the zincblend structures
class Zincblend:
    def __init__(self, material, sqci, tilt, edge, zmax):
        self.material = material
        if self.material == 'Si':
            self.aux = [-2.15168, 
                         4.22925, 
                         19.11650,
                         13.78950, 
                         -1.95933,
                         -4.24135,   
                         -1.52230,   
                         3.02562,     
                         3.15565,    
                         -2.28485,   
                         -0.80993,  
                         4.10364,    
                         -1.51801,   
                         -1.35554,   
                         2.38479,    
                         -1.68136,   
                         2.58880,    
                         -1.81400,
                         ]
            self.skparameters = array(self.aux, dtype=float)
            self.a0 = 5.431
            self.flag = 0
        if self.material == 'Ge':
            self.aux = [-1.95617,
                         5.30970,     
                         19.29600,   
                         13.58060,   
                         -1.39456,    
                         -3.56680,   
                         -2.01830,   
                         2.73135,     
                         2.68638,   
                         -2.64779,   
                         -1.12312,    
                         4.28921,     
                         -1.73707,   
                         -2.00115,   
                         2.10953,   
                         -1.32941,   
                         2.56261,    
                         -1.95120    
                         ]
            self.skparameters = array(self.aux, dtype=float)
            self.a0 = 5.6575
            self.flag = 0
            
        if self.material == 'InAs':
            self.aux = [ -5.500420,
                          4.151070,    
                          -0.581930,  
                          6.971630,   
                          19.710590,  
                          19.941380,   
                          13.031690,  
                          13.307090,  
                          -1.694350, 
                          -4.210450,  
                          -2.426740,   
                          -1.159870,  
                          2.598230,    
                          2.809360,   
                          2.067660,   
                          0.937340,   
                          -2.268370,  
                          -2.293090,  
                          -0.899370,  
                          -0.488990, 
                          4.310640,   
                          -1.288950,   
                          -1.731410,  
                          -1.978420,  
                          2.188860,  
                          2.456020,   
                          -1.584610,  
                          2.717930,  
                          -0.505090   
                          ]
            self.skparameters = array(self.aux, dtype=float)
            self.a0 = 6.0583
            self.flag = 1

        self.sqci=sqci;
        self.tilt=tilt;
        self.edge=edge;
        self.zmax=zmax;
        
        layers = int(4*self.zmax/(self.a0) + 1)
        
        if (rank==0):
            writeout("prima="), layers
                
        if layers%4==1:
            layers-=1
        elif layers%4==2:
            layers-=2
        elif layers%4==3:
            layers+=1
        
        if layers%4!=0:
            writeout("INTERRUPT AT WIRE"), material, parameters[0][i]
            writeout("NUMBER OF SLICES NOT MULTIPLE OF 4")
            quit()

        layers += 8
        self.L = (self.a0/4)*(layers-1)
        self.n_aux = int((4*self.edge/self.a0)*(4*self.edge/self.a0)) + 10;
        #forse se ci si leva il +10 non cambia nulla (provare)
        self.Nc_aux = int((4*self.zmax/self.a0)) + 10;
        self.zmax=self.L

        self.atoms=zeros(1);
        self.slices=zeros(1);
        self.max=zeros(1);
        self.rank=0;
        self.deltae=20.0;
        self.ics = zeros(self.n_aux*self.Nc_aux);
        self.ipsilon = zeros(self.n_aux*self.Nc_aux);
        self.zeta = zeros(self.n_aux*self.Nc_aux);
        self.H_aux=zeros( (self.Nc_aux*self.n_aux)*((self.Nc_aux*self.n_aux+1)/2)*(2+100));
        self.H=zeros((((self.Nc_aux*self.n_aux)*(self.Nc_aux*self.n_aux+1)/2),2+100));
        
        self.Zinc();
        self.n = int(self.atoms[0]);
        self.Nc= int(self.slices[0]);
        self.x = self.ics;
        self.y = self.ipsilon;
        self.z = self.zeta;

        self.Phi=zeros(self.n*self.Nc);
        self.Ei=zeros(self.n*self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.dE=0.001;
        self.eta=1e-8;
        self.mu1=0;
        self.mu2=0;
        self.Temp=300;
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(self.n*self.Nc);
        self.Egap=0;
    def gap(self):
        return 0;
    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);
    def charge_T(self):
        H_charge_T(self);
        self.E=array(self.E);
        self.T=array(self.T);
        self.charge=array(self.charge);
        return;
    def Zinc(self):

        writeout(self.skparameters)
#        quit()

        Zinc(self);
#        self.zeta = array(self.zeta);
#        ics1 = []
#        ipsilon1 = []
#        zeta1 = []
#        i = 0
#        j = 0
#        k = 0
#        temp = self.zeta[0]- self.a0
#        zeta1.append(temp)
#        aux = []
#        for ln in self.zeta:
#            if (self.zeta[i]- self.a0) == temp:
#                #temp = self.zeta[i]- self.a0
#                i = i + 1
#                j = j + 1
#            else:
#                zeta1.append(self.zeta[i]- self.a0)
#                temp = self.zeta[i]- self.a0
#                i = i + 1
#                aux.append(j)
#                j=1;
#
#        print aux
#        print self.zeta

        
#        for i in range (100):
        #print zeta1

#        print 'slices =', int(self.slices[0])
#        print 'atoms =', int(self.atoms[0])

#        zeta2 = []
#        for i in range (int(self.slices[0])):
#            for j in range(int(self.atoms[0])):
#                zeta2.append(zeta1[i])
#
#        print 'ECCOLO'
        #print zeta2

#        self.zeta = zeta2

        #print self.zeta

        H_back = []
       
        i = 0
        j = 0
        bound = int(self.max[0]/102)
        writeout(bound)
        
        for i in range ( bound ):
            row = []
            for j in range(102):
                row.append(self.H_aux[j + 102*i])
            H_back.append(row)
            #print row
            del row

            
            
        #print H_back[40]

        new = array(H_back, dtype=complex)

        self.H = new

#        print self.H[17]

#        quit()
       
        return;

def ciccione(vettore,n,Nc,z,a0):
        ics1 = []
        ipsilon1 = []
        zeta1 = []
        i = 0
        j = 0
        k = 0
        temp = z[0]- a0
        z1=[];
        z1.append(temp)
        aux = []
        for ln in arange(0,n*Nc):
            if (z[i]- a0) == temp:
                #temp = self.zeta[i]- self.a0
                i = i + 1
                j = j + 1
            else:
                z1.append(z[i]- a0)
                temp = z[i]- a0
                i = i + 1
                aux.append(j)
                j=1;

#       TODO: the following sum is equal to the total number of 
#       atoms, really present in the simulated nanowire
#
#       Ntot_atoms=sum(aux[:Nc])
#    
#
        array2 = []
        for i in range(Nc):
            k=0;
            if (aux[i]==n):
                for j in arange(sum(aux[:i]),sum(aux[:i])+n):
                    array2.append(vettore[j])
            else:
                for j in arange(sum(aux[:i]),sum(aux[:i])+aux[i]):
                    array2.append(vettore[j]);
                for j in arange(sum(aux[:i])+aux[i],sum(aux[:i])+n):
                    array2.append(0)        
           
        return array(array2, dtype=float); 

# This is the class for graphene within the top-of-the-barrier model
class graphene_TOB:
    if (newton_loaded=="false"):
        print ("scipy not installed")
        print ("Cannot go on")
        exit(0);
    def __init__(self,C1,C2,Vg1,Vg2,Vds,pot_ini):
        self.Vg1=Vg1;
        self.Vg2=Vg2;
        self.Vds=Vds;
        self.mu1=0;
        self.mu2=-Vds;
        self.Temp=300;
        self.thop=2.7
        self.acc=0.144e-9
        self.Pot_c=pot_ini;
        self.C1=C1;
        self.C2=C2;
    def current(self):
        vF=3*self.acc*self.thop*q/2/hbar;
        E=linspace(-3,3,6000)
        dE=E[1]-E[0];
        T=2*q*abs(E+self.Pot_c)/pi/hbar/vF
        vt=kboltz*self.Temp/q;
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*dE;
        return sum(arg)
    def rho(self):
        vF=3*self.acc*self.thop*q/2/hbar;
        return 0.5*(-2*(kboltz*self.Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,(self.mu1+self.Pot_c)/(kboltz*self.Temp/q))+2*(kboltz*self.Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,-(self.mu1+self.Pot_c)/(kboltz*self.Temp/q)))+0.5*(-2*(kboltz*self.Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,(self.mu2+self.Pot_c)/(kboltz*self.Temp/q))+2*(kboltz*self.Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,-(self.mu2+self.Pot_c)/(kboltz*self.Temp/q)))
    
    def charge_I(self):
        Pot_c=self.Pot_c; Vg1=self.Vg1; Vg2=self.Vg2; C1=self.C1; C2=self.C2; Temp=self.Temp; mu1=self.mu1; mu2=self.mu2;
        self.Pot_c=newton(eq1,self.Pot_c,fprime=None,args=(Vg1,Vg2,C1,C2,Temp,mu1,mu2),tol=1e-15,maxiter=10000)
        Id=self.current();
        return self.Pot_c,Id,self.rho()

def eq1(Pot_c,Vg1,Vg2,C1,C2,Temp,mu1,mu2):
    vF=3*0.144e-9*q*2.7/2/hbar;
    charge_S=-2*(kboltz*Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,(mu1+Pot_c)/(kboltz*Temp/q))+2*(kboltz*Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,-(mu1+Pot_c)/(kboltz*Temp/q));
    charge_D=-2*(kboltz*Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,(mu2+Pot_c)/(kboltz*Temp/q))+2*(kboltz*Temp)**2/pi/hbar**2/vF**2*Fermi_Integrals(1,-(mu2+Pot_c)/(kboltz*Temp/q));
    charge=(charge_S+charge_D)*0.5
    return C1*(Vg1-Pot_c)+C2*(Vg2-Pot_c)+q*charge;


class grid3D:
    def __init__(self,*args):
        # I initialize the rank
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
        else:
            rank=0;
        # args is a tuple and len(args) return 
        # the number of arguments
        # the number of arguments can be either 3 or 6
        # if 3, the first three inputs are the grid along the
        # x,y,z axis
        # if 6, the first three inputs are the grid along the
        # x,y,z axis, while the last three inputs are the x-y-z
        # coordinates of the atoms
        if (len(args)>3):
            xg=around(args[0],5);
            yg=around(args[1],5);
            zg=around(args[2],5);
            xC=around(args[3],5);
            yC=around(args[4],5);
            zC=around(args[5],5);
            npC=size(xC);
        else:
            xg=around(args[0],5);
            yg=around(args[1],5);
            zg=around(args[2],5);
            npC=0;

        #I create the grid
        if (npC!=0):
            #find the unique values for xC,yC and zC
            uxC=unique(xC);
            uyC=unique(yC);
            uzC=unique(zC);
            
            # I find the only the additional values which are in xg and not in uxC
            # the same for the other axis
            exg=intersect1d(setxor1d(xg,uxC),xg);
            eyg=intersect1d(setxor1d(yg,uyC),yg);
            ezg=intersect1d(setxor1d(zg,uzC),zg);

        if (npC!=0):
            x=unique(concatenate((uxC,xg),1));
            y=unique(concatenate((uyC,yg),1));
            z=unique(concatenate((uzC,zg),1));
        else:
            x=xg;
            y=yg;
            z=zg;


        # I start to compute the volume associated to each grid point
        X,Y=meshgrid(x,y);

        #Number of points
        nx=size(x);
        ny=size(y);
        nxy=nx*ny;
        nz=size(z);
        Np=nxy*nz;
        string="Number of grid points %s " %Np
        if (rank == 0): writeout(string)

        ####################################################################################
        #I create the Volume elements using the sorted grid
        xd=avervect(x);
        yd=avervect(y);
        zd=avervect(z);
        X,Y=meshgrid(x,y);
        X,Z=meshgrid(x,z);

        XD,ZD=meshgrid(xd,zd);
        surfxz=XD*ZD;
        YD,ZD=meshgrid(yd,zd);
        surfyz=YD*ZD;
        XD,YD=meshgrid(xd,yd);
        surfxy=XD*YD;

        #The volumes for the sorted grid are finally computed
        a,b=meshgrid((XD*YD).flatten(),zd);
        dVes=abs((a*b).flatten());

        if (rank == 0): writeout("Volumes computed")
        
        ####################################################################################
        # I create the dist vectors
        dists=zeros((Np,6));

        # I take care of dists[:,1]
        i=arange(0,nx);
        ip1=i+1;
        ip1[nx-1]=nx-1;
        xdistp=x[ip1]-x[i];
        dists[:,1]=meshgrid(meshgrid(xdistp,y)[0].flatten(),z)[0].flatten();
        del ip1,xdistp;

        # I take care of dists[:,0]
        im1=i-1;
        im1[0]=0;
        xdistm=x[i]-x[im1];
        dists[:,0]=meshgrid(meshgrid(xdistm,y)[0].flatten(),z)[0].flatten();
        del i,im1,xdistm;
        
        # I take care of dists[:,3]
        j=arange(0,ny);
        jp1=j+1;
        jp1[ny-1]=ny-1;
        ydistp=y[jp1]-y[j];
        dists[:,3]=meshgrid(meshgrid(x,ydistp)[1].flatten(),z)[0].flatten();
        del jp1,ydistp;

        # I take care of dists[:,2]
        jm1=j-1;
        jm1[0]=0;
        ydistm=y[j]-y[jm1];
        dists[:,2]=meshgrid(meshgrid(x,ydistm)[1].flatten(),z)[0].flatten();
        del j,jm1,ydistm;

        # I take care of dists[:,5]
        k=arange(0,nz);
        kp1=k+1;
        kp1[nz-1]=nz-1;
        zdistp=z[kp1]-z[k];
        dists[:,5]=meshgrid(meshgrid(x,y)[1].flatten(),zdistp)[1].flatten();
        del kp1,zdistp;

        # I take care of dists[:,4]
        km1=k-1;
        km1[0]=0;
        zdistm=z[k]-z[km1];
        dists[:,4]=meshgrid(meshgrid(x,y)[1].flatten(),zdistm)[1].flatten();
        del k,km1,zdistm;
        
        
        ####################################################################################
        #Now I work on the surfaces

        surfs=zeros((Np,6));

        #surf 0
        XD,YD=meshgrid(xd,yd)
        ##YD[:,0]=0;
        a,b=meshgrid(YD.flatten(),zd)
        surfs[:,0]=abs((a*b).flatten());
        #surf 1
        XD,YD=meshgrid(xd,yd)
        ##YD[:,nx-1]=0;
        a,b=meshgrid(YD.flatten(),zd)
        surfs[:,1]=abs((a*b).flatten());
        #surf 2
        XD,YD=meshgrid(xd,yd)
        ##XD[0,:]=0;
        a,b=meshgrid(XD.flatten(),zd)
        surfs[:,2]=abs((a*b).flatten());
        #surf 3
        XD,YD=meshgrid(xd,yd)
        ##XD[ny-1,:]=0;
        a,b=meshgrid(XD.flatten(),zd)
        surfs[:,3]=abs((a*b).flatten());
        #surf 4
        XD,YD=meshgrid(xd,yd)
        a,b=meshgrid((XD*YD).flatten(),z)
        surfs[:,4]=abs(a.flatten());
        ##surfs[0:nx*ny-1,4]=0;
        #surf 5
        XD,YD=meshgrid(xd,yd)
        a,b=meshgrid((XD*YD).flatten(),z)
        surfs[:,5]=abs(a.flatten());
        ##surfs[(nz-1)*(nx*ny):nz*nx*ny,5]=0;
        
        if (rank == 0): writeout("Surfaces created")
        
        
        ####################################################################################
        #Now I have to go back to the unsorted grid.
        #I create the sorted and unsorted coordinates
        #vectors as a function of the index
        
        #sorted positions
        x3Ds=meshgrid(meshgrid(x,y)[0].flatten(),z)[0].flatten();
        y3Ds=meshgrid(meshgrid(x,y)[1].flatten(),z)[0].flatten();
        z3Ds=meshgrid(meshgrid(x,y)[1].flatten(),z)[1].flatten();
        
        #unsorted positions
        
        if (npC!=0):
            xtemp=unique(concatenate((uxC,xg),1));
            ytemp=unique(concatenate((uyC,yg),1));
            ztemp=unique(concatenate((uzC,zg),1));

            if (rank == 0): writeout("I work on the swap array");
            NpC=size(xC);
            swap=array(arange(0,NpC),int);
            for i in range(0,NpC):
                ixC=nonzero(xtemp==xC[i])[0][0];
                iyC=nonzero(ytemp==yC[i])[0][0];
                izC=nonzero(ztemp==zC[i])[0][0];
                ii=ixC+iyC*nx+izC*nx*ny;
                swap[i]=ii;

            
        ####################################################################################
        # I now fill the attributes of the istance of the grid class
        self.x3D=x3Ds;
        self.y3D=y3Ds
        self.z3D=z3Ds
        self.dVe=dVes;
        self.surf=surfs;
        self.dist=dists;
        self.nx=nx;
        self.ny=ny;
        self.nz=nz;
        self.Np=Np;
        self.gridx=x;
        self.gridy=y;
        self.gridz=z;
        if (npC!=0):
            self.swap=swap;
        self.xmin=min(x);
        self.xmax=max(x);
        self.ymin=min(y);
        self.ymax=max(y);
        self.zmin=min(z);
        self.zmax=max(z);
        return;

class grid2D:
    def __init__(self,*args):

        # I initialize the rank
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
        else:
            rank=0;

        # args is a tuple and len(args) return 
        # the number of arguments
        # the number of arguments can be either 2 or 4
        # if 2, the first two inputs are the grid along the
        # x,y axis
        # if 4, the first two inputs are the grid along the
        # x,y axis, while the last two inputs are the x-y
        # coordinates of the atoms
        if (len(args)>2):
            xg=around(args[0],5);
            yg=around(args[1],5);
            xC=around(args[2],5);
            yC=around(args[3],5);
            npC=size(xC);
        else:
            xg=around(args[0],5);
            yg=around(args[1],5);
            npC=0;

        #I create the grid
        if (npC!=0):
            #find the unique values for xC,yC and zC
            uxC=unique(xC);
            uyC=unique(yC);
            
            # I find the only the additional values which are in xg and not in uxC
            # the same for the other axis
            exg=intersect1d(setxor1d(xg,uxC),xg);
            eyg=intersect1d(setxor1d(yg,uyC),yg);

        if (npC!=0):
            x=unique(concatenate((uxC,xg),1));
            y=unique(concatenate((uyC,yg),1));
        else:
            x=xg;
            y=yg;


        #Number of points
        nx=size(x);
        ny=size(y);
        nxy=nx*ny;
        Np=nxy;
        string="Number of grid points %s " %Np
        if (rank == 0): writeout(string)

        ####################################################################################
        #I create the Volume elements using the sorted grid
        xd=avervect(x);
        yd=avervect(y);
        X,Y=meshgrid(x,y);

        XD,YD=meshgrid(xd,yd);
        surfxy=XD*YD;

        if (rank == 0): writeout("Volumes computed")
        
        ####################################################################################
        # I create the dist vectors
        dists=zeros((Np,4));

        # I take care of dists[:,1]
        i=arange(0,nx);
        ip1=i+1;
        ip1[nx-1]=nx-1;
        xdistp=x[ip1]-x[i];
        dists[:,1]=meshgrid(xdistp,y)[0].flatten();
        del ip1,xdistp;

        # I take care of dists[:,0]
        im1=i-1;
        im1[0]=0;
        xdistm=x[i]-x[im1];
        dists[:,0]=meshgrid(xdistm,y)[0].flatten()
        del i,im1,xdistm;
        
        # I take care of dists[:,3]
        j=arange(0,ny);
        jp1=j+1;
        jp1[ny-1]=ny-1;
        ydistp=y[jp1]-y[j];
        dists[:,3]=meshgrid(x,ydistp)[1].flatten()
        del jp1,ydistp;

        # I take care of dists[:,2]
        jm1=j-1;
        jm1[0]=0;
        ydistm=y[j]-y[jm1];
        dists[:,2]=meshgrid(x,ydistm)[1].flatten();
        del j,jm1,ydistm;
        
        ####################################################################################
        #Now I work on the surface

        XD,YD=meshgrid(xd,yd)
        surfs=(XD*YD).flatten();
        
        if (rank == 0): writeout("Surface created")
        
        
        ####################################################################################
        #Now I have to go back to the unsorted grid.
        #I create the sorted and unsorted coordinates
        #vectors as a function of the index
        
        #sorted positions
        x2Ds=meshgrid(x,y)[0].flatten();
        y2Ds=meshgrid(x,y)[1].flatten();
        
        #unsorted positions
        
        if (npC!=0):
            xtemp=unique(concatenate((uxC,xg),1));
            ytemp=unique(concatenate((uyC,yg),1));

            if (rank == 0): writeout("I work on the swap array");
            NpC=size(xC);
            swap=array(arange(0,NpC),int);
            for i in range(0,NpC):
                ixC=nonzero(xtemp==xC[i])[0][0];
                iyC=nonzero(ytemp==yC[i])[0][0];
                ii=ixC+iyC*nx;
                swap[i]=ii;

            
        ####################################################################################
        # I now fill the attributes of the istance of the grid class
        self.x2D=x2Ds;
        self.y2D=y2Ds
        self.surf=surfs;
        self.dist=dists;
        self.nx=nx;
        self.ny=ny;
        self.Np=Np;
        self.gridx=x;
        self.gridy=y;
        if (npC!=0):
            self.swap=swap;
        self.xmin=min(x);
        self.xmax=max(x);
        self.ymin=min(y);
        self.ymax=max(y);
        return;

class grid1D:
    def __init__(self,*args):

        # I initialize the rank
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
        else:
            rank=0;

        # args is a tuple and len(args) return 
        # the number of arguments
        # the number of arguments can be either 1 or 2
        # if 1, the first input is the grid along the
        # x axis
        # if 2, the first input is the grid along the
        # x axis, while the second input is the x
        # coordinates of the atoms
        if (len(args)>1):
            xg=around(args[0],5);
            xC=around(args[1],5); # attenzione: modificato il 28/5/2011
            npC=size(xC);
        else:
            xg=around(args[0],5);
            npC=0;

        #I create the grid
        if (npC!=0):
            #find the unique values for xC
            uxC=unique(xC);
            
            # I find the only the additional values which are in xg and not in uxC
            exg=intersect1d(setxor1d(xg,uxC),xg);

        if (npC!=0):
            x=unique(concatenate((uxC,xg),1));
        else:
            x=xg;


        #Number of points
        nx=size(x);
        Np=nx;
        if (rank == 0): print(("Number of grid points ",Np))

        ####################################################################################
        # I create the dist vectors
        dists=zeros((Np,4));

        # I take care of dists[:,1]
        i=arange(0,nx);
        ip1=i+1;
        ip1[nx-1]=nx-1;
        xdistp=x[ip1]-x[i];
        dists[:,1]=xdistp;
        del ip1,xdistp;

        # I take care of dists[:,0]
        im1=i-1;
        im1[0]=0;
        xdistm=x[i]-x[im1];
        dists[:,0]=xdistm;
        del i,im1,xdistm;
        
        ####################################################################################
        #Now I have to go back to the unsorted grid.
        #I create the sorted and unsorted coordinates
        #vectors as a function of the index
        
        if (npC!=0):
            xtemp=unique(concatenate((uxC,xg),1));

            if (rank == 0): print("I work on the swap array");
            NpC=size(xC);
            swap=array(arange(0,NpC),int);
            for i in range(0,NpC):
                ixC=nonzero(xtemp==xC[i])[0][0];
                ii=ixC;
                swap[i]=ii;

            
        ####################################################################################
        # I now fill the attributes of the istance of the grid class
        self.x=x;
        self.dist=dists;
        self.nx=nx;
        self.Np=Np;
        self.gridx=x;
        if (npC!=0):
            self.swap=swap;
        self.xmin=min(x);
        self.xmax=max(x);
        return;

class region:
    def __init__(self,*args):
        self.name="none";
        self.geometry="hex";
        self.eps=3.9;
        self.rho=0;
        if (args[0]=="hex"):
            if (len(args)>5):
                self.xmin=args[1];
                self.xmax=args[2];
                self.ymin=args[3];
                self.ymax=args[4];
                self.zmin=args[5];
                self.zmax=args[6];
            elif ((len(args)>3)&(len(args)<=5)):
                self.xmin=args[1];
                self.xmax=args[2];
                self.ymin=args[3];
                self.ymax=args[4];
            elif (len(args)<=3):
                self.xmin=args[1];
                self.xmax=args[2];
    def set_material(self,material):
        if (material.lower()=="sio2"):
            self.eps=3.9;
            self.mel=0.5;
            self.met=0.5;
            self.Egap=8.05;
            self.chi=0.95;
            self.mhole=0.42
        if (material.lower()=="si"):
            self.eps=11.8;
            self.mel=0.916;
            self.met=0.19;
            self.Egap=1.124519;
            self.chi=4.05;
            self.mhole=0.549;

class gate:
    def __init__(self,*args):
        self.geometry="hex";
        self.Ef=0;
        self.wf=4.1;
        if (args[0]=="hex"):
            if (len(args)>5):
                self.xmin=args[1];
                self.xmax=args[2];
                self.ymin=args[3];
                self.ymax=args[4];
                self.zmin=args[5];
                self.zmax=args[6];
            elif ((len(args)>3)&(len(args)<=5)):
                self.xmin=args[1];
                self.xmax=args[2];
                self.ymin=args[3];
                self.ymax=args[4];
            elif (len(args)<=3):
                self.xmin=args[1];
                self.xmax=args[2];
        if (args[0]=="cyl"):
            self.xc=args[1];
            self.yc=args[2];
            self.radius=args[3];
            self.geometry="cyl"
        if (args[0]=="trapz"):
            self.xmin=args[1];
            self.xmax=args[2];
            self.y1=args[3];
            self.z1=args[4];
            self.y2=args[5];
            self.z2=args[6];
            self.y3=args[7];
            self.z3=args[8];
            self.y4=args[9];
            self.z4=args[10];
            self.geometry="trapz"

class interface3D:
    def __init__(self,*args):

        # I set the rank
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
            self.rank=rank;
        else:
            self.rank=0;

        # I compute the number of arguments (classes)
        Narg=size(args);
        # I first find the index of the class grid
        igrid=-10;
        for i in range(0,Narg):
            name=args[i].__class__.__name__
            if (name=="grid3D"):
                igrid=i;
        # If no grid class is specified I exit 
        if (igrid==-10):
            writeout("ERROR: grid not passed to structure")
            return;

        # I create the arrays to be used
        self.eps=zeros(args[igrid].Np);

        # I create the vector, where the boundary conditions
        # are specified:
        # if 2000   : inner point
        # if 1001   : Neumann 1
        # if 1002   : Neumann 2
        # if 1003   : Neumann 3
        # if 1004   : Neumann 4
        # if 1005   : Neumann 5
        # if 1006   : Neumann 6
        # if <= 1000: Fermi level of the gate

        # I start defining all the points as inner points
        self.boundary_conditions=2000*ones(args[igrid].Np);

        ###############################################################################################
        # Now I impose the Neumann Boundary conditions on 
        # the surfaces delimiting the 3D domain
        ###############################################################################################

        # I take care of Neumann1
        indexNeu1=nonzero(args[igrid].x3D==min(args[igrid].gridx));
        self.boundary_conditions[indexNeu1]=1001;

        # I take care of Neumann2
        indexNeu2=nonzero(args[igrid].x3D==max(args[igrid].gridx));
        self.boundary_conditions[indexNeu2]=1002;

        # I take care of Neumann3 
        indexNeu3=nonzero(args[igrid].y3D==min(args[igrid].gridy));
        self.boundary_conditions[indexNeu3]=1003;

        # I take care of Neumann4
        indexNeu4=nonzero(args[igrid].y3D==max(args[igrid].gridy));
        self.boundary_conditions[indexNeu4]=1004


        # I take care of Neumann5 and Neumann6
        indexNeu5=nonzero(args[igrid].z3D==min(args[igrid].gridz));
        self.boundary_conditions[indexNeu5]=1005;
        indexNeu6=nonzero(args[igrid].z3D==max(args[igrid].gridz));
        self.boundary_conditions[indexNeu6]=1006;

        ###############################################################################################
        # I check to which class the args belongs to 
        # and I proceed accordingly
        ###############################################################################################

        for i in range(0,Narg):
            name=args[i].__class__.__name__
            # I check if the class is a gate
            if (name=="gate"):
                #I check if the geometry is an hexahedron
                if (args[i].geometry=="hex"):
                    # I find the indexes of the 3D grid which belongs to the gate
                    # with hex geometry
                    index=nonzero((args[i].xmin<=args[igrid].x3D)&(args[i].xmax>=args[igrid].x3D)&
                                  (args[i].ymin<=args[igrid].y3D)&(args[i].ymax>=args[igrid].y3D)&
                                  (args[i].zmin<=args[igrid].z3D)&(args[i].zmax>=args[igrid].z3D));
                    self.boundary_conditions[index]=args[i].Ef;
                    args[i].index=index;
                if (args[i].geometry=="trapz"):
                    # I find the indexes of the 2D grid which belongs to the gate
                    # with trapezoidal geometry
                    if (args[i].y2==args[i].y1):
                        m1=(args[i].z2-args[i].z1)/(args[i].y2-args[i].y1+1e-3)
                    else:
                        m1=(args[i].z2-args[i].z1)/(args[i].y2-args[i].y1)
                    if (args[i].y3==args[i].y2):
                        m2=(args[i].z3-args[i].z2)/(args[i].y3-args[i].y2+1e-3)
                    else:
                        m2=(args[i].z3-args[i].z2)/(args[i].y3-args[i].y2)
                    if (args[i].y4==args[i].y3):
                        m3=(args[i].z4-args[i].z3)/(args[i].y4-args[i].y3+1e-3)
                    else:
                        m3=(args[i].z4-args[i].z3)/(args[i].y4-args[i].y3)
                    if (args[i].y4==args[i].y1):
                        m4=(args[i].z4-args[i].z1)/(args[i].y4-args[i].y1+1e-3)
                    else:
                        m4=(args[i].z4-args[i].z1)/(args[i].y4-args[i].y1)
    
                    index=nonzero((args[igrid].z3D>=(m1*(args[igrid].y3D-args[i].y1)+args[i].z1))&
                                  (args[igrid].z3D>=(m2*(args[igrid].y3D-args[i].y2)+args[i].z2))&
                                  (args[igrid].z3D<=(m3*(args[igrid].y3D-args[i].y3)+args[i].z3))&
                                  (args[igrid].z3D<=(m2*(args[igrid].y3D-args[i].y1)+args[i].z1))&
                                  (args[i].xmin<=args[igrid].x3D)&(args[i].xmax>=args[igrid].x3D));
                    self.boundary_conditions[index]=args[i].Ef;
                    args[i].index=index;

            elif (name=="region"):
                if (args[i].geometry=="hex"):
                    # I find the indexes of the 3D grid which belongs to the gate
                    # with hex geometry
                    index=nonzero((args[i].xmin<=args[igrid].x3D)&(args[i].xmax>=args[igrid].x3D)&
                                  (args[i].ymin<=args[igrid].y3D)&(args[i].ymax>=args[igrid].y3D)&
                                  (args[i].zmin<=args[igrid].z3D)&(args[i].zmax>=args[igrid].z3D));
                    self.eps[index]=args[i].eps;
            elif (name=="grid3D"):
                #dummy line 
                name;
            else:
                writeout("ERROR: Unrecognized input")
                return;

        ###############################################################################################
        # I fill the field of the interface class
        ###############################################################################################

        #self.boundary already filled
        #self.eps already filled
        self.Phiold=zeros(args[igrid].Np)
        self.Phi=zeros(args[igrid].Np);
        self.normpoisson=1e-3;
        self.tolldomn=1e-1;
        self.underel=0;
        self.free_charge=zeros(args[igrid].Np);
        self.fixed_charge=zeros(args[igrid].Np);
        self.normd=5e-2;
        self.modespace="no"
        self.MPI="no"
        self.MPI_kt="no"
        return;

class interface2D:
    def __init__(self,*args):

        # I set the rank
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
            self.rank=rank;
        else:
            self.rank=0;

        # I compute the number of arguments (classes)
        Narg=size(args);
        # I first find the index of the class grid
        igrid=-10;
        for i in range(0,Narg):
            name=args[i].__class__.__name__
            if (name=="grid2D"):
                igrid=i;
        # If no grid class is specified I exit 
        if (igrid==-10):
            writeout("ERROR: grid not passed to structure")
            return;

        # I create the arrays to be used
        self.eps=zeros(args[igrid].Np);

        # I create the vector, where the boundary conditions
        # are specified:
        # if 2000   : inner point
        # if 1001   : Neumann 1
        # if 1002   : Neumann 2
        # if 1003   : Neumann 3
        # if 1004   : Neumann 4
        # if <= 1000: Fermi level of the gate

        # I start defining all the points as inner points
        self.boundary_conditions=2000*ones(args[igrid].Np);

        ###############################################################################################
        # Now I impose the Neumann Boundary conditions on 
        # the surfaces delimiting the 3D domain
        ###############################################################################################

        # I take care of Neumann1
        indexNeu1=nonzero(args[igrid].x2D==min(args[igrid].gridx));
        self.boundary_conditions[indexNeu1]=1001;

        # I take care of Neumann2
        indexNeu2=nonzero(args[igrid].x2D==max(args[igrid].gridx));
        self.boundary_conditions[indexNeu2]=1002;

        # I take care of Neumann3 
        indexNeu3=nonzero(args[igrid].y2D==min(args[igrid].gridy));
        self.boundary_conditions[indexNeu3]=1003;

        # I take care of Neumann4
        indexNeu4=nonzero(args[igrid].y2D==max(args[igrid].gridy));
        self.boundary_conditions[indexNeu4]=1004


        ###############################################################################################
        # I check to which class the args belongs to 
        # and I proceed accordingly
        ###############################################################################################

        for i in range(0,Narg):
            name=args[i].__class__.__name__
            # I check if the class is a gate
            if (name=="gate"):
                #I check if the geometry is an hexahedron
                if (args[i].geometry=="hex"):
                    # I find the indexes of the 2D grid which belongs to the gate
                    # with hex geometry
                    index=nonzero((args[i].xmin<=args[igrid].x2D)&(args[i].xmax>=args[igrid].x2D)&
                                  (args[i].ymin<=args[igrid].y2D)&(args[i].ymax>=args[igrid].y2D));
                    self.boundary_conditions[index]=args[i].Ef;
                    args[i].index=index;
                #I check if the geometry is an cylindrical
                if (args[i].geometry=="cyl"):
                    # I find the indexes of the 2D grid which belongs to the gate
                    # with cyl geometry
                    index=nonzero(((args[i].xc-args[igrid].x2D)**2+(args[i].yc-args[igrid].y2D)**2)<(args[i].radius)**2);
                    self.boundary_conditions[index]=args[i].Ef;
                    args[i].index=index;
            elif (name=="region"):
                if (args[i].geometry=="hex"):
                    # I find the indexes of the 2D grid which belongs to the gate
                    # with hex geometry
                    index=nonzero((args[i].xmin<=args[igrid].x2D)&(args[i].xmax>=args[igrid].x2D)&
                                  (args[i].ymin<=args[igrid].y2D)&(args[i].ymax>=args[igrid].y2D));
                    self.eps[index]=args[i].eps;
            elif (name=="grid2D"):
                #dummy line 
                name;
            else:
                writeout("ERROR: Unrecognized input")
                return;

        ###############################################################################################
        # I fill the field of the interface class
        ###############################################################################################

        #self.boundary already filled
        #self.eps already filled
        self.Phiold=zeros(args[igrid].Np)
        self.Phi=zeros(args[igrid].Np);
        self.normpoisson=1e-3;
        self.tolldomn=1e-1;
        self.underel=0;
        self.free_charge=zeros(args[igrid].Np);
        self.fixed_charge=zeros(args[igrid].Np);
        self.normd=5e-2;
        self.modespace="no"
        self.MPI="no"
        self.MPI_kt="no"
        return;

class interface1D:
    def __init__(self,*args):

        # I set the rank
        if (mpi4py_loaded):
            rank = MPI.COMM_WORLD.Get_rank()
            self.rank=rank;
        else:
            self.rank=0;

        # I compute the number of arguments (classes)
        Narg=size(args);
        # I first find the index of the class grid
        igrid=-10;
        for i in range(0,Narg):
            name=args[i].__class__.__name__
            if (name=="grid1D"):
                igrid=i;
        # If no grid class is specified I exit 
        if (igrid==-10):
            print("ERROR: grid not passed to structure")
            return;

        # I create the arrays to be used
        self.eps=zeros(args[igrid].Np);
        self.mel=zeros(args[igrid].Np);
        self.met=zeros(args[igrid].Np);
        self.chi=zeros(args[igrid].Np);
        self.Egap=zeros(args[igrid].Np);
        self.fixed_charge=zeros(args[igrid].Np);
        self.mhole=zeros(args[igrid].Np);

        # I create the vector, where the boundary conditions
        # are specified:
        # if 2000   : inner point
        # if 1001   : Neumann 1
        # if 1002   : Neumann 2
        # if <= 1000: Fermi level of the gate

        # I start defining all the points as inner points
        self.boundary_conditions=2000*ones(args[igrid].Np);

        ###############################################################################################
        # Now I impose the Neumann Boundary conditions on 
        # the surfaces delimiting the 3D domain
        ###############################################################################################

        # I take care of Neumann1
        indexNeu1=nonzero(args[igrid].x==min(args[igrid].gridx));
        self.boundary_conditions[indexNeu1]=1001;

        # I take care of Neumann2
        indexNeu2=nonzero(args[igrid].x==max(args[igrid].gridx));
        self.boundary_conditions[indexNeu2]=1002;


        ###############################################################################################
        # I check to which class the args belongs to 
        # and I proceed accordingly
        ###############################################################################################

        for i in range(0,Narg):
            name=args[i].__class__.__name__
            # I check if the class is a gate
            if (name=="gate"):
                #I check if the geometry is an hexahedron
                if (args[i].geometry=="hex"):
                    # I find the indexes of the 2D grid which belongs to the gate
                    # with hex geometry
                    index=nonzero((args[i].xmin<=args[igrid].x)&(args[i].xmax>=args[igrid].x));
                    self.boundary_conditions[index]=args[i].Ef;
                    args[i].index=index;
            elif (name=="region"):
                if (args[i].geometry=="hex"):
                    dist=avervect(args[igrid].x)*1e-9;
                    # I find the indexes of the 2D grid which belongs to the gate
                    # with hex geometry
                    index=nonzero((args[i].xmin<=args[igrid].x)&(args[i].xmax>=args[igrid].x));
                    self.eps[index]=args[i].eps;
                    self.mel[index]=args[i].mel;
                    self.met[index]=args[i].met;
                    self.chi[index]=args[i].chi;
                    self.Egap[index]=args[i].Egap;
                    self.fixed_charge[index]=args[i].rho*dist[index];
                    self.mhole[index]=args[i].mhole;
                    
            elif (name=="grid1D"):
                #dummy line 
                name;
            else:
                print("ERROR: Unrecognized input")
                return;

        ###############################################################################################
        # I fill the field of the interface class
        ###############################################################################################

        #self.boundary already filled
        #self.eps already filled
        self.Phiold=zeros(args[igrid].Np)
        self.Phi=zeros(args[igrid].Np);
        self.normpoisson=1e-3;
        self.tolldomn=1e-1;
        self.underel=0;
        self.free_charge=zeros(args[igrid].Np);
        self.normd=5e-2;
        self.modespace="no"
        self.MPI="no"
        self.MPI_kt="no"
        return;



def dope_reservoir(grid,interface,channel,molar_fraction,bbox):
    name=grid.__class__.__name__;
    if (name=="grid3D"):
        xmin=bbox[0];
        xmax=bbox[1];
        ymin=bbox[2];
        ymax=bbox[3];
        zmin=bbox[4];
        zmax=bbox[5];

        index=nonzero((xmin<=grid.x3D[grid.swap])&(xmax>=grid.x3D[grid.swap])&
                      (ymin<=grid.y3D[grid.swap])&(ymax>=grid.y3D[grid.swap])&
                      (zmin<=grid.z3D[grid.swap])&(zmax>=grid.z3D[grid.swap]))
        interface.fixed_charge[grid.swap[index]]=molar_fraction;
    elif (name=="grid2D"):
        xmin=bbox[0];
        xmax=bbox[1];
        ymin=bbox[2];
        ymax=bbox[3];

        index=nonzero((xmin<=grid.x2D[grid.swap])&(xmax>=grid.x2D[grid.swap])&
                      (ymin<=grid.y2D[grid.swap])&(ymax>=grid.y2D[grid.swap]))
        interface.fixed_charge[grid.swap[index]]=molar_fraction/channel.delta*1e9;
    elif (name=="grid1D"):
        xmin=bbox[0];
        xmax=bbox[1];

        index=nonzero((xmin<=grid.x[grid.swap])&(xmax>=grid.x[grid.swap]));
        interface.fixed_charge[grid.swap[index]]=molar_fraction/(channel.deltaz*channel.deltay)*1e18; 
# MODIFICATO IL 6/6/2011: aggiunto il deltay e deltaz

    return index;

class Device:
    def __init__(self):
        self.Nregions=1;
        self.regions=[];
        self.E=zeros(NEmax);
    def test(self):
        return self.E;
        
def test_var_args(farg, *args):
    writeout("formal arg:"), size(args)
    for arg in args:
        writeout("another arg:"), arg

def avervect(x):
    # This function compute the length of
    # the Voronoi segment of a one-dimensional array x
    nx=size(x);
    xd=zeros(nx);
    xini=x[0];
    xd[0]=abs(x[0]-x[1])*0.5;
    for i in range(1,nx-1):
        xd[i]=abs((x[i+1]-x[i-1])*0.5);
    xd[nx-1]=abs(x[nx-1]-x[nx-2])*0.5
    return xd;

def save_format_xyz(outputfile,x,y,z,atom):
    
    if sys.version > '3':
        import subprocess;
    else:
        import subprocess

    out=[x*10,y*10,z*10]
    fp=open(outputfile,"w");
    fp.write(str(size(x)));
    fp.write("\n");
    fp.write("\n");
    for i in range(0,size(x)):
        string="%s %s %s %s" %(atom,out[0][i],out[1][i],out[2][i]);
        fp.write(string);
        fp.write(" ");
        fp.write("\n");
    fp.close()
    return;

"""def convert_pdb(filename,thop):
    fp=open(filename,"r");
    hh=[];
    atoms=0;
    i=0;
    x=[];
    y=[];
    z=[];
    h=[];
    h.append([1,0,0]);
    for line in fp:
        hh.append(line);
        atoms=atoms+(hh[i].split()).count('HETATM');
        if (((hh[i].split()).count('HETATM')==1)|((hh[i].split()).count('ATOM')==1)):
            x.append((hh[i].split())[5]);
            y.append((hh[i].split())[6]);
            z.append((hh[i].split())[7]);
            h.append([int((hh[i].split())[1]),int((hh[i].split())[1]),0]);
        if ((hh[i].split()).count('CONECT')==1):
            a=(hh[i].split());
            NPV=size(a)-1
            for j in range(0,NPV):
                a1=int(a[1]);
                if (a1<int(a[j+1])):
                    h.append([a1,int(a[j+1]),thop])
        if ((hh[i].split()).count('CRYST1')==1):
            a=(hh[i].split());
            if (double(a[1])>=100):
                deltax=0.0;
            else:
                deltax=double(a[1])/10.0;
            if (double(a[2])>=100):
                deltay=0.0;
            else:
                deltay=double(a[2])/10.0;
            if (double(a[3])>=100):
                deltaz=0.0;
            else:
                deltaz=double(a[3])/10.0;            
    

        i=i+1;
    fp.close()
    H=array(h,dtype(complex));
    x=array(x,dtype(float))/10.0;
    y=array(y,dtype(float))/10.0;
    z=array(z,dtype(float))/10.0;
    return H,x,y,z,deltax,deltay,deltaz;"""

def create_H_from_xyz(x,y,z,orbitals,onsite,thop,d_bond,Nbond):
    # WE ASSUME THAT:
    #
    # 1) TRANSPORT IS IN THE Z DIRECTION
    # 2) THE STRUCTURE IS COMPOSED BY THE SAME TYPE OF ATOMS
    # 3) ALONG THE Z-DIRECTION THE STRUCTURE IS PERIODIC WITH PERIOD EQUAL TO 4 SLICES
    #

    # I find the minimum and maximum coordinates at the border
    # so to take care of the passivation of the atoms at the borders
    xmin=min(x);
    xmax=max(x);
    ymin=min(y);
    ymax=max(y);
    zmin=min(z);
    zmax=max(z);

    # I compute the number of slices (ASSUMPTION 2)
    Nc=int(size(unique(z)));
    #  I have already computed n at the beginning
    #    n=int(size(nonzero(z==zmin)));
    # I compute the number of atoms in the first 4 slices
    temp=unique(z);
    Natom_slices=size(nonzero(z<=temp[3]));
    del temp;

    # I check the maximum number of atoms on each slice;
    u=unique(z);
    Nuz=size(u);
    n=-1;
    for i in range(0,Nuz):
        nnew=size(nonzero(z==u[i]));
        if (nnew>=n):
            n=nnew;
    del i;

    # Now I start doing though stuff
    # I fill x,y and z with dummy atoms
    # If it is a dummy atom, the coordinate is equal to dummy_coord
    
    dummy_coord=10000;
    xa=[];
    ya=[];
    za=[];
    k=0;
    for i in range(0,Nuz):
#        print ya
        nnew=size(nonzero(z==u[i]));
        for j in range(0,nnew):
            xa.append(x[k]);
            ya.append(y[k]);
            za.append(z[k]);
            k=k+1;
        if (nnew<n):
            for j in range(nnew,n):
                xa.append(dummy_coord);
                ya.append(dummy_coord);
                za.append(dummy_coord);
#                k=k+1;
            
    del x,y,z,u,i
    x=array(xa,dtype(float));
    y=array(ya,dtype(float));
    z=array(za,dtype(float));    
#    del xa,ya,za

    Np=size(x);
    Ncol_max=10;
    NN=zeros((Np,Ncol_max),dtype(int));
    border=[]
    # I first find the Nearest Neighbour
    for i in range(0,Np):
        ind=nonzero((sqrt((x-x[i])**2+(y-y[i])**2+(z-z[i])**2)<=d_bond)&(sqrt((x-x[i])**2+(y-y[i])**2+(z-z[i])**2)>1e-10))[0];
        if (size(ind)>Ncol_max):
            print()
            writeout("ERROR IN create_H_from_xyz subroutine in NanoTCAD_ViDES.py file")
            writeout("Use a larger value for Ncol_max")
            print()
            exit(0);
#        print i
        NN[i,0]=i+1;
        NN[i,1:size(ind)+1]=ind+1;
        NPV=size(nonzero(NN[i,:]))-1;
        if (NPV<Nbond):
            border.append(i);
        
        
    # Now I work on the Hamiltonian
    atoms=0;
    i=0;
    h=[];

    # I fill the h list with the number of orbitals
    ll=[orbitals,0];
    fill=zeros(orbitals**2);
    h.append(ll+list(fill))
    del ll,i

    


    # I take care of the diagonal elements
    for i in range(0,Np):
        if ((x[i]<dummy_coord)):
            if (orbitals>1):
                # (ASSUMPTION 1)
                if i in border:
                    xfn=zeros(4);
                    yfn=zeros(4);
                    zfn=zeros(4);
                    if (z[i]==zmin):
                        NPV=size(nonzero(NN[i+4*n,:]))-1;
                        xfn=x[NN[i+n*4,1:NPV+1]-1];
                        yfn=y[NN[i+n*4,1:NPV+1]-1];
                        zfn=z[NN[i+n*4,1:NPV+1]-1];
                        xp=x[i+n*4];
                        yp=y[i+n*4];
                        zp=z[i+n*4];
                    elif (z[i]==zmax):
                        NPV=size(nonzero(NN[i-4*n,:]))-1;
                        xfn=x[NN[i-n*4,1:NPV+1]-1];
                        yfn=y[NN[i-n*4,1:NPV+1]-1];
                        zfn=z[NN[i-n*4,1:NPV+1]-1];
                        xp=x[i-n*4];
                        yp=y[i-n*4];
                        zp=z[i-n*4];
                    else:
                        NPV=size(nonzero(NN[i,:]))-1;
                        xfn=x[NN[i,1:NPV+1]-1];
                        yfn=y[NN[i,1:NPV+1]-1];
                        zfn=z[NN[i,1:NPV+1]-1];
                        xp=x[i];
                        yp=y[i];
                        zp=z[i];

                    deltae=20.0;
                    tempM=Sipassivation(xp,yp,zp,NPV,xfn,yfn,zfn,deltae);
    #                print tempM
    #                print x[i],y[i],z[i]
    #                print xfn
    #                print yfn
    #                print zfn
    #                exit(0);
                    B=zeros((10,10));
                    B[:4,:4]=tempM.reshape(4,4);
                    h.append([i+1,i+1]+list((diag(onsite)+B).flatten()));
    #                h.append([i+1,i+1]+list((diag(onsite)).flatten()));
                    del B,tempM,xfn,yfn,zfn;
                else:
                    h.append([i+1,i+1]+list((diag(onsite)).flatten()));
            else:
                h.append([i+1,i+1]+list(fill));
        else:
            # If the atom is dummy then I mark it with the 77777 value
            # Right now it works only for one orbital
            h.append([i+1,i+1]+list(77777*ones(orbitals**2)));


    # I take care of the off-diagonal elements
    for i in range(0,Np):
        NPV=size(nonzero(NN[i,:]))-1;
        for j in range(0,NPV):
            a1=int(NN[i,0]);
            if (a1<int(NN[i,j+1])):
                if (orbitals>1):
                    # I compute the cosine
                    module=sqrt(((double(x[a1-1])-double(x[int(NN[i,j+1])-1]))**2)+(double(y[a1-1])-double(y[int(NN[i,j+1])-1]))**2+(double(z[a1-1])-double(z[int(NN[i,j+1])-1]))**2);
                    cosx=(-double(x[a1-1])+double(x[int(NN[i,j+1])-1]))/module;
                    cosy=(-double(y[a1-1])+double(y[int(NN[i,j+1])-1]))/module;
                    cosz=(-double(z[a1-1])+double(z[int(NN[i,j+1])-1]))/module;
#                    print a1,int(NN[i,j+1]),cosx,cosy,cosz,module
#                    input=hstack((array([cosx,cosy,cosy]),thop));
#                    print input
#                    matrix_thop=Simatrix(input);
                    matrix_thop=Simatrix(cosx,cosy,cosz,thop);
#                    print matrix_thop
#                    print "----------------"
                    h.append([a1,int(NN[i,j+1])]+list(matrix_thop));
                else:
                    h.append([a1,int(NN[i,j+1]),thop])
    H=array(h,dtype=complex);
    return H,n,Nc,x,y,z;

def get_xyz_from_file(filename):
    fp=open(filename,"r");
    xa=[]
    ya=[]
    za=[]
    for line in fp:
        if (size(line.split())>3):
            xa.append((line.split())[1]);
            ya.append((line.split())[2]);
            za.append((line.split())[3]);
    x=array(xa,dtype(float));
    y=array(ya,dtype(float));
    z=array(za,dtype(float));
    del xa,ya,za
    return x,y,z;


def convert_pdb(filename,orbitals,thop):
    # ASSUMPTION: ALL THE ATOMS ARE OF THE SAME MATERIAL
    
    # I first read the atoms coordinates
    hh=[];
    deltax=0;
    deltay=0;
    deltaz=0;
    x=[];
    y=[];
    z=[];
    i=0;
    fp=open(filename,"r");
    for line in fp:
        hh.append(line);
        if (((hh[i].split()).count('HETATM')==1)|((hh[i].split()).count('ATOM')==1)):
#            ATOM_TYPE=(hh[i].split())[2];
            x.append((hh[i].split())[5]);
            y.append((hh[i].split())[6]);
            z.append((hh[i].split())[7]);
        i=i+1;
    fp.close()
    del hh;

    # Now I work on the Hamiltonian
    hh=[];
    atoms=0;
    i=0;
    h=[];

    # I fill the h list with the number of orbitals
    ll=[orbitals,0];
    fill=zeros(orbitals**2);
    h.append(ll+list(fill))
    del ll

    # I fill the rest of the h list
    fp=open(filename,"r");
    for line in fp:
        hh.append(line);
        atoms=atoms+(hh[i].split()).count('HETATM');
        if (((hh[i].split()).count('HETATM')==1)|((hh[i].split()).count('ATOM')==1)):
             if (orbitals>1):
                 h.append([int((hh[i].split())[1]),int((hh[i].split())[1])]+list((diag(onsite)).flatten()));
             else:
                 h.append([int((hh[i].split())[1]),int((hh[i].split())[1])]+list(fill));
        if ((hh[i].split()).count('CONECT')==1):
            a=(hh[i].split());
            NPV=size(a)-1
            for j in range(0,NPV):
                a1=int(a[1]);
                if (a1<int(a[j+1])):
                    if (orbitals>1):
                        # I compute the cosine
                        module=sqrt(((double(x[a1-1])-double(x[int(a[j+1])-1]))**2)+(double(y[a1-1])-double(y[int(a[j+1])-1]))**2+(double(z[a1-1])-double(z[int(a[j+1])-1]))**2);
                        cosx=(double(x[a1-1])-double(x[int(a[j+1])-1]))/module;
                        cosy=(double(y[a1-1])-double(y[int(a[j+1])-1]))/module;
                        cosz=(double(z[a1-1])-double(z[int(a[j+1])-1]))/module;
                        cosx=1;cosy=1;cosz=1;
                        input=hstack((array([cosx,cosy,cosy]),thop));
                        matrix_thop=Simatrix(input);
                        h.append([a1,int(a[j+1])]+list(matrix_thop));
                    else:
                        h.append([a1,int(a[j+1]),thop])

        if ((hh[i].split()).count('CRYST1')==1):
            a=(hh[i].split());
            if (double(a[1])>=100):
                deltax=0.0;
            else:
                deltax=double(a[1])/10.0;
            if (double(a[2])>=100):
                deltay=0.0;
            else:
                deltay=double(a[2])/10.0;
            if (double(a[3])>=100):
                deltaz=0.0;
            else:
                deltaz=double(a[3])/10.0;            


        i=i+1;
    fp.close()

    H=array(h,dtype(complex));
    x=array(x,dtype(float))/10.0;
    y=array(y,dtype(float))/10.0;
    z=array(z,dtype(float))/10.0;
    return H,x,y,z,deltax,deltay,deltaz;

def Hamiltonian_per(H,x,y,z,deltax,deltay,deltaz,aCC,thop,k):
    Np=size(x);
    Hnew=H.copy();
    conn_per=[]
    for ii in range(0,Np):
        xc=x[ii];
        yc=y[ii];
        zc=z[ii];
        # Here I compare with 1.05*aCC in order to take into account numerical tollerances
        indp=nonzero(sqrt((x-xc+deltax)**2+(y-yc+deltay)**2+(z-zc+deltaz)**2)<aCC*1.05)[0]+1;
        indm=nonzero(sqrt((x-xc-deltax)**2+(y-yc-deltay)**2+(z-zc-deltaz)**2)<aCC*1.05)[0]+1;
        if (size(indp)>0):
            for j in range(0,size(indp)):
                conn_per.append([ii+1,indp[j]]);
        if (size(indm)>0):
            for j in range(0,size(indm)):
                conn_per.append([ii+1,indm[j]]);


    del ii
    Nconn=len(conn_per);
    for ii in range(Nconn):
        ind=nonzero((H[:,0]==conn_per[ii][0])&(H[:,1]==conn_per[ii][1]))[0]
        if (size(ind)>0):
            if (deltax>0):
                segno=sign(x[int(abs(H[ind,0]))-1]-x[int(abs(H[ind,1]))-1]);
                Hnew[ind,2]=H[ind,2]+thop*exp(-segno*k*deltax*1j);
            elif (deltay>0):
                segno=sign(y[int(abs(H[ind,0]))-1]-y[int(abs(H[ind,1]))-1]);
                Hnew[ind,2]=H[ind,2]+thop*exp(-segno*k*deltay*1j);
            else:
                segno=sign(z[int(abs(H[ind,0]))-1]-z[int(abs(H[ind,1]))-1]);
                Hnew[ind,2]=H[ind,2]+thop*exp(-segno*k*deltaz*1j);
        else:
            if (conn_per[ii][0]<conn_per[ii][1]):
                if (deltax>0):
                    segno=sign(x[conn_per[ii][0]-1]-x[conn_per[ii][1]-1]);
                    temp=array([conn_per[ii][0],conn_per[ii][1],thop*exp(-segno*k*deltax*1j)]);
                elif (deltay>0):
                    segno=sign(y[conn_per[ii][0]-1]-y[conn_per[ii][1]-1]);
                    temp=array([conn_per[ii][0],conn_per[ii][1],thop*exp(-segno*k*deltay*1j)]);
                else:
                    segno=sign(z[conn_per[ii][0]-1]-z[conn_per[ii][1]-1]);
                    temp=array([conn_per[ii][0],conn_per[ii][1],thop*exp(-segno*k*deltaz*1j)]);
                Hnew=vstack([Hnew,temp]);
    del ii
    return Hnew


class nanoribbon_fast_ohmic:
    acc=0.144;
    def __init__(self,n,L):
        self.Nc=int(4*(floor((floor(L/nanoribbon_fast_ohmic.acc)-1)/3)));
        self.n=n;
        self.Phi=zeros(n*self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.dE=1e-3;
        self.thop=-2.7;
        self.eta=1e-8;
        self.mu1=0;
        self.mu2=0;
        self.Temp=300;
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(self.n*self.Nc);
        self.rank=0;
        self.atoms_coordinates();
        self.defects_list=[]
        self.onsite_E=-1.5;
    def atoms_coordinates(self):
        GNR_atoms_coordinates(self);
        self.x=array(self.x);
        self.y=array(self.y);
        self.z=array(self.z);
        return;
    def gap(self):
        return GNRgap(self);
    def charge_T(self):
        
        M=self.Nc;
        N=self.n;
        t=self.thop;
        Energy = 0.0
        Ene = 0.0
        p = 0.0
        d = 0.0
        orbitals = [1, 0]
        hamiltonian = []
        zeroes = [0, 0, 0, 0]
        ene = [Energy, 0, 0, Ene]
        coupling1 = [t, 0, 0, p]
        coupling2 = [t*1.12, 0, 0, p]
        orbitals = orbitals + zeroes
        hamiltonian.append(orbitals)

        for j in range(M):
            for i in range(N):
                n = i + 1 + j*N
                p = [n,n]
                p =  p + ene
                hamiltonian.append(p)

        

        for j in range(1, M-1, +4):
            for i in range(1, N):
                n = i + 1 + j*N
                m = i + (j+1)*N 
                p = [n,m]
                p =  p + coupling1
                hamiltonian.append(p)
          #      hamiltonian.append([m, n, t, p, d])

        for j in range(3, M-1, +4):
            for i in range(0, N-1):
                n = i + 1 + j*N
                m = i + 2 + (j+1)*N 
                p = [n,m]
                p =  p + coupling1
                hamiltonian.append(p)
           #     hamiltonian.append([m, n, t, p, d])


        # nell'if ripristinare il fattore t*1.12  
        for j in range(0, M-1, +4):
            for i in range(N):
                n = i + 1 + j*N
                m = i + 1 + (j+1)*N
                if i == 0:
                      p = [n,m]
                      p =  p + coupling2
                      hamiltonian.append(p)
            #          hamiltonian.append([m, n, t*1.12, p, d])
                else :
                      p = [n,m]
                      p =  p + coupling1
                      hamiltonian.append(p)
             #         hamiltonian.append([m, n, t, p, d])

        for j in range(1, M-1, +4):
            for i in range(N):
                n = i + 1 + j*N
                m = i + 1 + (j+1)*N 
                p = [n,m]
                p =  p + coupling1
                hamiltonian.append(p)
#            hamiltonian.append([m, n, t, p, d])

        # nell'if ripristinare il fattore t*1.12 
        for j in range(2, M-1, +4):
            for i in range(N):
                n = i + 1 + j*N
                m = i + 1 + (j+1)*N
                if i == (N-1):
                      p = [n,m]
                      p =  p + coupling2
                      hamiltonian.append(p)
                #                 hamiltonian.append([m, n, t*1.12, p, d])
                else :
                      p = [n,m]
                      p =  p + coupling1
                      hamiltonian.append(p)
#                hamiltonian.append([m, n, t, p, d])

        for j in range(3, M-1, +4):
            for i in range(N):
                n = i + 1 + j*N
                m = i + 1 + (j+1)*N 
                p = [n,m]
                p =  p + coupling1
                hamiltonian.append(p)
#         hamiltonian.append([m, n, t, p, d])

        H = Hamiltonian(N,M)
        # I work on the defects
        ind=array(self.defects_list,dtype=int);
        H.H=array(hamiltonian,dtype=complex)
        H.H[ind,2]=self.onsite_E;
        H.Eupper = self.Eupper;
        H.Elower = self.Elower;
        H.rank=self.rank;
        H.dE=self.dE;
        H.Phi=self.Phi;
        H.Ei=-self.Phi;
        H.eta=self.eta;
        H.mu1=self.mu1;
        H.mu2=self.mu2;
        H.Egap=self.gap();

        H.charge_T()

        self.E=array(H.E);
        self.T=array(H.T);
        self.charge=array(H.charge);
        del hamiltonian,H
        return;

    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);
        

# This is the class for the solution of the 1D drift-diffusion
class multisubband1D:
    def __init__(self, nx, ny, Neig):
        self.ny=ny;
        self.nx=nx;
        self.x=zeros(nx);
        self.y=zeros(ny);
        self.Phi=zeros(nx*self.ny);
        self.Ei=zeros(nx*self.ny);
        self.Egap=zeros(nx*self.ny);
        self.Temp=300;
        self.charge=zeros(nx*self.ny);
        self.rank=0;
        self.Neig=Neig;
        self.Psi=zeros((nx*ny,Neig));
        self.eig=zeros((ny,Neig));
        self.mass=zeros((nx,ny));
        self.mu=100e-4*ones(self.ny);
        self.genric=zeros(self.ny);
        self.n1d=zeros(self.ny);
        self.ecs=zeros(self.ny);
        self.charge_left_contact=0;
        self.charge_right_contact=0;
        self.tolljay=1e-3;

# This is the class for the solution of the QM 1D 
class QM1D:
    def __init__(self, nx, Neig,gridx,p=None,charge_T=None):
        if charge_T is not None:
            self.charge_T=types.MethodType(charge_T,self);
        self.nx=nx;
        self.x=zeros(nx);
        self.ny=1;
        ny=1;
        self.Phi=zeros(nx*self.ny);
        self.Ei=zeros(nx*self.ny);
        self.Temp=300;
        self.charge=zeros(nx*self.ny);
        self.rank=0;
        self.Neig=Neig;
        self.Psi=zeros((nx*ny,Neig));
        self.eig=zeros((ny,Neig));
        if p is not None:
            self.Egap=p.Egap;
            self.massl=p.mel
            self.masst=p.met;
            self.massh=p.mhole
            self.chi=p.chi
            self.mass=p.mel;
        else:
            self.Egap=zeros(nx*self.ny)
            self.massl=zeros(nx*self.ny)
            self.masst=zeros(nx*self.ny)
            self.massh=zeros(nx*self.ny)
            self.chi=zeros(nx*self.ny)
            self.mass=zeros(nx*self.ny)
        self.Ef=0;
        self.x=gridx;
        self.ecs=zeros(self.ny);
    def charge_T(self):
        del self.charge
        self.charge=zeros(self.nx*self.ny);
        self.Ei=-self.Phi;
        # I compute the confined electrons
        dist=avervect(self.x)
#        self.Ei=4.05-self.Phi-self.chi-self.Egap*0.5
        self.mass=self.massl;
        solve_schroedinger_1D(self);
        vt=self.Temp*kboltz/q;
        for i in range(0,self.Neig):
            self.charge=self.charge-2*dist*1e-9*(self.Psi[:,i])**2*self.masst*m0*kboltz*self.Temp/pi/hbar**2*log(1+exp(-(self.eig[0,i]-self.Ef)/vt));
        self.mass=self.masst;
        solve_schroedinger_1D(self);
        vt=self.Temp*kboltz/q;
        for i in range(0,self.Neig):
            self.charge=self.charge-4*dist*1e-9*(self.Psi[:,i])**2*self.massl*m0*kboltz*self.Temp/pi/hbar**2*log(1+exp(-(self.eig[0,i]-self.Ef)/vt));   
        # I now add the holes
        for i in range(0,size(self.charge)):
            self.charge[i]=self.charge[i]+dist[i]*1e-9*(2/sqrt(pi))*2*(vt/(2*pi)*(self.massh[i]*m0/hbar)*(q/hbar))**1.5*fphalf((self.Ei[i]-self.Egap[i]*0.5-self.Ef)/vt)

        return;
    def current(self):
        return 0;
