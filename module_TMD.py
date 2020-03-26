#Copyright (c) 2004-2015, G. Fiori, G. Iannaccone, University of Pisa.
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#    - Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.
#    - Redistributions in binary form must reproduce the above copyright
#        notice, this list of conditions and the following disclaimer in the
#        documentation and/or other materials provided with the distribution.
#    - All advertising materials mentioning features or use of this software
#        must display the following acknowledgement:
#        This product includes software developed by G.Fiori and G.Iannaccone at
#        University of Pisa.
#    - Neither the name of the University of Pisa nor the
#        names of its contributors may be used to endorse or promote products
#        derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY
#EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE AUTHORS
#BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#                       SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#                       INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
#THE POSSIBILITY OF SUCH DAMAGE.


from NanoTCAD_ViDES import *
class TMD:
    
    def __init__(self,L,n_or_p):
        #degeneracy of 2 is already taken into account
        #in the Hamiltonian, since the minimum is in Kf,
        #so by itself degeneracy is 2 as in graphene
        #if you specify deg=2, then total degeneracy is 4
        me=0.572;
        mh=0.659
        self.Egap=1.78
        self.deg=1;
        self.n=1;
        ymin=-20;
        ymax=L+20;
        self.acc=0.2;
        self.Nc=int(4*(floor((floor(L/self.acc)-1)/3)));
        self.Phi=zeros(self.Nc);
        self.Ei=zeros(self.Nc);
        self.Eupper=1000.0;
        self.Elower=-1000.0;
        self.delta=sqrt(3)*self.acc;
        self.kmax=pi/self.delta;
        self.kmin=0;
        self.dk=0.1;
        self.dE=1e-3;
        self.thop=-2.59;
        if (n_or_p=="n"):
            self.thopBN=-sqrt(2*q*self.Egap/(3*(self.acc*1e-9*sqrt(3))**2*m0*me))*hbar/q;
        else:
            self.thopBN=-sqrt(2*q*self.Egap/(3*(self.acc*1e-9*sqrt(3))**2*m0*mh))*hbar/q;
        self.eta=1e-5;
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
        self.ymin=ymin;
        self.ymax=ymax;
        self.BC_MX2=0.2;
        self.BV_MX2=0.2-self.Egap;
        self.atoms_coordinates();
    def atoms_coordinates(self):
        GNR_atoms_coordinates(self);
        self.y=array(self.z);
        self.x=zeros(size(self.y));
        return;
    def gap(self):
        return self.Egap;
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
            if ((self.y[i-1]>=self.ymin)&(self.y[i-1]<=self.ymax)):
                if ((i%2)==0):
                    h[i][2]=self.BC_MX2;
                else:
                    h[i][2]=self.BV_MX2;

        kk=1;
        for ii in range(slices+1,2*slices):
            if ((ii%2)==1):
                h[ii][0]=kk;
                h[ii][1]=kk+1;
                if ((self.y[kk-1]>=self.ymin)&(self.y[kk-1]<=self.ymax)):
                    h[ii][2]=self.thopBN;
                else:
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
            if (self.rank==0): print("----------------------------------")
            if (self.rank==0): print(("    kx range: [%s,%s] ") %(self.kmin,self.kmax));
            if (self.rank==0): print(("    iteration %s ") %i);
            if (self.rank==0): print("----------------------------------")
            flaggo=0;
            kk=1;
            # I fill the Hamiltonian for the actual wavevector k in the cycle
            for ii in range(slices+1,2*slices):
                if ((ii%2)==0):
                    h[ii][0]=kk;
                    h[ii][1]=kk+1;
                    if ((flaggo%2)==0):
                        if ((self.y[kk-1]>=self.ymin)&(self.y[kk-1]<=self.ymax)):
                            h[ii][2]=self.thopBN+self.thopBN*exp(k*self.delta*1j);
                        else:
                            h[ii][2]=self.thop+self.thop*exp(k*self.delta*1j);
                    else:
                        if ((self.y[kk-1]>=self.ymin)&(self.y[kk-1]<=self.ymax)):
                            h[ii][2]=self.thopBN+self.thopBN*exp(-k*self.delta*1j);
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
                self.T=self.deg*H.T*(2*self.dk/(2*pi));
                self.charge=self.deg*H.charge*(2*self.dk/(2*pi));
            else:
                # the factor 2 is because I integrate over kx>0
                self.T=self.T+self.deg*H.T*(2*self.dk/(2*pi));
                self.charge=self.charge+self.deg*H.charge*(2*self.dk/(2*pi));

            if (self.T2D=="yes"):
                print(size(Z[:,i]),size(H.T),size(EE));
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
        #        self.charge=array(self.charge)*1e9;
        temporary=zeros(size(self.charge));
        #Let's give back the average charge within 2 near atoms
        for ii in range(0,self.Nc,2):
            temporary[ii]=(self.charge[ii]+self.charge[ii+1])*0.5*1e9;
            temporary[ii+1]=(self.charge[ii]+self.charge[ii+1])*0.5*1e9;
        self.charge=temporary
        #        savetxt("ncar.temp",self.charge)
        #    exit(0);
        del kvect,h;
        return;

    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);


