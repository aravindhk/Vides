#!/usr/bin/python
from NanoTCAD_ViDES import *
from subprocess import call
from os import chdir


# L = a*(N-1)/4
# D = (a*sqrt(2)/2)*N
 
rank = 0


fields = [1]


#Testing the cicle with the smallest nanowires
#SI
parameters1 = [3.85,6.00,1.00]
#InAs
parameters2 = [4.29,4.465433,1.00]
#Ge
parameters3 = [4.00,5.2,1.07]

rangematerial = ['Si','InAs','Ge' ]

#rangematerial = ['Si', 'InAs']

fp=open("nanowire_test_field","w");
for material in rangematerial:

    

        if material == 'Si':
            parameters = parameters1
        elif material == 'InAs':
            parameters = parameters2
        elif material == 'Ge':
            parameters = parameters3

        for Field in fields:
            Egap = parameters[1]
            Ei = parameters[2]
            Voltage = Egap + 0.3
            #Voltage = 0.0
#            if material == 'Si':
#                a0 = 5.431
#            if material == 'Ge':
#                a0 = 5.6575
#            if material == 'InAs':
#                a0 = 6.0583

            L = Voltage/(Field)
#            layers = int(4*L/(a0) + 1)
#
#            if (rank==0):
#                print "prima=", layers
#
#            if layers%4==1:
#                layers-=1
#            elif layers%4==2:
#                layers-=2
#            elif layers%4==3:
#                layers+=1
#
#            print layers
#
#            if layers%4!=0:
#                print "INTERRUPT AT WIRE", material, parameters[0][i]
#                print "NUMBER OF SLICES NOT MULTIPLE OF 4"
#                quit()
#
#            layers += 8


            S = 0
            Or = 1
            D = parameters[0]
#            L = (a0/4)*(layers-1)



            H = Zincblend(material, S, Or, D, L)

            print H.sqci, H.tilt, H.edge, H.zmax

            print H.n
            print H.Nc
#            quit()


            deltaE = 0.001

            V = (Voltage)/2
            H.Elower = Ei
            H.Eupper = Ei+0.1
	    
            H.dE = deltaE
            H.eta = 1e-8
            n = H.n
            Nc = H.Nc

            print H.n, H.Nc, n, Nc


            print (H.Nc_aux*H.n_aux)*((H.Nc_aux*H.n_aux+1)/2)*(2+100)

            print H.max[0]


            H.mu1 = V
            H.mu2 = -V


     
	    print 'length of the wire from LMAX', H.L

            effF = round((2*V)/(H.L-2*H.a0),3)

	    print 'active length', H.L-2*H.a0
	   

            print 'applied field', effF

	    i = 0;

	    print H.z[:H.n*H.Nc];

#	    exit()

	    #for ln in H.z:
		    
	#	    H.Phi[i]=H.mu1*(1-2*(H.z[i]/(H.L)))
		#    i = i + 1

	    H.Phi=H.mu1*(1-2*((H.z-H.a0)/(H.L)))

	    #print H.Phi[:H.n*H.Nc];
  
		    
#            for i in range (4):
#                for j in range(n):
#                    k = i*n + j
                #Constant potential, test for correct energy shift
                #H.Phi[k] = V 
                #Linear potential drop
                #H.Phi[k] = H.mu1*(1-2*(H.z[k]/length))
#                    H.Phi[k] = H.mu1
                #H.Phi[k] = H.mu1*(1-2*(H.z[k]/length))
                #           plot(H.z,H.Phi)   
                #          show()

#            for i in range (Nc-4, Nc):
#                for j in range(n):
#                    k = i*n + j
#                    H.Phi[k] = H.mu2

#            for i in range (4, Nc-4):
#                for j in range(n):
#                    k = i*n + j
#                    H.Phi[k] = H.mu1*(1-2*((H.z[k]-H.a0)/(H.L-2*H.a0)))
	    
	    H.Phi = ciccione(H.Phi,H.n,H.Nc,H.z,H.a0);
	    
	    H.z = ciccione(H.z,H.n,H.Nc,H.z,H.a0);

#	    print H.Phi[:H.n*H.Nc];

	    for i in range(4*n):
		    H.Phi[i]=H.Phi[0]

            for i in range (Nc-4, Nc):
                for j in range(n):
                    k = i*n + j
                    H.Phi[k] = H.mu2

            for i in range (4, Nc-4):
                for j in range(n):
                    k = i*n + j
                    H.Phi[k] = H.mu1*(1-2*((H.z[i*n]-2*H.a0)/(H.L-2*H.a0)))
		    print H.z[k]

		    #H.Phi[k] = H.mu1*(1-2*((H.z[k]-H.z[4*n])/(length-2*H.z[4*n])))

	    print 'lenght', H.z[(Nc-1)*(n) ]
	    print  H.z[4*n]

	    print H.L

	    print H.z
	

	    print H.Phi;

	    	    
	    #exit()

            print H.charge_T()       
	    print "sono qua dopo python"

            string="transmission-test-%s-field.dat" %material;
            a=loadtxt(string);
            if (max(abs(H.T[:size(a[:,1])]-a[:,1]))<1e-5):
                string="%s PASSED \n" %material;
            else:
                string="%s NOT PASSED \n" %material;
            fp.write(string);

            print (max(abs(H.T[:size(a[:,1])]-a[:,1])))
	    
#            plot(H.E,H.T);
#            hold
#            plot(a[:,0],a[:,1],'o');
#            show()


#            if (rank==0):
#                pippo=[H.E,H.T]
#                savetxt("transmission.dat", transpose(pippo))     
#                c="cp transmission.dat transmission-test-%s.dat" %material
#                commands.getoutput(c);
 
fp.close()
