#!/usr/bin/python
from NanoTCAD_ViDES import *
from subprocess import call
from os import chdir


# L = a*(N-1)/4
# D = (a*sqrt(2)/2)*N
 
rank = 0


fields = [0.1]


#Testing the cicle with the smallest nanowires
#SI
parameters1 = [[3.85, 7.85],[6.00, 4.00],[1.00, 0.74]]
#InAs
parameters2 = [[4.29, 7.85],[6.00, 4.00],[1.00, 0.74]]
#Ge
parameters3 = [[4.00, 7.85],[6.00, 4.00],[1.00, 0.74]]


fp=open("nanowire_test","w");
materialrange=['Si','Ge','InAs'];
for material in materialrange:

    if (material=="Si"):
        parameters=parameters1;
    elif (material=="InAs"):
        parameters=parameters2;
    elif (material=="Ge"):
        parameters=parameters3;

    for index in range(0,1,+1):
        for Field in fields:
            Egap = parameters[1][index]
            Ei = parameters[2][index]
            #Voltage = Egap + 0.3
            Voltage = 0.0
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
            D = parameters[0][index]
#            L = (a0/4)*(layers-1)



            H = Zincblend(material, S, Or, D, L)

            print H.sqci, H.tilt, H.edge, H.zmax


            deltaE = 0.1

            V = (Voltage)/2
            H.Elower = 0
            H.Eupper = 10
            H.dE = deltaE
            H.eta = 1e-8
            n = H.n
            Nc = H.Nc

            print H.n, H.Nc, n, Nc


            print (H.Nc_aux*H.n_aux)*((H.Nc_aux*H.n_aux+1)/2)*(2+100)

            print H.max[0]


            H.mu1 = V
            H.mu2 = -V


            length = H.z[(Nc-1)*(n) + (n-1)]



            effF = round((2*V)/(length-2*H.z[4*n]),3)

            print effF

            for i in range (4):
                for j in range(n):
                    k = i*n + j
                #Constant potential, test for correct energy shift
                #H.Phi[k] = V 
                #Linear potential drop
                #H.Phi[k] = H.mu1*(1-2*(H.z[k]/length))
                    H.Phi[k] = H.mu1
                #H.Phi[k] = H.mu1*(1-2*(H.z[k]/length))
                #           plot(H.z,H.Phi)   
                #          show()

            for i in range (Nc-4, Nc):
                for j in range(n):
                    k = i*n + j
                    H.Phi[k] = H.mu2

            for i in range (4, Nc-4):
                for j in range(n):
                    k = i*n + j
                    H.Phi[k] = H.mu1*(1-2*((H.z[k]-H.z[4*n])/(length-2*H.z[4*n])))


            H.charge_T()       

            string="transmission-test-%s.dat" %material;
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


#            print H.T
#            print a[:,1];
            
#            show()
fp.close()

            

#            if (rank==0):
#                pippo=[H.E,H.T]
#                savetxt("transmission.dat", transpose(pippo))     
#                c="cp transmission.dat transmission-test.dat" 
#                commands.getoutput(c);
 
        

