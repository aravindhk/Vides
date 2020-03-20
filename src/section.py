# ======================================================================
#  Copyright (c) 2010, G. Fiori, University of Pisa
#  
#  This file is released under the BSD license.
#  See the file "license.txt" for information on usage and
#  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# ====================================================================== 

from numpy import *
import os
import sys
if sys.version > '3':
    import subprocess;
else:
    import subprocess

def section(slicedir,quantity,coordslice,grid):
    if (slicedir=="x"):
        # I find the index of the closest 
        # coordinate to coordslice
        index=nonzero(abs(grid.gridx-coordslice)==
                      min(abs(grid.gridx-coordslice)))[0][0]+1;
        count=0;
        swapx=argsort(grid.gridx);
        swapy=argsort(grid.gridy);
        swapz=argsort(grid.gridz);
        swapx3D=meshgrid(meshgrid(swapx,swapy)[0].flatten(),swapz)[0].flatten();
        swapy3D=meshgrid(meshgrid(swapx,swapy)[1].flatten(),swapz)[0].flatten();
        swapz3D=meshgrid(meshgrid(swapx,swapy)[1].flatten(),swapz)[1].flatten();    
        indexsort=swapx3D+swapy3D*grid.nx+swapz3D*grid.nx*grid.ny;
        quantitys=quantity[indexsort];
        out=zeros(grid.ny*grid.nz);
        yy=zeros(grid.ny*grid.nz);
        zz=zeros(grid.ny*grid.nz);
        fp=open("newplot","w");
        for k in range(0,grid.nz):
            for j in range(0,grid.ny):
                ix=index+j*grid.nx+k*grid.nx*grid.ny;
                fp.write("%s %s %s \n" %(grid.gridy[swapy[j]],grid.gridz[swapz[k]],
                                         quantitys[ix]))
            fp.write("\n");
        fp.close();

        fp=os.popen("gnuplot","w")
        fp.write("set title 'SECTIONX' \n")
        fp.write("set mouse \n")
        fp.write("splot 'newplot' w l \n")
        fp.flush()
        eval(input("waiting ...... Press any key to continue"))
        fp.close()
    elif (slicedir=="y"):
        # I find the index of the closest 
        # coordinate to coordslice
        index=nonzero(abs(grid.gridy-coordslice)==
                      min(abs(grid.gridy-coordslice)))[0][0];
        count=0;
        out=zeros(grid.nx*grid.nz);
        xx=zeros(grid.nx*grid.nz);
        zz=zeros(grid.nx*grid.nz);
        fp=open("newplot","w");
        for k in range(0,grid.nz):
            for i in range(0,grid.nx):
                ix=i+index*grid.nx+k*grid.nx*grid.ny;
                fp.write("%s %s %s \n" %(grid.gridx[i],grid.gridz[k],
                                         quantity[ix]))
            fp.write("\n");
        fp.close();

        fp=os.popen("gnuplot","w")
        fp.write("set title 'SECTIONY' \n")
        fp.write("set mouse \n")
        fp.write("splot 'newplot' w l \n")
        fp.flush()
        eval(input("waiting ...... Press any key to continue"))
        fp.close()
    elif (slicedir=="z"):
        # I find the index of the closest 
        # coordinate to coordslice
        index=nonzero(abs(grid.gridz-coordslice)==
                      min(abs(grid.gridz-coordslice)))[0][0];
        count=0;
        out=zeros(grid.nx*grid.ny);
        xx=zeros(grid.nx*grid.ny);
        yy=zeros(grid.nx*grid.ny);
        fp=open("newplot","w");
        for j in range(0,grid.ny):
            for i in range(0,grid.nx):
                ix=i+j*grid.nx+index*grid.nx*grid.ny;
                fp.write("%s %s %s \n" %(grid.gridx[i],grid.gridy[j],
                                         quantity[ix]))
            fp.write("\n");
        fp.close();

        fp=os.popen("gnuplot","w")
        fp.write("set title 'SECTIONZ' \n")
        fp.write("set mouse \n")
        fp.write("splot 'newplot' w l \n")
        fp.flush()
        eval(input("waiting ...... Press any key to continue"))
        fp.close()
    return;
        
            
def plot2D(quantity,grid,filename):

    fp=open(filename,"w");
    for j in range(0,grid.ny):
        for i in range(0,grid.nx):
            s="%s %s %s \n" %(grid.gridx[i],grid.gridy[j],quantity[i+j*grid.nx]);
            fp.write(s);
        s="\n";
        fp.write(s);
    fp.close();

    fp=os.popen("gnuplot","w")
    fp.write("set mouse \n")
    string="splot '%s' w l \n" %filename
    print(string)
    fp.write(string)
    fp.flush()
    eval(input("waiting ...... Press any key to continue"))
    fp.close()
