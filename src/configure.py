#!/usr/bin/python
import os
import numpy

def cambia(ingresso,uscita,stringain,numero):
    #la funzione legge il file di ingresso "ingresso" e da'
    # in uscita il file "uscita". Ingresso ed uscita sono
    # delle stringhe. La funzione cambia, sostituisce la
    # stringa stringain con il numero (convertito internamente
    # in stringa) numero.
    fp=open(ingresso,"r")
    f2=open(uscita,"w")
    s=fp.read()
    temp='%s'%numero
    ss=s.replace(stringain,temp)
    f2.write(ss)
    f2.close()
    return

############################################################################
# I find the python version
import sys
v=sys.version[:3];

#I find where is Python.h
c="python%s-config --includes>tt" %v
os.system(c);
fp=open("tt","r");
a=fp.read();
aa=a.split()
INCPY=aa[0][2:];

#I find where is arrayobject.h
s1=numpy.get_include();
s2="numpy"
s3="%s/%s" %(s1,s2);
INCNPY=s3
prefix=sys.prefix;

#I find the link to the libraries
c="python%s-config --libs>tt" %v
os.system(c);
fp=open("tt","r");
a=fp.read();
aa=a.split()
LIB=aa[-1][:];

cambia("Makefile.form","M1","{fff}",v);
cambia("M1","M2","{INCPY}",INCPY);
cambia("M2","M3","{INCNPY}",INCNPY);
cambia("M3","M4","{LIB}",LIB);
cambia("M4","Makefile","{prefix}",prefix);

c="rm tt;rm M1;rm M2;rm M3;rm M4"
os.system(c);
