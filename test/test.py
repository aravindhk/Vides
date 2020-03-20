import subprocess
import os
from numpy import *
import sys

def default_version():
    c="python -c 'import sys; v=sys.version[:3]; print(v)'> python_version_temp"
    os.system(c)
    fp=open("python_version_temp","r")
    v=fp.read();
    fp.close()
    c="rm python_version_temp"
    os.system(c)
    v=v[:3]
    return v;

dv=default_version();

if (sys.version[:3]==dv):
    pitone="python"
else:
    pitone="python%s" %sys.version[:3]

print("\n\n")
print("**********************************")
print("        NanoTCAD ViDES test")
print("**********************************")
print("\n\n")

print("-------------------------------")
print("      Starting test on GNR")
print("-------------------------------")

c="%s ./script_GNR_auto.py" %pitone
os.system(c)

fp=open("GNR_test");
string=fp.read()
print(string)
fp.close()

print("-------------------------------")
print("   Starting test on graphene")
print("-------------------------------")
c="%s ./test_demetrio_graphene_auto.py" %pitone
os.system(c)

fp=open("graphene_test");
string2=fp.read()
print(string2)
fp.close()

print("-------------------------------")
print("   Starting test on bilayer")
print("-------------------------------")
c="%s ./test_bilayer_class_auto.py" %pitone
os.system(c)

fp=open("bilayer_test");
string3=fp.read()
print(string3)
fp.close()


print("-------------------------------")
print("   Starting test on nanowire")
print("-------------------------------")
print("Zero electric field")
c="%s ./test_SNW.py" %pitone
os.system(c)

fp=open("nanowire_test");
string4=fp.read()
print(string4)
fp.close()

string5="\n \n --------------------- \n %s --------------------- \n %s --------------------- \n %s --------------------- \n TEST on nanowire \n \n TEST on nanowire with non-zero electric field \n %s" %(string,string2,string3,string4);

fp=open("test.log","w");
fp.write(string5);
fp.close()


