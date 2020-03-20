CC	=	gcc
FF	=	gfortran
CFLAGS	=	-fPIC -O3
FFLAGS	=	-fPIC -O3 -lpython2.6
NOMAIN  =       
CCMPI   =	mpicc
FFMPI   =	mpif77
INCPY	=	/usr/include/python2.6
INCNPY	= 	/usr/lib/pymodules/python2.6/numpy/core/include/numpy

OBJECTS =	nrutil.o nonuniformgridmod.o complass.o cmatrix.o cvectorm.o max.o create_updiag.o create_lowdiag.o cfree_cmatrix.o cfree_cvectorm.o min.o selfanalitical.o selfschottky.o LDOS.o create_beta1.o create_beta2mode.o create_beta2transpmode.o VAVdaga.o ccmul.o ccadd.o ccsub.o cdmul.o ccsqrt.o ccdiv.o cIm.o cRe.o cmatsub.o rgfblock.o spectralfun.o cfree_ctensor4.o transmission.o cmatdaga.o Fermi_Dirac.o csub.o csum.o ctensor4.o cmatinv.o cmatmul.o cmatsum.o cdabs.o cmatrixm.o ccvector.o cfree_ccvector.o cmatmul3.o cmatRe.o create_beta2.o create_beta2transp.o NanoTCAD_ViDESmod.o LDOSmode.o create_updiagmode.o create_lowdiagmode.o selfanaliticalmode.o VdagaAV.o energysort.o spectralfunmode.o dsort.f create_updiagGNR.o create_lowdiagGNR.o selfGNR.o gappo.o create_beta1GNR.o create_beta2GNR.o create_beta2transpGNR.o Gzerozero.o cmatnorm2diff.o  nanotube.o inversion.o create_J.c ffunpoisson.o preparej.o Jbuild.o domn.o abbs.o indice.o norma2.o destroy_J.o stampaout.o graphenedown.o

OBJECTS_REL_1.1 =	 Hamiltonian_py.o H_charge_T.o selfH_W.o GzerozeroH_W.o zeig_S.o slamch.o csscal.o cswap.o ccopy.o cdotc.o cdotu.o cgemm.o cgemv.o  cscal.o  ctrmm.o ctrmv.o icamax.o scnrm2.o  caxpy.o cgerc.o scabs1.o ilaclc.o ilaclr.o icmax1.o zgeev.o isamax.o scasum.o sscal.o ctrsv.o zdscal.o dlamch.o idamax.o dznrm2.o zcopy.o zaxpy.o ilazlc.o ilazlr.o zgerc.o dzasum.o dscal.o zdotu.o zdotc.o ztrsv.o zggev.o zhgeqz.o zgghrd.o ztgevc.o zggbak.o zlanhs.o zggbal.f zgeqrf.o zgeqr2.o zlarfp.o zgges.o ztgsen.o zlacn2.o ztgexc.o ztgex2.o ztgsyl.o ztgsy2.o zlatdf.o zgetc2.o zgecon.o zdrscl.o dzsum1.o zgesc2.o izmax1.o ConvertPycomplex_to_Ccomplex.o phonon_GNR.o free_dvectorm.o free_ivectorm.o ceig.o phonon_graphene.o electron_GNR.o Bose_Einstein.o rateACABS.o rateACEM.o rateOPTABS.o rateOPTEM.o zhpevx.o solvePoisson2D.o ffunpoisson2D.o preparej2D.o Jbuild2D.o indice2D.o create_J2D.o destroy_J2D.o sign.o c3tensor.o cfree_c3tensor.o

OBJECTS_REL_1.2 =	 Zinc.o nanowire.o hamiltonian.o hamil.o

all:	$(OBJECTS) $(OBJECTS_REL_1.1) $(OBJECTS_REL_1.2) NanoTCAD_ViDESmod

$(OBJECTS) : 	

NanoTCAD_ViDESmod   : $(OBJECTS) $(OBJECTS_REL_1.1) $(OBJECTS_REL_1.2) NanoTCAD_ViDESmod.o
	$(FF) $(CFLAGS) -shared -o NanoTCAD_ViDESmod.so $(OBJECTS) $(OBJECTS_REL_1.1) $(OBJECTS_REL_1.2) -lm $(NOMAIN) -lpython2.6

comp  :
	rm NanoTCAD_ViDESmod.o
	make

clean :		
	rm -f *.o
	rm -f *.so	
	rm -f core

# I make a check if I want to compile the parallelized code
ifeq ($(MAKECMDGOALS),ViDESMPI)
.f.o  : 
	$(FFMPI) $(FFLAGS) -c $<
else
.f.o  :
	$(FF) $(FFLAGS) -I$(INCPY) -I$(INCNPY) -c $<
endif

# I make a check if I want to compile the parallelized code
ifeq ($(MAKECMDGOALS),ViDESMPI)
.c.o  : 
	$(CCMPI) $(CFLAGS) -c $<
else
.c.o  :
	$(CC) $(CFLAGS) -I$(INCPY) -I$(INCNPY) -c $<
endif

install : NanoTCAD_ViDESmod
	cp NanoTCAD_ViDES.py /usr/lib/python2.6
	cp NanoTCAD_ViDESmod.so /usr/lib/python2.6
	cp section.py /usr/lib/python2.6

uninstall : NanoTCAD_ViDESmod
	  rm /usr/lib/python2.6/NanoTCAD_ViDES.py
	  rm /usr/lib/python2.6/NanoTCAD_ViDESmod.so
	  rm /usr/lib/python2.6/section.py

test : 
	cp NanoTCAD_ViDES.py ../test
	cp NanoTCAD_ViDESmod.so ../test
	cp section.py ../test
	cd ../test;python test.py
	rm ../test/NanoTCAD_ViDES.py
	rm ../test/NanoTCAD_ViDESmod.so
	rm ../test/section.py
