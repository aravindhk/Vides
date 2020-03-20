c ======================================================================
c  Copyright (c) 2004-2008,  G. Fiori, G. Iannaccone  University of Pisa
c
c  This file is released under the BSD license.
c  See the file "license.txt" for information on usage and
c  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
c ====================================================================== 

	subroutine nonuniformgridmod(grid,Np,meshx,deltax,numregioni)
	implicit none
c
c	nmax e' la dimensione max del vettore x() che conterra' la griglia
c		
	integer*4 nmax,id,Npmax	
	parameter(nmax=20000,Npmax=10000)	 
c       
	integer*4 numregioni,sommapunti(3),i,initial,totpunti(3),count
	integer*4 numpassi(3),j,terminal(3),regione,punti(3),flg3,flagcnt
	integer*4 n,Nring,indmin,indmax,Nnewx,Nnewy,flagGNR,rank
	integer*4 Nnewz,nxn,nyn,nzn,Numaxis,Np
	real xiniziale,xfinale,deltainit,deltafin
	real fraz,x(Npmax),ratio(3),coeff,NUMSOLIDI
	double precision xxnew(Npmax),yynew(Npmax),zznew(Npmax)
	double precision    grid(Np)
	double precision meshx(nmax),meshy(nmax),meshz(nmax)
	double precision deltax(nmax),deltay(nmax),deltaz(nmax)
	double precision xins,xfins,yins,yfins,zins,zfins,x1
	double precision scambiox(Npmax),scambioy(Npmax),scambioz(Npmax)
          character*20   dirname_in
          character*20   dirname_out 
          character*7    preno
          character*9    init
          character*1    eq
          character*1    t1
          character*1    t2 
          character*1    vir
          character*2    n1
          character*2    n2
          character*2    n3
          character*28   tot
          character*4    esten
          character*1    lab(3)
          character*12   name_files(3)
          init='parameter'
          eq='='
          t1='('
          t2=')'
          vir=','
          n1='nx'
          n2='ny'      
          n3='nz'
          preno='griglia'
          esten='.out'
          lab(1)='x'
          lab(2)='y'
          lab(3)='z'

c
c	apertura file di ingresso
c	                       
c	
c       'griglia.dat' contiene il file di ingresso
c           che e' composto da:
c	   - parametri per la simulazione
c          - per ogni direzione da: 
c          - numero di regioni
c         
c           e per ogni regione da:            
c          -  indice di regione     
c          -  x iniziale       
c          -  x finale 
c          -  passo iniziale
c          -  passo finale
c	
	  
	Numaxis=1;
	do 100 id=1,Numaxis  
	sommapunti(id) = 1
	punti(id) = 1
	
c       lettura dei dati per ogni regione
c       
	do 40 i=1,numregioni
	   
	      initial = sommapunti(id)
c		      
	      xiniziale=meshx(i)

	      xfinale=meshx(i+1)
	      
	      deltainit=deltax(i)
	      deltafin=deltax(i+1)
	      

c	Se deltaint = deltafin la procedura e' diversa => goto 10
c	
	if (deltainit .EQ. deltafin) goto 10
c                                
c   	dati xiniziale, xfinale, il deltax all'inizio della zona e 
c       il deltax alla fine della zona, e' possibile determinare:
c	 -la ragione (ratio) con cui varia il deltax
c        -il numero di passi necessari per coprire la regione 
c                                 
        ratio(id)=(deltainit-xfinale+xiniziale)
     &   /(deltafin-xfinale+xiniziale)
	fraz = deltafin/deltainit
	numpassi(id) = nint((log(fraz)/log(ratio(id)))+1)
c       
c	N.B. il numero di passi non risulta un numero intero!
c	=> e' necessario arrotondarlo (con NINT) e poi ridefinire
c	il deltax iniziale in funzione del nuovo valore.
c	e' normale quindi che deltax iniziale e finale risultino
c	leggermente differenti (2-3%) dai valori impostati.
c                                 	
	deltainit = (xfinale-xiniziale)
     &  *((1-ratio(id))/(1-ratio(id)**numpassi(id)))
c	
	goto 20
c	
c	procedura nel caso deltainit = deltafin
c
 10	ratio(id) = 1
	numpassi(id) = nint((xfinale-xiniziale)/deltainit)
	deltainit = (xfinale-xiniziale)/numpassi(id)
c	
 20	continue
c	
c	Dopo aver calcolato i parametri necessari 
c	=> calcolo i singoli punti e li inserisco nel vettore
c
c	inserimento del punto iniziale di ogni zona 
c       nel vettore x() della griglia
c
	x(initial) = xiniziale
c
c	apertura file di uscita e                             
c	inserimento di x(iniziale) nel file di uscita
c
cc        encode(12,888,name_files(id))
cc     &  preno,lab(id),esten
	
cc	open(unit=10+id,file=name_files(id),status='unknown')

	grid(initial)=x(initial)

	punti(id) = punti(id)+1 
c	
c	terminal contiene il penultimo punto della zona in esame
c       N.B. l'ultimo punto di ogni zona concide col primo della zona seguente 
c       => non e' necessario calcolarlo!
c   
	terminal(id) = sommapunti(id)+numpassi(id)-1
c    
c	calcolo dei punti della griglia per la zona in esame
c	e loro inserimento nel file di uscita
c	N.B. dal secondo al penultimo!
c	 
	do 30 j=initial+1,terminal(id)
c		 
	   grid(j)=grid(j-1)+deltainit*ratio(id)**(j-initial-1)
	   punti(id) = punti(id)+1
30	continue
c	
c	aggiornamento del contatore per il vettore x()
c	  	
        sommapunti(id) = sommapunti(id)+numpassi(id)
40	continue
 
c	
c	chiusura del file di ingresso
c
c	
c	inserimento dell'ultimo punto della struttura nel vettore x()
c	e poi nel file di uscita
c	 
	grid(terminal(id)+1)=xfinale
c	
c	totpunti contiene il numero di punti che costituiscono la griglia
c		
	totpunti(id) = punti(id)
	if (id.eq.1) then
	   nxn=totpunti(id)
	elseif (id.eq.2) then
	   nyn=totpunti(id)
	elseif (id.eq.3) then
	   nzn=totpunti(id)
	endif
c	
c	chiusura file di uscita
c	
 100	continue	
	

	
 7777	continue

	Np=nxn
	
	return
	end
	
	
	
