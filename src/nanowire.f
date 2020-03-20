c ======================================================================
c  Copyright (c) 2012,  P. Marconcini, G. Fiori, G. Iannaccone  University of Pisa
c
c  This file is released under the BSD license.
c  See the file "license.txt" for information on usage and
c  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
c ====================================================================== 


      subroutine nanowire( a, sqci, tilt, edge, zmax, Nc, 
     $     ics, ipsilon, zeta, n_aux, Nc_aux  )

      
c     
c    
c    
c
c     .. Integer Arguments ..
      implicit none
      integer sqci, tilt, n, Nc, n_aux, Nc_aux, k_in
      double precision zeta( n_aux*Nc_aux ), ics( n_aux*Nc_aux ),
     $     ipsilon( n_aux*Nc_aux )
   
c     .. double Arguments ..
      double precision edge, zmax
c     
c     atomic coordinates of a silicon nanowire with [001] axis
c     
c     implicit none
      integer np,na
      parameter(np=330,na=10000)
c     maximum number of points along x, y or z; maximum number of atoms
      integer i,j,k,nx,ny,nz,i1,j1,k1,i2,j2,k2,
     #ic,jc,kc,ic1,jc1,kc1,inside,ato(np,np,np),ato1(np,np,np),
     #ia(na),ja(na),ka(na),typea(na),nat,nat1,nalayer(np),
     #nea1(na),nea2(na),nea3(na),nea4(na),
     #count,count1,count2,nl,natom,nuat
c     #,inside1
      double precision a,a4,pi,inf,xmax,ymax,xma,yma,zma
c     #,a2
      if (np**3.lt.na) then
         write(*,*) 'error: np**3<na !'
         goto 10
      end if
c      a=5.431d0
c     in Angstrom
c      a2=a/2.0d0
      a4=a/4.0d0
      pi=acos(-1.0d0)
      inf=0.02d0*a4
c     infinitesimal quantity
c
c     COMMENT THE OLD PARAMETERS
c     READED FROM FILE, NOW PASSED AS PARAMETERS

c     open(unit=10,file='input.dat',status='unknown')
c     read(10,*) sqci
c     square(0), circular(1) or triangular(2) section
c     read(10,*) tilt
c     if square, it indicates if not tilted(0) or with a 45 degree tilt(1);
c     if triangular, it indicates if the long edge is along an axis(0) or 
c     along a diagonal(1)
c     read(10,*) edge
c     in Angstrom: 
c     edge of the square section, diameter of the circular section or 
c     short edge of the triangular section
c     read(10,*) zmax
c     in Angstrom: greatest z (length of the nanowire)
c      close(10)
c
      if (sqci.eq.0) then
         if (tilt.eq.0) then 
            xmax=edge
c           greatest value of x
            ymax=edge
c           greatest value of y
         end if
         if (tilt.eq.1) then 
            xmax=2.0d0*(edge/sqrt(2.0d0))
            ymax=2.0d0*(edge/sqrt(2.0d0))
         end if
      end if
      if (sqci.eq.1) then
         xmax=edge
         ymax=edge
      end if
      if (sqci.eq.2) then
         if (tilt.eq.0) then 
            xmax=2.0d0*(edge/sqrt(2.0d0))
c           greatest value of x
            ymax=edge/sqrt(2.0d0)
c           greatest value of y
         end if
         if (tilt.eq.1) then 
            xmax=edge
            ymax=edge
         end if
      end if
c
      xma=xmax+2*a
      yma=ymax+2*a
      zma=zmax+2*a+inf
      nx=int(xma/a4)
      ny=int(yma/a4)
      nz=int(zma/a4)
c     note: on the plane z=0 there are immediately anions
      if ((nx+1.gt.np).or.(ny+1.gt.np).or.(nz+1.gt.np)) then
         write(*,*) 'error: nx+1 or ny+1 or nz+1 is >np !'
         goto 10
      end if
c
c     initialization
c
      do k=1,nz+1
         nalayer(k)=0
         do j=1,ny+1
            do i=1,nx+1
               ato(i,j,k)=0
               ato1(i,j,k)=0
            end do
         end do
      end do
c
c     positions of the anions and cations inside the nanowire:
c
      do k=1,nz+1
         k1=k-1
c        i.e. k1=z/a4
         if (mod(k1,2).eq.0) then
            k2=k1/2
c           i.e. k2=z/a2
            do j=1,ny+1
               j1=j-4-1
c              i.e. j1=y/a4
               if (mod(j1,2).eq.0) then
                  j2=j1/2
c                 i.e. j2=y/a2
                  do i=1,nx+1
                     i1=i-4-1
c                    i.e. i1=x/a4
                     if (mod(i1,2).eq.0) then
                        i2=i1/2
c                       i.e. i2=x/a2
                        if (mod(abs(i2+j2+k2),2).eq.0) then
                           if (inside(i1,j1,k1,sqci,tilt,edge,zma, a)
     #                        .eq.1) ato(i,j,k)=1
c                          anion
                           ic1=i1+1
                           jc1=j1+1
                           kc1=k1+1
                           ic=i+1
                           jc=j+1
                           kc=k+1
                           if (inside(ic1,jc1,kc1,sqci,tilt,edge,zma, a)
     #                        .eq.1) ato(ic,jc,kc)=2
c                          cation
                        end if      
                     end if
                  end do
               end if
            end do
         end if
      end do 
c
c     passivation atoms:
c
      do k=1,nz+1
         do j=1,ny+1
            do i=1,nx+1
               if (ato(i,j,k).eq.1) then
c                  if (inside1(i+1,j+1,k+1).eq.1) then
                  if ((i+1.le.nx+1).and.(j+1.le.ny+1)
     #                .and.(k+1.le.nz+1)) then
                     if (ato(i+1,j+1,k+1).eq.0) ato(i+1,j+1,k+1)=-2
c                    cation-like passivation atom
                  end if
c                  if (inside1(i-1,j-1,k+1).eq.1) then
                  if ((i-1.ge.1).and.(j-1.ge.1)
     #                .and.(k+1.le.nz+1)) then
                     if (ato(i-1,j-1,k+1).eq.0) ato(i-1,j-1,k+1)=-2
                  end if
c                  if (inside1(i+1,j-1,k-1).eq.1) then
                  if ((i+1.le.nx+1).and.(j-1.ge.1)
     #                .and.(k-1.ge.1)) then
                     if (ato(i+1,j-1,k-1).eq.0) ato(i+1,j-1,k-1)=-2
                  end if
c                  if (inside1(i-1,j+1,k-1).eq.1) then
                  if ((i-1.ge.1).and.(j+1.le.ny+1)
     #                .and.(k-1.ge.1)) then
                     if (ato(i-1,j+1,k-1).eq.0) ato(i-1,j+1,k-1)=-2
                  end if

               end if
               if (ato(i,j,k).eq.2) then
c                  if (inside1(i-1,j+1,k+1).eq.1) then
                  if ((i-1.ge.1).and.(j+1.le.ny+1)
     #                .and.(k+1.le.nz+1)) then
                     if (ato(i-1,j+1,k+1).eq.0) ato(i-1,j+1,k+1)=-1
c                    anion-like passivation atom
                  end if
c                  if (inside1(i+1,j-1,k+1).eq.1) then
                  if ((i+1.le.nx+1).and.(j-1.ge.1)
     #                .and.(k+1.le.nz+1)) then
                     if (ato(i+1,j-1,k+1).eq.0) ato(i+1,j-1,k+1)=-1
                  end if
c                  if (inside1(i-1,j-1,k-1).eq.1) then
                  if ((i-1.ge.1).and.(j-1.ge.1)
     #                .and.(k-1.ge.1)) then
                     if (ato(i-1,j-1,k-1).eq.0) ato(i-1,j-1,k-1)=-1
                  end if
c                  if (inside1(i+1,j+1,k-1).eq.1) then
                  if ((i+1.le.nx+1).and.(j+1.le.ny+1)
     #                .and.(k-1.ge.1)) then
                     if (ato(i+1,j+1,k-1).eq.0) ato(i+1,j+1,k-1)=-1
                  end if
               end if
            end do
         end do
      end do
c
      nuat=0
      do k=1,nz+1
         do j=1,ny+1
            do i=1,nx+1
               if (ato(i,j,k).gt.0) then
                  nuat=nuat+1
               end if
            end do
         end do
      end do
      if (nuat.gt.na) then
         write(*,*) 'error: nuat>na !'
         goto 10
      end if
c
c     data for each atom:
c
      nat=0
      nat1=0
      do k=1,nz+1
         do j=1,ny+1
            do i=1,nx+1
               if (ato(i,j,k).gt.0) then
                  nalayer(k)=nalayer(k)+1
c                 number of atoms in each atom layer
                  nat=nat+1
                  if ((k.gt.4).and.(k.lt.nz-2)) nat1=nat1+1
                  ato1(i,j,k)=nat
                  ia(nat)=i
c                 in such a way that xa=(ia-4-1)*a4
                  ja(nat)=j
c                 in such a way that ya=(ja-4-1)*a4
                  ka(nat)=k
c                 in such a way that za=(ka-1)*a4
                  typea(nat)=ato(i,j,k)
               end if
               if (ato(i,j,k).lt.0) then
                  ato1(i,j,k)=ato(i,j,k)
               end if
            end do
         end do
      end do
c
      open(unit=13,file='nalayer.dat',status='unknown')
      write(13,*) nat
      write(13,*) nz+1
      do i=1,nz+1
         write(13,*) nalayer(i)
      end do
      close(13)
c
      open(unit=18,file='nalayer1.dat',status='unknown')
      write(18,*) nat1
      write(18,*) nz+1-8
      do i=1+4,nz+1-4
         write(18,*) nalayer(i)
      end do
      close(18)

      Nc = nz+1-8
      
c
c     nearest neighbors:
c     
c     nn-type (nearest neighbour type):
c     with respect to an anion, 1,2,3, and 4 are respectively in the
c     directions (1,1,1), (-1,-1,1), (1,-1,-1) and (-1,1,-1);
c     with respect to a cation, 1,2,3, and 4 are respectively in the
c     directions (-1,1,1), (1,-1,1), (-1,-1,-1) and (1,1,-1)
c
      do i=1,nat
         if (typea(i).eq.1) then
c            if (inside1(ia(i)+1,ja(i)+1,ka(i)+1).eq.1) then
            if ((ia(i)+1.le.nx+1).and.(ja(i)+1.le.ny+1)
     #          .and.(ka(i)+1.le.nz+1)) then
               nea1(i)=ato1(ia(i)+1,ja(i)+1,ka(i)+1)
            else
               nea1(i)=0
            end if
c            if (inside1(ia(i)-1,ja(i)-1,ka(i)+1).eq.1) then
            if ((ia(i)-1.ge.1).and.(ja(i)-1.ge.1)
     #          .and.(ka(i)+1.le.nz+1)) then
               nea2(i)=ato1(ia(i)-1,ja(i)-1,ka(i)+1)
            else
               nea2(i)=0
            end if
c            if (inside1(ia(i)+1,ja(i)-1,ka(i)-1).eq.1) then
            if ((ia(i)+1.le.nx+1).and.(ja(i)-1.ge.1)
     #          .and.(ka(i)-1.ge.1)) then
               nea3(i)=ato1(ia(i)+1,ja(i)-1,ka(i)-1)
            else
               nea3(i)=0
            end if
c            if (inside1(ia(i)-1,ja(i)+1,ka(i)-1).eq.1) then
            if ((ia(i)-1.ge.1).and.(ja(i)+1.le.ny+1)
     #          .and.(ka(i)-1.ge.1)) then
               nea4(i)=ato1(ia(i)-1,ja(i)+1,ka(i)-1)
            else
               nea4(i)=0
            end if
         end if
         if (typea(i).eq.2) then
c            if (inside1(ia(i)-1,ja(i)+1,ka(i)+1).eq.1) then
            if ((ia(i)-1.ge.1).and.(ja(i)+1.le.ny+1)
     #          .and.(ka(i)+1.le.nz+1)) then
               nea1(i)=ato1(ia(i)-1,ja(i)+1,ka(i)+1)
            else
               nea1(i)=0
            end if
c            if (inside1(ia(i)+1,ja(i)-1,ka(i)+1).eq.1) then
            if ((ia(i)+1.le.nx+1).and.(ja(i)-1.ge.1)
     #          .and.(ka(i)+1.le.nz+1)) then
               nea2(i)=ato1(ia(i)+1,ja(i)-1,ka(i)+1)
            else
               nea2(i)=0
            end if
c            if (inside1(ia(i)-1,ja(i)-1,ka(i)-1).eq.1) then
            if ((ia(i)-1.ge.1).and.(ja(i)-1.ge.1)
     #          .and.(ka(i)-1.ge.1)) then
               nea3(i)=ato1(ia(i)-1,ja(i)-1,ka(i)-1)
            else
               nea3(i)=0
            end if
c            if (inside1(ia(i)+1,ja(i)+1,ka(i)-1).eq.1) then
            if ((ia(i)+1.le.nx+1).and.(ja(i)+1.le.ny+1)
     #          .and.(ka(i)-1.ge.1)) then
               nea4(i)=ato1(ia(i)+1,ja(i)+1,ka(i)-1)
            else
               nea4(i)=0
            end if
         end if
      end do
c
c     output as a function of the number of the layer and of the atom 
c
      open(unit=15,file='infoatoms.dat',status='unknown')
      open(unit=16,file='posatoms.dat',status='unknown')
      do i=1,nat
         count=0
         j=1
 20      count1=count
         if (j.ne.1) then
            count2=count1-nalayer(j-1)
         end if
         count=count1+nalayer(j)
         if (i.le.count) then
            nl=j
            natom=i-count1
c
            write(15,*) nl,natom
            write(15,*) typea(i)
c           atom specified by its layer and its number over that layer
c           and its at-type (atom type): Sia=1,Sic=2
            write(16,*) nl,natom,ia(i)-4-1,ja(i)-4-1,ka(i)-1
c           position (in a4 units) of the atom specified by its layer 
c           and its number over that layer; the position of each row
c           in this file represents the original number of this atom
c
            if (nea1(i).gt.0) then
               write(15,*) nea1(i)-count,typea(nea1(i))
c              neighbor in the superior layer and its at-type
            else
               write(15,*) nea1(i),nea1(i)
            end if
c
            if (nea2(i).gt.0) then
               write(15,*) nea2(i)-count,typea(nea2(i))
c              other neighbor in the superior layer and its at-type
            else
               write(15,*) nea2(i),nea2(i)
            end if
c
            if (nea3(i).gt.0) then
               write(15,*) nea3(i)-count2,typea(nea3(i))
c              neighbor in the inferior layer and its at-type
            else
               write(15,*) nea3(i),nea3(i)
            end if
c
            if (nea4(i).gt.0) then
               write(15,*) nea4(i)-count2,typea(nea4(i))
c              other neighbor in the inferior layer and its at-type
            else
               write(15,*) nea4(i),nea4(i)
            end if
c
         else
            j=j+1
            goto 20
         end if
      end do
      close(15)
      close(16)
c
c     output to the files nanowire.xyz and nanowire1.xyz
c
c      open(unit=11,file='nanowire.xyz',status='unknown')
c      write(11,*) nat
c      write(11,*) ' '
c      do i=1,nat
c         write(11,*) 'Si',(ia(i)-4-1)*a4,(ja(i)-4-1)*a4,(ka(i)-1)*a4
c        in Angstrom
c      end do
c      close(11)
c
c      open(unit=17,file='nanowire1.xyz',status='unknown')
c      write(17,*) nat1
c      write(17,*) ' '
      k_in = 0 
      do i=1,nat
         if ((ka(i).gt.4).and.(ka(i).lt.nz-2)) then
c            write(17,*) 'Si', (ia(i)-4-1)*a4,(ja(i)-4-1)*a4,(ka(i)-1)*a4
            k_in = k_in + 1
            ics(k_in) = (ia(i)-4-1)*a4
            ipsilon(k_in) = (ja(i)-4-1)*a4
            zeta(k_in) = (ka(i)-1)*a4
c           in Angstrom
         end if
      end do
c      close(17)
      
      return
c
 10   end
      
c      return
      
c      end

      integer function inside(i1,j1,k1,sqci,tilt,edge,zma, a)
      implicit none
      integer i1,j1,k1,sqci,tilt
      double precision edge,zma,a,a4,inf,hedge,sqedge
c
c      a=5.431d0
c     in Angstrom
      a4=a/4.0d0
      inf=0.02d0*a4
      hedge=edge/2.0d0
      sqedge=edge/sqrt(2.0d0)
c
      if (sqci.eq.0) then
         if (tilt.eq.0) then
            if ((abs(i1*a4-hedge).le.hedge+inf).and.
     #          (abs(j1*a4-hedge).le.hedge+inf).and.
     #          (k1*a4.le.zma)) then
               inside=1
            else
               inside=0
            end if
         end if
         if (tilt.eq.1) then
            if ((abs(j1*a4-sqedge).le.
     #           -abs(i1*a4-sqedge)+(sqedge+inf)).and.
     #          (k1*a4.le.zma)) then
               inside=1
            else
               inside=0
            end if
         end if
      end if
c
      if (sqci.eq.1) then
         if ((sqrt((i1*a4-hedge)**2+(j1*a4-hedge)**2).le.(hedge+inf))
     #      .and.(k1*a4.le.zma)) then
            inside=1
         else
            inside=0
         end if
      end if
c
      if (sqci.eq.2) then
         if (tilt.eq.0) then
            if ((j1*a4.ge.-inf).and.
     #          (j1*a4.le.-abs(i1*a4-sqedge)+(sqedge+inf)).and.
     #          (k1*a4.le.zma)) then
               inside=1
            else
               inside=0
            end if
         end if
         if (tilt.eq.1) then
            if ((i1*a4.ge.-inf).and.(j1*a4.ge.-inf).and.
     #          (j1*a4.le.-i1*a4+(edge+inf)).and.
     #          (k1*a4.le.zma)) then
               inside=1
            else
               inside=0
            end if
         end if
      end if
c
      end


c      integer function inside1(i,j,k,nx,ny,nz)
c      implicit none
c      integer i,j,k,nx,ny,nz
c      if ((i.ge.1).and.(i.le.nx+1).and.
c     #    (j.ge.1).and.(j.le.ny+1).and.
c     #    (k.ge.1).and.(k.le.nz+1)) then
c         inside1=1
c      else
c         inside1=0
c      end if
c      end
