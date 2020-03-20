c ======================================================================
c  Copyright (c) 2012,  P. Marconcini, G. Fiori, G. Iannaccone  University of Pisa
c
c  This file is released under the BSD license.
c  See the file "license.txt" for information on usage and
c  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
c ====================================================================== 



      subroutine passivationsiavo(x,y,z,numfn,xfn,yfn,zfn,deltae,hpass)
c
c     x, y, z (input parameters) are the coordinates of the
c     considered silicon atom;
c     numf (input parameter) is the number of its first-neighbour
c     silicon atoms;
c     the first numf elements of the vectors xfn, yfn, zfn
c     (input parameters) are the coordinates of its numf
c     first-neighbour silicon atoms;
c     deltae (input parameter) is the energy by which we increase
c     the onsite energy of the hybridized orbitals pointing
c     outside the nanowire;
c     hpass (output parameter) is the matrix to sum to the
c     hamiltonian to emulate the passivation.
c
      implicit none
      integer i,j,numfn,ancat,l1,l2,l3,l4
c    #,upp1,upp2,doo1,doo2
      double precision x,y,z,xfn(4),yfn(4),zfn(4),dx,dy,dz,
     #deltae,hpass(4,4),const
c
      ancat=0
c     type of atom: anion (1) or cation (2)
      l1=1
      l2=1
      l3=1
      l4=1
c
      do i=1,numfn
c
         dx=xfn(i)-x
         dy=yfn(i)-y
         dz=zfn(i)-z
c
         if ((dx.gt.0).and.(dy.gt.0).and.(dz.gt.0)) then
c        silicon atom in the direction [1 1 1] with respect to an anion
            if (ancat.eq.2) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=1
               l1=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.lt.0).and.(dy.lt.0).and.(dz.gt.0)) then
c        silicon atom in the direction [-1 -1 1] with respect to an anion
            if (ancat.eq.2) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=1
               l2=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.gt.0).and.(dy.lt.0).and.(dz.lt.0)) then
c        silicon atom in the direction [1 -1 -1] with respect to an anion
            if (ancat.eq.2) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=1
               l3=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.lt.0).and.(dy.gt.0).and.(dz.lt.0)) then
c        silicon atom in the direction [-1 1 -1] with respect to an anion
            if (ancat.eq.2) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=1
               l4=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.lt.0).and.(dy.gt.0).and.(dz.gt.0)) then
c        silicon atom in the direction [-1 1 1] with respect to a cation
            if (ancat.eq.1) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=2
               l1=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.gt.0).and.(dy.lt.0).and.(dz.gt.0)) then
c        silicon atom in the direction [1 -1 1] with respect to a cation
            if (ancat.eq.1) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=2
               l2=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.lt.0).and.(dy.lt.0).and.(dz.lt.0)) then
c        silicon atom in the direction [-1 -1 -1] with respect to a cation
            if (ancat.eq.1) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=2
               l3=0
c              no passivation in that direction
            end if
         end if
c
         if ((dx.gt.0).and.(dy.gt.0).and.(dz.lt.0)) then
c        silicon atom in the direction [1 1 -1] with respect to a cation
            if (ancat.eq.1) then
               write(*,*) 'error in ancat!'
               return
            else
               ancat=2
               l4=0
c              no passivation in that direction
            end if
         end if
c
      end do
c
c      if (ancat.eq.1) then
c         if (l1.eq.1) upp1=-2
c         if (l2.eq.1) upp2=-2
c         if (l3.eq.1) doo1=-2
c         if (l4.eq.1) doo2=-2
c      else
c         if (l1.eq.1) upp1=-1
c         if (l2.eq.1) upp2=-1
c         if (l3.eq.1) doo1=-1
c         if (l4.eq.1) doo2=-1
c      end if
c
      if (ancat.eq.1) then
c
         hpass(1,1)=l1+l2+l3+l4
         hpass(1,2)=l1-l2+l3-l4
         hpass(1,3)=l1-l2-l3+l4
         hpass(1,4)=l1+l2-l3-l4
         hpass(2,1)=l1-l2+l3-l4
         hpass(2,2)=l1+l2+l3+l4
         hpass(2,3)=l1+l2-l3-l4
         hpass(2,4)=l1-l2-l3+l4
         hpass(3,1)=l1-l2-l3+l4
         hpass(3,2)=l1+l2-l3-l4
         hpass(3,3)=l1+l2+l3+l4
         hpass(3,4)=l1-l2+l3-l4
         hpass(4,1)=l1+l2-l3-l4
         hpass(4,2)=l1-l2-l3+l4
         hpass(4,3)=l1-l2+l3-l4
         hpass(4,4)=l1+l2+l3+l4
c
      else
c
         hpass(1,1)=l1+l2+l3+l4
         hpass(1,2)=-l1+l2-l3+l4
         hpass(1,3)=l1-l2-l3+l4
         hpass(1,4)=l1+l2-l3-l4
         hpass(2,1)=-l1+l2-l3+l4
         hpass(2,2)=l1+l2+l3+l4
         hpass(2,3)=-l1-l2+l3+l4
         hpass(2,4)=-l1+l2+l3-l4
         hpass(3,1)=l1-l2-l3+l4
         hpass(3,2)=-l1-l2+l3+l4
         hpass(3,3)=l1+l2+l3+l4
         hpass(3,4)=l1-l2+l3-l4
         hpass(4,1)=l1+l2-l3-l4
         hpass(4,2)=-l1+l2+l3-l4
         hpass(4,3)=l1-l2+l3-l4
         hpass(4,4)=l1+l2+l3+l4
c
      end if
c
      const=deltae/4.0d0
c
      do i=1,4
         do j=1,4
            hpass(i,j)=hpass(i,j)*const
         end do
      end do   
c
      end
