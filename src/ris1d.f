      subroutine ris1d(ecs,xs,ms,
     $     nxs,eigval,eigfun,nauto)
      implicit none
      integer nmaxs1
      integer*4 nxs,i,j,k,ix,l,i1,j1,ii
      integer ierr,nauto
      double precision ecs(nxs),xs(nxs),m0,
     $     ms(nxs),eigfun(nxs*nauto),zzz(nxs-2,nxs-2),
     $     a(nxs-2),b(nxs-2),hbar,q,deltai(nxs),eigval(nauto),
     $     di,dip,bi,xs2(nxs),ms2(nxs)
      hbar=1.05459e-34
      q=1.60219e-19
      m0=9.1095e-31
      
      do i=1,nxs*nauto
         eigfun(i)=0
      enddo

      do i=1,nxs
         xs2(i)=xs(i)*1e-9
         ms2(i)=ms(i)*m0
      enddo
      
      
c$$$      open(unit=67,file='BC2.out',status='unknown')
c$$$      do i=1,nxs
c$$$         write(67,*) xs2(i),ecs(i)
c$$$      enddo
c$$$      close(67)

c     inizializzo i deltai
        

      do i=1,nxs
         if (i.eq.1) then
            deltai(i)=sqrt(2/(xs2(i+1)-xs2(i)))
         elseif (i.eq.nxs) then
            deltai(i)=sqrt(2/(xs2(i)-xs2(i-1)))
         else
            deltai(i)=sqrt(2/(xs2(i+1)-xs2(i-1)))
         endif
      enddo

      
      do i=1,nxs
         if ((i.ne.1).and.(i.ne.nxs)) then
            a(i-1)=hbar**2/(2*q)*
     $           (1/((ms2(i))
     $           *(xs2(i+1)-xs2(i)))+
     $           1/((ms2(i-1))
     $           *(xs2(i)-xs2(i-1))))
     $           +0.5*(xs2(i+1)-xs2(i-1))*ecs(i)
            a(i-1)=a(i-1)*2/(xs2(i+1)-xs2(i-1))
            
            
            if (i.ne.nxs-1) then
               
               b(i)=-hbar**2/(2*q*(
     $              ms2(i))
     $              *(xs2(i+1)-xs2(i)))
               di=deltai(i)
               dip=deltai(i+1)
               bi=b(i)
               b(i)=deltai(i)*b(i)*deltai(i+1)
            endif
         endif
      enddo


c$$$      open(unit=67,file='a10.out',status='unknown')
c$$$      open(unit=68,file='b10.out',status='unknown')
c$$$      do i=1,nxs-2
c$$$         write(67,*) a(i)
c$$$         write(68,*) b(i)
c$$$      enddo
c$$$      close(67)
c$$$      close(68)

c     rendo identica la matrice iniziale
      
      do i1=1,nxs-2
         do j1=1,nxs-2
            if(i1.eq.j1) then
               zzz(i1,j1)=1
            else
               zzz(i1,j1)=0
            endif
         enddo
      enddo
      
      ierr=0
      b(1)=0
      call imtql2(nxs-2,nxs-2,a,b,zzz,ierr)
      if (ierr.ne.0) then
         print*,'IMTQL2 ha dato un errore',ierr
         pause
      endif


      do l=1,nauto
         do i=1,nxs
            if ((i.ne.1).and.(i.ne.nxs)) then
               eigfun(i+(l-1)*nxs)=zzz(i-1,l)
     $              *sqrt(2/(xs2(i+1)-xs2(i-1))) 
            else
               eigfun(i+(l-1)*nxs)=0
            endif
         enddo
         eigval(l)=a(l)
      enddo
      

c$$$      open(unit=67,file='eigve.out',status='unknown')
c$$$      do i=1,nxs
c$$$         write(67,*) xs2(i),eigfun(i+nxs)
c$$$      enddo
c$$$      close(67)

      return
      end
      


      subroutine imtql2(nm,n,d,e,z,ierr) 
c 
      integer*4 i,j,k,l,m,n,ii,nm,mml 
      integer ierr 
      double precision d(n),e(n),z(nm,n) 
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag 
c 
c     this subroutine is a translation of the algol procedure imtql2, 
c     num. math. 12, 377-383(1968) by martin and wilkinson, 
c     as modified in num. math. 15, 450(1970) by dubrulle. 
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971). 
c 
c     this subroutine finds the eigenvalues and eigenvectors 
c     of a symmetric tridiagonal matrix by the implicit ql method. 
c     the eigenvectors of a full symmetric matrix can also 
c     be found if  tred2  has been used to reduce this 
c     full matrix to tridiagonal form. 
c 
c     on input 
c 
c        nm must be set to the row dimension of two-dimensional 
c          array parameters as declared in the calling program 
c          dimension statement. 
c 
c        n is the order of the matrix. 
c 
c        d contains the diagonal elements of the input matrix. 
c 
c        e contains the subdiagonal elements of the input matrix 
c          in its last n-1 positions.  e(1) is arbitrary. 
c 
c        z contains the transformation matrix produced in the 
c          reduction by  tred2, if performed.  if the eigenvectors 
c          of the tridiagonal matrix are desired, z must contain 
c          the identity matrix. 
c 
c      on output 
c 
c        d contains the eigenvalues in ascending order.  if an 
c          error exit is made, the eigenvalues are correct but 
c          unordered for indices 1,2,...,ierr-1. 
c 
c        e has been destroyed. 
c 
c        z contains orthonormal eigenvectors of the symmetric 
c          tridiagonal (or full) matrix.  if an error exit is made, 
c          z contains the eigenvectors associated with the stored 
c          eigenvalues. 
c 
c        ierr is set to 
c          zero       for normal return, 
c          j          if the j-th eigenvalue has not been 
c                     determined after 30 iterations. 
c 
c     calls pythag for  dsqrt(a*a + b*b) . 
c 
c     questions and comments should be directed to burton s. garbow, 
c     mathematics and computer science div, argonne national laboratory 
c 
c     this version dated august 1983. 
c 
c     ------------------------------------------------------------------ 
c 
      ierr = 0 
      if (n .eq. 1) go to 1001 
c 
      do 100 i = 2, n 
  100 e(i-1) = e(i) 
c 
      e(n) = 0.0d0 
c 
      do 240 l = 1, n 
         j = 0 
c     .......... look for small sub-diagonal element .......... 
  105    do 110 m = l, n 
            if (m .eq. n) go to 120 
            tst1 = dabs(d(m)) + dabs(d(m+1)) 
            tst2 = tst1 + dabs(e(m)) 
            if (tst2 .eq. tst1) go to 120 
  110    continue 
c 
  120    p = d(l) 
         if (m .eq. l) go to 240 
         if (j .eq. 30) go to 1000 
         j = j + 1 
c     .......... form shift .......... 
         g = (d(l+1) - p) / (2.0d0 * e(l)) 
         r = pythag(g,1.0d0) 
         g = d(m) - p + e(l) / (g + dsign(r,g)) 
         s = 1.0d0 
         c = 1.0d0 
         p = 0.0d0 
         mml = m - l 
c     .......... for i=m-1 step -1 until l do -- .......... 
         do 200 ii = 1, mml 
            i = m - ii 
            f = s * e(i) 
            b = c * e(i) 
            r = pythag(f,g) 
            e(i+1) = r 
            if (r .eq. 0.0d0) go to 210 
            s = f / r 
            c = g / r 
            g = d(i+1) - p 
            r = (d(i) - g) * s + 2.0d0 * c * b 
            p = s * r 
            d(i+1) = g + p 
            g = c * r - b 
c     .......... form vector .......... 
            do 180 k = 1, n 
               f = z(k,i+1) 
               z(k,i+1) = s * z(k,i) + c * f 
               z(k,i) = c * z(k,i) - s * f 
  180       continue 
c 
  200    continue 
c 
         d(l) = d(l) - p 
         e(l) = g 
         e(m) = 0.0d0 
         go to 105 
c     .......... recover from underflow .......... 
  210    d(i+1) = d(i+1) - p 
         e(m) = 0.0d0 
         go to 105 
  240 continue 
c     .......... order eigenvalues and eigenvectors .......... 
      do 300 ii = 2, n 
         i = ii - 1 
         k = i 
         p = d(i) 
c 
         do 260 j = ii, n 
            if (d(j) .ge. p) go to 260 
            k = j 
            p = d(j) 
  260    continue 
c 
         if (k .eq. i) go to 300 
         d(k) = d(i) 
         d(i) = p 
c 
         do 280 j = 1, n 
            p = z(j,i) 
            z(j,i) = z(j,k) 
            z(j,k) = p 
  280    continue 
c 
  300 continue 
c 
      go to 1001 
c     .......... set error -- no convergence to an 
c                eigenvalue after 30 iterations .......... 
 1000 ierr = l 
 1001 return 
      end 
 
 
      double precision function pythag(a,b) 
      double precision a,b 
c 
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow 
c 
      double precision p,r,s,t,u 
      p = dmax1(dabs(a),dabs(b)) 
      if (p .eq. 0.0d0) go to 20 
      r = (dmin1(dabs(a),dabs(b))/p)**2 
   10 continue 
         t = 4.0d0 + r 
         if (t .eq. 4.0d0) go to 20 
         s = r/t 
         u = 1.0d0 + 2.0d0*s 
         p = u*p 
         r = (s/u)**2 * r 
      go to 10 
   20 pythag = p 
      return 
      end 
