      subroutine sub_jay1d(n1d,autovalori,gy,ny,mu,genric,
     $     tolljay3,vt,nlims,nlimd,Dn)
      implicit none
      integer*4 i,j,k,kk,nz,ny,chess(ny),ciclo
      double precision gy(ny),effo(ny),n1d(ny),autovalori(ny)
      double precision normaf,normaf2,a1,a2,a3,tolljay3,
     $     mu(ny),genric(ny),q,bern,vt,nlims,nlimd,Dn(ny),beta(ny)
      q=1.60219e-19
c     creo il vettore chess
      do i=1,ny
         chess(i)=mod(i,2)
         beta(i)=vt
c     Dn(i)/mu(i)
      enddo      
C     SOLUZIONE INIZIALE DEL POTENZIALE
      open(unit=33,file="outt",status="unknown")
      do j=1,ny
         write(33,*) gy(j),autovalori(j),mu(j),genric(j),n1d(j)
      enddo
      close(33)

      normaf=10
      ciclo=1
 103  k=1
 101  do i=1,ny
         if (chess(i).eq.k) then
            if (i.eq.1) then
               n1d(i)=nlims
            elseif (i.eq.ny) then
               n1d(i)=nlimd
            else
               a1=q*0.5*(mu(i)+mu(i-1))*beta(i)*
     $              bern(-(autovalori(i-1)-autovalori(i))/beta(i))/
     $              (gy(i)-gy(i-1))
               a2=q*0.5*(mu(i)+mu(i+1))*beta(i)*
     $              bern(-(autovalori(i+1)-autovalori(i))/beta(i))/
     $              (gy(i+1)-gy(i))
               a3=q*0.5*(mu(i)+mu(i-1))*beta(i)*
     $              bern((autovalori(i-1)-autovalori(i))/beta(i))/
     $              (gy(i)-gy(i-1))
     $              +q*0.5*(mu(i)+mu(i+1))*beta(i)*
     $              bern((autovalori(i+1)-autovalori(i))/beta(i))/
     $              (gy(i+1)-gy(i))
               
               n1d(i)=1/a3*(a1*n1d(i-1)+a2*n1d(i+1))
     $              -genric(i)*0.5*(gy(i+1)-gy(i-1))
            endif
         endif
      enddo
      if (k.eq.0) goto 102
      k=0
      goto 101
 102  continue
      if (mod(ciclo,1).eq.0) write(6,*) 'Iterazione N.',ciclo
      ciclo=ciclo+1
      call norma2jay3(n1d,normaf2,ny)
      if (abs(normaf-normaf2).ge.tolljay3) then
         if (mod(ciclo,100000).eq.0) 
     $        write(6,*) 'Differenza norma',abs(normaf-normaf2)
         normaf=normaf2
         goto 103
      endif
c$$$      open(unit=33,file='out',status='unknown')
c$$$      do i=1,ny
c$$$         write(33,*) gy(i),n1d(i),autovalori(i)
c$$$      enddo
c$$$      close(33)
c$$$      pause

c$$$      open(unit=33,file='out',status='unknown')
c$$$      do i=2,ny-1
c$$$         a2=q*0.5*(mu(i)+mu(i+1))*beta(i)*
c$$$     $        bern(-(autovalori(i+1)-autovalori(i))/beta(i))/
c$$$     $        (gy(i+1)-gy(i))*n1d(i+1)
c$$$         a3=q*0.5*(mu(i)+mu(i+1))*beta(i)*
c$$$     $        bern((autovalori(i+1)-autovalori(i))/beta(i))/
c$$$     $        (gy(i+1)-gy(i))*n1d(i)
c$$$         write(33,*) gy(i),a2-a3
c$$$      enddo
c$$$      close(33)
c$$$      pause
      end
      
      subroutine norma2jay3(vet,norma,Np)
      implicit none
      integer*4 i,Np
      double precision vet(Np),sommaq,norma
      sommaq=0
      do i=1,Np
         sommaq=sommaq+vet(i)**2
      end do
      norma=dsqrt(sommaq)
      end
      
      double precision function bern(x)
      implicit none
      double precision x,tmp,f,y,d
      d=8.333333333333333e-2
      if (abs(x).le.1e-6) then
         y=1.0+x*((x*d)-0.5)
      else
         tmp=exp(x)
         f=1.0/(tmp-1.0)
         y=x*f
      endif
      bern=y
      return
      end
