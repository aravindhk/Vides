      SUBROUTINE hamiltonian( deltae, number, 
     $     H_aux, n_aux, Nc_aux, k_in, skparameters  )
c     
c     hamiltonian of a silicon nanowire
c     
      implicit none
      
c     passed parameters
      double precision deltae, skparameters(18)
      integer number, n_aux, Nc_aux
      double precision H_aux((Nc_aux*n_aux)*
     $     ((Nc_aux*n_aux+1)/2)*(2+100))

      integer nlay,natol, k_in
c     maximum number of layers, maximum number of atoms per layer 
      parameter(nlay=50,natol=20)
      integer i,j,nat,nla,nalayer(nlay),nl,natom,ty(nlay,natol),
     #up1(nlay,natol),up2(nlay,natol),tup1(nlay,natol),
     #tup2(nlay,natol),do1(nlay,natol),do2(nlay,natol),
     #tdo1(nlay,natol),tdo2(nlay,natol),
     #prev,prevorb(nlay,natol),norb(nlay),i1,j1,prorbl(nlay),nat1,
     #num1,num2,num1a,num1o,num2a,num2o,napre5,na(nlay,natol),pr,pr1,
     #prorbl1(nlay),nna,nna1,narid(nlay*natol),naext(nlay*natol),
     #nalamax,ordh,norbh(4),na1,na2,rem,nc,preorbla(4),pre,pre1,
     #ll,snumk,info,ordh1,ordhe,ordhe1,preorblae(4)
      double precision es,ep,ew,ed,ss1,ww1,sw1,sp1,wp1,sd1,wd1,
     #pp1,pp2,pd1,pd2,dd1,dd2,dd3,l,m,n,
     #q,pot(nlay,natol),diag(nlay,natol*10,natol*10),
c     #subdiag(nlay-1,natol*10,natol*10),
     #superdiag(nlay-1,natol*10,natol*10),val,
     #a,pi,dk,k,en(4*natol*10),val1,val2,rwork(3*4*natol*10)
c     #,hami1((nlay-6)*natol*10,(nlay-6)*natol*10),
c     #hamie1((nlay-6)*natol*10,(nlay-6)*natol*10)
      double complex ui,work(2*4*natol*10)
c     #,hami(4*natol*10,4*natol*10),hamie(4*natol*10,4*natol*10)
      character*14 str1*14,str2*14,str1a*30,str2a*30,
     #str1ae*31,str2ae*31
c
      q=1.60217733d-19
c
      open(unit=10,file='nalayer.dat',status='unknown')
      read(10,*) nat
c     total number of atoms
      read(10,*) nla
c     number of layers
      if (nla.gt.nlay) then
         write(*,*) 'error: nla>nlay !'
         goto 10
      end if
      do i=1,nla
         read(10,*) nalayer(i)
c        number of atoms in each layer
         if (nalayer(i).gt.natol) then
            write(*,*) 'error: nalayer(i)>natol !'
            goto 10
         end if
      end do
      close(10)
      napre5=nalayer(1)+nalayer(2)+nalayer(3)+nalayer(4)
c
      open(unit=11,file='infoatoms.dat',status='unknown')
      do i=1,nat
         read(11,*) nl,natom
         na(nl,natom)=i
         read(11,*) ty(nl,natom)
c        each atom of the nanowire
         read(11,*) up1(nl,natom),tup1(nl,natom)
         read(11,*) up2(nl,natom),tup2(nl,natom)
         read(11,*) do1(nl,natom),tdo1(nl,natom)
         read(11,*) do2(nl,natom),tdo2(nl,natom)
c        its 4 first neighbors
      end do
      close(11)
c
      do i=1,nla
         prev=0
         do j=1,nalayer(i)
            prevorb(i,j)=prev
c           number of atomic orbitals to be considered in the 
c           corresponding layer before the specified atom
            prev=prev+10
         end do
         norb(i)=prev
c        total number of atomic orbitals to be considered in each layer
      end do
c
c---------------------------------------------------
c
      nalamax=0
      do i=4,nla-3
         if (nalayer(i).gt.nalamax) then
            nalamax=nalayer(i)
         end if
      end do
c
      open(unit=80,file='nalamax.dat',status='unknown')
      write(80,*) nalamax
      close(80)

      number = nalamax
c
      open(unit=81,file='correspondence.dat',status='unknown')
      nna1=-4*nalamax
      nna=0
      do i=1,nla
         do j=1,nalamax
            nna1=nna1+1
            if (j.le.nalayer(i)) then
               nna=nna+1
               narid(nna)=nna-napre5
               naext(nna)=nna1
               write(81,*) narid(nna),naext(nna)
            end if
         end do
      end do
      close(81)
c
c---------------------------------------------------
c
c      open(unit=12,file='potential.dat',status='unknown')
c      do i=1,nla
c         do j=1,nalayer(i)
c            read(12,*) pot(i,j)
cc            pot(i,j)=pot(i,j)*q
ccc           in Joule
c         end do
c      end do
c      close(12)
c

      es = skparameters(1)
      ep= skparameters(2)
      ew= skparameters(3)
      ed= skparameters(4)
      ss1= skparameters(5)
      ww1= skparameters(6)
      sw1= skparameters(7)
      sp1= skparameters(8)
      wp1= skparameters(9)
      sd1= skparameters(10)
      wd1= skparameters(11)
      pp1= skparameters(12)
      pp2= skparameters(13)
      pd1= skparameters(14)
      pd2= skparameters(15)
      dd1= skparameters(16)
      dd2= skparameters(17)
      dd3= skparameters(18)

c      open(unit=13,file='skparameters.dat',status='unknown')
c      read(13,*) es
c      read(13,*) ep
c      read(13,*) ew
c      read(13,*) ed
c      read(13,*) ss1
c      read(13,*) ww1
c      read(13,*) sw1
c      read(13,*) sp1
c      read(13,*) wp1
c      read(13,*) sd1
c      read(13,*) wd1
c      read(13,*) pp1
c      read(13,*) pp2
c      read(13,*) pd1
c      read(13,*) pd2
c      read(13,*) dd1
c      read(13,*) dd2
c      read(13,*) dd3
c     w is s*, 1 is sigma, 2 is pi, 3 is delta
c      close(13)
c
c COMMENT TO THE OLD WAY OF PASSING DELTAE
c      open(unit=15,file='deltae.dat',status='unknown')
c      read(15,*) deltae
c      close(15)
c
c      es=es*q
c      ep=ep*q
c      ew=ew*q
c      ed=ed*q
c      ss1=ss1*q
c      ww1=ww1*q
c      sw1=sw1*q
c      sp1=sp1*q
c      wp1=wp1*q
c      sd1=sd1*q
c      wd1=wd1*q
c      pp1=pp1*q
c      pp2=pp2*q
c      pd1=pd1*q
c      pd2=pd2*q
c      dd1=dd1*q
c      dd2=dd2*q
c      dd3=dd3*q
c      deltae=deltae*q
cc     in Joule
c
      call matrix(*10,nla,nalayer,ty,up1,tup1,up2,tup2,
     #do1,tdo1,do2,tdo2,es,ep,ew,ed,ss1,ww1,sw1,sp1,
     #wp1,sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,deltae,
     #pot,prevorb,norb,diag,
c     #subdiag,
     #superdiag,nlay,natol)
cc
cc
c      open(unit=50,file='hamrid.dat',status='unknown')
cc
c      do i=5,nla-4
c         do j=1,nalayer(i)
c            pr=prevorb(i,j)
cc
c            do i1=1,10
c               do j1=1,10
c                  write(50,90) narid(na(i,j)),'(',i1,')',
c     #                         narid(na(i,j)),'(',j1,')',
c     #                                diag(i,pr+i1,pr+j1)
c               end do
c            end do
cc
c            if (i.ne.(nla-4)) then
cc
c               if (up1(i,j).lt.up2(i,j)) then
cc
c                  if (up1(i,j).gt.0) then
c                     pr1=prevorb(i+1,up1(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(50,90) narid(na(i,j)),'(',i1,')',
c     #                         narid(na(i+1,up1(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c                  if (up2(i,j).gt.0) then
c                     pr1=prevorb(i+1,up2(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(50,90) narid(na(i,j)),'(',i1,')',
c     #                         narid(na(i+1,up2(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c               else
cc
c                  if (up2(i,j).gt.0) then
c                     pr1=prevorb(i+1,up2(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(50,90) narid(na(i,j)),'(',i1,')',
c     #                         narid(na(i+1,up2(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c                  if (up1(i,j).gt.0) then
c                     pr1=prevorb(i+1,up1(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(50,90) narid(na(i,j)),'(',i1,')',
c     #                         narid(na(i+1,up1(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c               end if
cc
c            end if
cc
c         end do
c      end do
cc
c      close(50)
cc
cc
      open(unit=51,file='hamrid1.dat',status='unknown')
c
      write(51,91) '['
c
      write(51,92) '[',10,', ',0
      do i1=1,10
         do j1=1,10
            write(51,93) ', ',0.0d0
         end do
      end do
      write(51,91) ']'
c
      do i=5,nla-4
         do j=1,nalayer(i)
            pr=prevorb(i,j)
c
            write(51,92) '[',narid(na(i,j)),', ',
     #                       narid(na(i,j))
            do i1=1,10
               do j1=1,10
                  write(51,93) ', ',diag(i,pr+i1,pr+j1)
               end do
            end do
            write(51,91) ']'
c
            if (i.ne.(nla-4)) then
c
               if (up1(i,j).lt.up2(i,j)) then
c
                  if (up1(i,j).gt.0) then
                     pr1=prevorb(i+1,up1(i,j))
                     write(51,92) '[',narid(na(i,j)),', ',
     #                                narid(na(i+1,up1(i,j)))
                     do i1=1,10
                        do j1=1,10
                           write(51,93) ', ',superdiag(i,pr+i1,pr1+j1)
                        end do
                     end do
                     write(51,91) ']'
                  end if
c
                  if (up2(i,j).gt.0) then
                     pr1=prevorb(i+1,up2(i,j))
                     write(51,92) '[',narid(na(i,j)),', ',
     #                                narid(na(i+1,up2(i,j)))
                     do i1=1,10
                        do j1=1,10
                           write(51,93) ', ',superdiag(i,pr+i1,pr1+j1)
                        end do
                     end do
                     write(51,91) ']'
                  end if
c
               else
c
                  if (up2(i,j).gt.0) then
                     pr1=prevorb(i+1,up2(i,j))
                     write(51,92) '[',narid(na(i,j)),', ',
     #                                narid(na(i+1,up2(i,j)))
                     do i1=1,10
                        do j1=1,10
                           write(51,93) ', ',superdiag(i,pr+i1,pr1+j1)
                        end do
                     end do
                     write(51,91) ']'
                  end if
c
                  if (up1(i,j).gt.0) then
                     pr1=prevorb(i+1,up1(i,j))
                     write(51,92) '[',narid(na(i,j)),', ',
     #                                narid(na(i+1,up1(i,j)))
                     do i1=1,10
                        do j1=1,10
                           write(51,93) ', ',superdiag(i,pr+i1,pr1+j1)
                        end do
                     end do
                     write(51,91) ']'
                  end if
c
               end if
c
            end if
c
         end do
      end do
c
      write(51,91) ']'
c
      close(51)
c
c
      open(unit=40,file='norb.dat',status='unknown')
      do i=1,nla
         write(40,*) norb(i)
c        total number of atomic orbitals to be considered in each layer
      end do
      close(40)
c
      open(unit=44,file='norb1.dat',status='unknown')
      do i=4,nla-3
         write(44,*) norb(i)
c        total number of atomic orbitals to be considered in each layer
      end do
      close(44)
c
c      open(unit=45,file='diag.dat',status='unknown')
c      do i=1,nla
c         do i1=1,norb(i)
c            do j1=1,norb(i)
c               write(45,*) diag(i,i1,j1)
c            end do
c         end do
c      end do
c      close(45)
c
      open(unit=41,file='diag1.dat',status='unknown')
      do i=4,nla-3
         do i1=1,norb(i)
            do j1=1,norb(i)
               write(41,*) diag(i,i1,j1)
            end do
         end do
      end do
      close(41)
c
c      open(unit=46,file='superdiag.dat',status='unknown')
c      do i=1,nla-1
c         do i1=1,norb(i)
c            do j1=1,norb(i+1)
c               write(46,*) superdiag(i,i1,j1)
c            end do
c         end do
c      end do
c      close(46)
c
      open(unit=42,file='superdiag1.dat',status='unknown')
      do i=4,nla-4
         do i1=1,norb(i)
            do j1=1,norb(i+1)
               write(42,*) superdiag(i,i1,j1)
            end do
         end do
      end do
      close(42)
cc
cc      open(unit=47,file='subdiag.dat',status='unknown')
cc      do i=1,nla-1
cc         do i1=1,norb(i+1)
cc            do j1=1,norb(i)
cc               write(47,*) subdiag(i,i1,j1)
cc            end do
cc         end do
cc      end do
cc      close(47)
cc
c      open(unit=43,file='subdiag1.dat',status='unknown')
c      do i=4,nla-4
c         do i1=1,norb(i+1)
c            do j1=1,norb(i)
c               write(43,*) subdiag(i,i1,j1)
c            end do
c         end do
c      end do
c      close(43)
cc
c      ordh1=0
c      do i=4,nla-3
c         ordh1=ordh1+norb(i)
c      end do
c      open(unit=94,file='ordh1.dat',status='unknown')
c      write(94,*) ordh1
c      close(94)
cc
c      do i=1,ordh1
c         do j=1,ordh1
c            hami1(i,j)=0.0d0
c         end do
c      end do
cc
c      prorbl(4)=0
c      do i=4,nla-4
c         prorbl(i+1)=prorbl(i)+norb(i)
cc        number of atomic orbitals to be considered before the i-th layer
c      end do
cc
c      open(unit=48,file='hamiltonian.dat',status='unknown')
c      do i=4,nla-3
c         do i1=1,norb(i)
c            do j=4,nla-3
c               do j1=1,norb(j)
c                  if (i.eq.j+1) then
c                     write(48,*) prorbl(i)+i1,prorbl(j)+j1,
c     #                           subdiag(j,i1,j1)
c                     hami1(prorbl(i)+i1,prorbl(j)+j1)=subdiag(j,i1,j1)
c                  else if (i.eq.j) then
c                     write(48,*) prorbl(i)+i1,prorbl(j)+j1,
c     #                           diag(i,i1,j1)
c                     hami1(prorbl(i)+i1,prorbl(j)+j1)=diag(i,i1,j1)
c                  else if (i.eq.j-1) then
c                     write(48,*) prorbl(i)+i1,prorbl(j)+j1,
c     #                           superdiag(i,i1,j1)
c                     hami1(prorbl(i)+i1,prorbl(j)+j1)=superdiag(i,i1,j1)
c                  else
c                     write(48,*) prorbl(i)+i1,prorbl(j)+j1,
c     #                           0.0d0
c                     hami1(prorbl(i)+i1,prorbl(j)+j1)=0.0d0
c                  end if
c               end do
c            end do
c         end do
c      end do
c      close(48)
cc
cc     check if H is symmetrical
cc
c      open(unit=88,file='test1.dat',status='unknown')
c      do i=1,ordh1
c         do j=1,ordh1
c            if (hami1(i,j).ne.hami1(j,i)) then
c               write(88,*) i,j
c               write(88,*) 'i,j',hami1(i,j)
c               write(88,*) 'j,i',hami1(j,i)
c               write(88,*) 'diff',hami1(i,j)-hami1(j,i)
c               write(88,*) ''
c            end if
c         end do
c      end do
c      close(88)
cc
c      nat1=0
c      do i=4,nla-3
c         nat1=nat1+nalayer(i)
c      end do
cc
c      open(unit=48,file='hamiltonian.dat',status='unknown')
c      open(unit=49,file='hamiltonianbis.dat',status='unknown')
c      do i=1,nat1*10
c         do j=1,nat1*10
c            read(48,*) num1,num2,val
c            num1a=1+(num1-1)/10-nalayer(4)
c            num1o=num1-(num1a-1)*10
c            if (num1o.eq.1) then
c               str1='s             '
c            else if (num1o.eq.2) then
c               str1='px            '
c            else if (num1o.eq.3) then
c               str1='py            '
c            else if (num1o.eq.4) then
c               str1='pz            '
c            else if (num1o.eq.5) then
c               str1='dxy           '
c            else if (num1o.eq.6) then
c               str1='dyz           '
c            else if (num1o.eq.7) then
c               str1='dzx           '
c            else if (num1o.eq.8) then
c               str1='d(x**2-y**2)  '
c            else if (num1o.eq.9) then
c               str1='d(3z**2-r**2) '
c            else if (num1o.eq.10) then
c               str1='s*            '
c            end if
c            num2a=1+(num2-1)/10-nalayer(4)
c            num2o=num2-(num2a-1)*10
c            if (num2o.eq.1) then
c               str2='s             '
c            else if (num2o.eq.2) then
c               str2='px            '
c            else if (num2o.eq.3) then
c               str2='py            '
c            else if (num2o.eq.4) then
c               str2='pz            '
c            else if (num2o.eq.5) then
c               str2='dxy           '
c            else if (num2o.eq.6) then
c               str2='dyz           '
c            else if (num2o.eq.7) then
c               str2='dzx           '
c            else if (num2o.eq.8) then
c               str2='d(x**2-y**2)  '
c            else if (num2o.eq.9) then
c               str2='d(3z**2-r**2) '
c            else if (num2o.eq.10) then
c               str2='s*            '
c            end if
c            write(49,*) num1a,str1,num2a,str2,val
c         end do
c      end do
c      close(48)
c      close(49)
cc
cc/////////////////////////////////////////////////
cc
c      ordh=0
c      do i=5,8
c         ordh=ordh+norb(i)
c      end do
c      open(unit=83,file='ordh.dat',status='unknown')
c      write(83,*) ordh
c      close(83)
cc
c      open(unit=87,file='snumk.dat',status='unknown')
c      read(87,*) snumk
cc     half-number of k's to consider in the nanowire Brillouin zone
c      close(87)      
cc
c      do i=1,4
c         norbh(i)=norb(i+4)
c      end do
c      preorbla(1)=0
c      do i=1,3
c         preorbla(i+1)=preorbla(i)+norbh(i)
c      end do
cc
c      a=5.431d0
cc     in Angstrom
c      pi=acos(-1.0d0)
c      dk=pi/(a*snumk)
cc     in Angstrom**(-1)
c      ui=(0.0d0,1.0d0)
cc     imaginary unit
cc
c      do ll=1,snumk+1
c         k=(ll-1)*dk
cc
c         do i=1,ordh
c            do j=1,ordh
c               hami(i,j)=0.0d0
c            end do
c         end do
cc
c         do i=1,4
c            pre=preorbla(i)
c            do i1=1,norbh(i)
c               do j1=1,norbh(i)
c                  hami(pre+i1,pre+j1)=diag(i+4,i1,j1)
c               end do
c            end do
c         end do
c         do i=1,3
c            pre=preorbla(i)
c            pre1=preorbla(i+1)
c            do i1=1,norbh(i)
c               do j1=1,norbh(i+1)
c                  hami(pre+i1,pre1+j1)=superdiag(i+4,i1,j1)
c                  hami(pre1+j1,pre+i1)=subdiag(i+4,j1,i1)
c               end do
c            end do
c         end do
c         pre=preorbla(4)
c         pre1=preorbla(1)
c         do i1=1,norbh(4)
c            do j1=1,norbh(1)
c               hami(pre+i1,pre1+j1)=superdiag(8,i1,j1)*exp(ui*k*a)
c               hami(pre1+j1,pre+i1)=subdiag(4,j1,i1)*exp(-ui*k*a)
c            end do
c         end do
cc
cc        check if H is Hermitian
cc
c         open(unit=84,file='test2.dat',status='unknown',access='append')
c         do i=1,ordh
c            do j=1,ordh
c               if (hami(i,j).ne.dconjg(hami(j,i))) then
c                  write(84,*) ll,i,j
c                  write(84,*) 'i,j',dreal(hami(i,j)),dimag(hami(i,j))
c                  write(84,*) 'j,i',dreal(hami(j,i)),dimag(hami(j,i))
c                  write(84,*) 'diff',dreal(hami(i,j)-dconjg(hami(j,i))),
c     #                               dimag(hami(i,j)-dconjg(hami(j,i)))
c                  write(84,*) ''
c               end if
c            end do
c         end do
c         close(84)
cc
c         call zheev('N','U',ordh,hami,4*natol*10,en,work,2*ordh,
c     #               rwork,info)
cc        pay attention to put the correct leading dimension
cc        (i.e. the dimension of the variable declaration)!
cc
c         do i=1,ordh
c            if (i.lt.10) then
c               str1a=char(ichar('0')+i)
c               nc=1
c            else
c               na1=i
c               na2=na1/10
c               rem=mod(na1,10)
c               str1a=char(ichar('0')+rem)
c               nc=1
c               na1=na2
c 24            if (na2.lt.10) goto 34
c               na2=na1/10
c               rem=mod(na1,10)
c               str1a=char(ichar('0')+rem)//str1a(1:nc)
c               nc=nc+1
c               na1=na2
c               goto 24
c 34            str1a=char(ichar('0')+na1)//str1a(1:nc)
c               nc=nc+1
c            end if 
c            if (30.lt.(nc+8)) then
c                write(*,*) 'change 30: 30<(nc+8) !!'
c                goto 10
c            end if
c            str2a='band'//str1a(1:nc)//'.dat'
c            open(unit=84,file=str2a,status='unknown',access='append')
cc            write(84,*) k/(pi/a),en(i)/q
c            write(84,*) k/(pi/a),en(i)
cc           k in units of pi/a; energy in eV
c            close(84)
c         end do
c         if (ll.eq.1) then
c            open(unit=85,file='enk0.dat',status='unknown')
c            do i=1,ordh
c               write(85,*) en(i)
c            end do
c            close(85)
c         end if
cc
cc
c      end do
cc
cc
c      open(unit=86,file='bands.dat',status='unknown')
c      do i=1,ordh
c         if (i.lt.10) then
c            str1a=char(ichar('0')+i)
c            nc=1
c         else
c            na1=i
c            na2=na1/10
c            rem=mod(na1,10)
c            str1a=char(ichar('0')+rem)
c            nc=1
c            na1=na2
c 44         if (na2.lt.10) goto 54
c            na2=na1/10
c            rem=mod(na1,10)
c            str1a=char(ichar('0')+rem)//str1a(1:nc)
c            nc=nc+1
c            na1=na2
c            goto 44
c 54         str1a=char(ichar('0')+na1)//str1a(1:nc)
c            nc=nc+1
c         end if 
c         if (30.lt.(nc+8)) then
c            write(*,*) 'change 30: 30<(nc+8) !!'
c            goto 10
c         end if
c         str2a='band'//str1a(1:nc)//'.dat'
c         open(unit=85,file=str2a,status='unknown')
c         do ll=1,snumk+1
c            read(85,*) val1,val2
c            write(86,*) val1,val2
c         end do
c         close(85)
c         write(86,*) ' '
c      end do
c      close(86)
cc
cc/////////////////////////////////////////////////
cc
c      open(unit=70,file='exthamrid.dat',status='unknown')
cc
c      do i=5,nla-4
c         do j=1,nalayer(i)
c            pr=prevorb(i,j)
cc
c            do i1=1,10
c               do j1=1,10
c                  write(70,90) naext(na(i,j)),'(',i1,')',
c     #                         naext(na(i,j)),'(',j1,')',
c     #                                diag(i,pr+i1,pr+j1)
c               end do
c            end do
cc
c            if (i.ne.(nla-4)) then
cc
c               if (up1(i,j).lt.up2(i,j)) then
cc
c                  if (up1(i,j).gt.0) then
c                     pr1=prevorb(i+1,up1(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(70,90) naext(na(i,j)),'(',i1,')',
c     #                         naext(na(i+1,up1(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c                  if (up2(i,j).gt.0) then
c                     pr1=prevorb(i+1,up2(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(70,90) naext(na(i,j)),'(',i1,')',
c     #                         naext(na(i+1,up2(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c               else
cc
c                  if (up2(i,j).gt.0) then
c                     pr1=prevorb(i+1,up2(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(70,90) naext(na(i,j)),'(',i1,')',
c     #                         naext(na(i+1,up2(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c                  if (up1(i,j).gt.0) then
c                     pr1=prevorb(i+1,up1(i,j))
c                     do i1=1,10
c                        do j1=1,10
c                           write(70,90) naext(na(i,j)),'(',i1,')',
c     #                         naext(na(i+1,up1(i,j))),'(',j1,')',
c     #                                   superdiag(i,pr+i1,pr1+j1)
c                        end do
c                     end do
c                  end if
cc
c               end if
cc
c            end if
cc
c         end do
c      end do
cc
c      close(70)
cc
cc
      open(unit=71,file='exthamrid1.dat',status='unknown')
c
c     write(71,91) '['
c
      
      k_in = 0 
      write(71,95) 10,', ',0
      k_in = k_in + 1
      H_aux(k_in) = 10
      k_in = k_in + 1
      H_aux(k_in) = 0

      do i1=1,10
         do j1=1,10
            write(71,93) ', ',0.0d0
            k_in = k_in + 1
            H_aux(k_in) = 0
         end do
      end do
      write(71,91) ']'
c
      do i=5,nla-4
         do j=1,nalayer(i)
            pr=prevorb(i,j)
c
            write(71,92) '[',naext(na(i,j)),', ',
     #                       naext(na(i,j))

            k_in = k_in + 1
            H_aux(k_in) = naext(na(i,j))
c            write(*,*) naext(na(i,j))
c            write(*,*) H_aux(k_in)
            k_in = k_in + 1
            H_aux(k_in) = naext(na(i,j))

            do i1=1,10
               do j1=1,10
                  write(71,93) ', ',diag(i,pr+i1,pr+j1)
                  
                  k_in = k_in + 1
                  H_aux(k_in) = diag(i,pr+i1,pr+j1)

               end do
            end do
            if (i.eq.nla-4.and.j.eq.nalayer(nla-4)) then
               goto 100
            end if
            write(71,91) ']'
 100        continue

c
            if (i.ne.(nla-4)) then
c
               if (up1(i,j).lt.up2(i,j)) then
c
                  if (up1(i,j).gt.0) then
                     pr1=prevorb(i+1,up1(i,j))
                     write(71,92) '[',naext(na(i,j)),', ',
     #                                naext(na(i+1,up1(i,j)))

                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i,j))
                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i+1,up1(i,j)))

                     do i1=1,10
                        do j1=1,10
                           write(71,93) ', ',superdiag(i,pr+i1,pr1+j1)

                           k_in = k_in + 1
                           H_aux(k_in) = superdiag(i,pr+i1,pr1+j1)

                        end do
                     end do
                     write(71,91) ']'
                  end if
c
                  if (up2(i,j).gt.0) then
                     pr1=prevorb(i+1,up2(i,j))
                     write(71,92) '[',naext(na(i,j)),', ',
     #                                naext(na(i+1,up2(i,j)))

                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i,j))
                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i+1,up2(i,j)))

                     do i1=1,10
                        do j1=1,10
                           write(71,93) ', ',superdiag(i,pr+i1,pr1+j1)

                           k_in = k_in + 1
                           H_aux(k_in) = superdiag(i,pr+i1,pr1+j1)

                        end do
                     end do
                     write(71,91) ']'
                  end if
c
               else
c
                  if (up2(i,j).gt.0) then
                     pr1=prevorb(i+1,up2(i,j))
                     write(71,92) '[',naext(na(i,j)),', ',
     #                                naext(na(i+1,up2(i,j)))

                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i,j))
                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i+1,up2(i,j)))

                     do i1=1,10
                        do j1=1,10
                           write(71,93) ', ',superdiag(i,pr+i1,pr1+j1)

                           k_in = k_in + 1
                           H_aux(k_in) = superdiag(i,pr+i1,pr1+j1)

                        end do
                     end do
                     write(71,91) ']'
                  end if
c
                  if (up1(i,j).gt.0) then
                     pr1=prevorb(i+1,up1(i,j))
                     write(71,92) '[',naext(na(i,j)),', ',
     #                                naext(na(i+1,up1(i,j)))

                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i,j))
                     k_in = k_in + 1
                     H_aux(k_in) = naext(na(i+1,up1(i,j)))

                     do i1=1,10
                        do j1=1,10
                           write(71,93) ', ',superdiag(i,pr+i1,pr1+j1)

                           k_in = k_in + 1
                           H_aux(k_in) = superdiag(i,pr+i1,pr1+j1)

                        end do
                     end do
                     write(71,91) ']'
                  end if
c
               end if
c
            end if
c
         end do
      end do
c
c      write(71,91) ']'
c
      close(71)
c
c
c      open(unit=65,file='extdiag.dat',status='unknown')
c      do i=1,nla
c         do i1=1,nalamax*10
c            do j1=1,nalamax*10
c               if ((i1.le.norb(i)).and.(j1.le.norb(i))) then
c                  write(65,*) diag(i,i1,j1)
c               else
c                  write(65,*) 0.0d0
c               end if
c            end do
c         end do
c      end do
c      close(65)
c
      open(unit=61,file='extdiag1.dat',status='unknown')
      do i=4,nla-3
         do i1=1,nalamax*10
            do j1=1,nalamax*10
               if ((i1.le.norb(i)).and.(j1.le.norb(i))) then
                  write(61,*) diag(i,i1,j1)
               else
                  write(61,*) 0.0d0
               end if
            end do
         end do
      end do
      close(61)
c
c      open(unit=66,file='extsuperdiag.dat',status='unknown')
c      do i=1,nla-1
c         do i1=1,nalamax*10
c            do j1=1,nalamax*10
c               if ((i1.le.norb(i)).and.(j1.le.norb(i+1))) then
c                  write(66,*) superdiag(i,i1,j1)
c               else
c                  write(66,*) 0.0d0
c               end if
c            end do
c         end do
c      end do
c      close(66)
c
      open(unit=62,file='extsuperdiag1.dat',status='unknown')
      do i=4,nla-4
         do i1=1,nalamax*10
            do j1=1,nalamax*10
               if ((i1.le.norb(i)).and.(j1.le.norb(i+1))) then
                  write(62,*) superdiag(i,i1,j1)
               else
                  write(62,*) 0.0d0
               end if
            end do
         end do
      end do
      close(62)
cc
cc      open(unit=67,file='extsubdiag.dat',status='unknown')
cc      do i=1,nla-1
cc         do i1=1,nalamax*10
cc            do j1=1,nalamax*10
cc               if ((i1.le.norb(i+1)).and.(j1.le.norb(i))) then
cc                  write(67,*) subdiag(i,i1,j1)
cc               else
cc                  write(67,*) 0.0d0
cc               end if
cc            end do
cc         end do
cc      end do
cc      close(67)
cc
c      open(unit=63,file='extsubdiag1.dat',status='unknown')
c      do i=4,nla-4
c         do i1=1,nalamax*10
c            do j1=1,nalamax*10
c               if ((i1.le.norb(i+1)).and.(j1.le.norb(i))) then
c                  write(63,*) subdiag(i,i1,j1)
c               else
c                  write(63,*) 0.0d0
c               end if
c            end do
c         end do
c      end do
c      close(63)
cc
c      ordhe1=6*nalamax*10
c      open(unit=95,file='ordhe1.dat',status='unknown')
c      write(95,*) ordhe1
c      close(95)
cc
c      do i=1,ordhe1
c         do j=1,ordhe1
c            hamie1(i,j)=0.0d0
c         end do
c      end do
cc
c      do i=4,nla-3
c         prorbl1(i)=(i-4)*nalamax*10
c      end do
cc
c      open(unit=68,file='exthamiltonian.dat',status='unknown')
c      do i=4,nla-3
c         do i1=1,nalamax*10
c            do j=4,nla-3
c               do j1=1,nalamax*10
c                  if ((i1.le.norb(i)).and.(j1.le.norb(j))) then
c                     if (i.eq.j+1) then
c                        write(68,*) prorbl1(i)+i1,prorbl1(j)+j1,
c     #                              subdiag(j,i1,j1)
c                        hamie1(prorbl1(i)+i1,prorbl1(j)+j1)=
c     #                       subdiag(j,i1,j1)
c                     else if (i.eq.j) then
c                        write(68,*) prorbl1(i)+i1,prorbl1(j)+j1,
c     #                              diag(i,i1,j1)
c                        hamie1(prorbl1(i)+i1,prorbl1(j)+j1)=
c     #                       diag(i,i1,j1)
c                     else if (i.eq.j-1) then
c                        write(68,*) prorbl1(i)+i1,prorbl1(j)+j1,
c     #                              superdiag(i,i1,j1)
c                        hamie1(prorbl1(i)+i1,prorbl1(j)+j1)=
c     #                       superdiag(i,i1,j1)   
c                     else
c                        write(68,*) prorbl1(i)+i1,prorbl1(j)+j1,
c     #                              0.0d0
c                        hamie1(prorbl1(i)+i1,prorbl1(j)+j1)=
c     #                       0.0d0
c                     end if
c                  else
c                     write(68,*) prorbl1(i)+i1,prorbl1(j)+j1,0.0d0
c                     hamie1(prorbl1(i)+i1,prorbl1(j)+j1)=0.0d0
c                  end if
c               end do
c            end do
c         end do
c      end do
c      close(68)
cc
cc
cc     check if H is symmetrical
cc
c      open(unit=89,file='test3.dat',status='unknown')
c      do i=1,ordhe1
c         do j=1,ordhe1
c            if (hamie1(i,j).ne.hamie1(j,i)) then
c               write(89,*) i,j
c               write(89,*) 'i,j',hamie1(i,j)
c               write(89,*) 'j,i',hamie1(j,i)
c               write(89,*) 'diff',hamie1(i,j)-hamie1(j,i)
c               write(89,*) ''
c            end if
c         end do
c      end do
c      close(89)
c
c      open(unit=68,file='exthamiltonian.dat',status='unknown')
c      open(unit=69,file='exthamiltonianbis.dat',status='unknown')
c      do i=1,(nla-6)*nalamax*10
c         do j=1,(nla-6)*nalamax*10
c            read(68,*) num1,num2,val
c            num1a=1+(num1-1)/10-nalamax
c            num1o=num1-(num1a-1)*10
c            if (num1o.eq.1) then
c               str1='s             '
c            else if (num1o.eq.2) then
c               str1='px            '
c            else if (num1o.eq.3) then
c               str1='py            '
c            else if (num1o.eq.4) then
c               str1='pz            '
c            else if (num1o.eq.5) then
c               str1='dxy           '
c            else if (num1o.eq.6) then
c               str1='dyz           '
c            else if (num1o.eq.7) then
c               str1='dzx           '
c            else if (num1o.eq.8) then
c               str1='d(x**2-y**2)  '
c            else if (num1o.eq.9) then
c               str1='d(3z**2-r**2) '
c            else if (num1o.eq.10) then
c               str1='s*            '
c            end if
c            num2a=1+(num2-1)/10-nalamax
c            num2o=num2-(num2a-1)*10
c            if (num2o.eq.1) then
c               str2='s             '
c            else if (num2o.eq.2) then
c               str2='px            '
c            else if (num2o.eq.3) then
c               str2='py            '
c            else if (num2o.eq.4) then
c               str2='pz            '
c            else if (num2o.eq.5) then
c               str2='dxy           '
c            else if (num2o.eq.6) then
c               str2='dyz           '
c            else if (num2o.eq.7) then
c               str2='dzx           '
c            else if (num2o.eq.8) then
c               str2='d(x**2-y**2)  '
c            else if (num2o.eq.9) then
c               str2='d(3z**2-r**2) '
c            else if (num2o.eq.10) then
c               str2='s*            '
c            end if
c            write(69,*) num1a,str1,num2a,str2,val
c         end do
c      end do
c      close(68)
c      close(69)
cc
cc/////////////////////////////////////////////////
cc
c      ordhe=4*nalamax*10
c      open(unit=96,file='ordhe.dat',status='unknown')
c      write(96,*) ordhe
c      close(96)
cc
c      preorblae(1)=0
c      do i=1,3
c         preorblae(i+1)=preorblae(i)+nalamax*10
c      end do
cc
c      do ll=1,snumk+1
c         k=(ll-1)*dk
cc
c         do i=1,ordhe
c            do j=1,ordhe
c               hamie(i,j)=0.0d0
c            end do
c         end do
cc
c         do i=1,4
c            pre=preorblae(i)
c            do i1=1,nalamax*10
c               do j1=1,nalamax*10
c                  if ((i1.le.norb(i)).and.(j1.le.norb(i))) then
c                     hamie(pre+i1,pre+j1)=diag(i+4,i1,j1)
c                  else
c                     hamie(pre+i1,pre+j1)=0.0d0
c                  end if
c               end do
c            end do
c         end do
c         do i=1,3
c            pre=preorblae(i)
c            pre1=preorblae(i+1)
c            do i1=1,nalamax*10
c               do j1=1,nalamax*10
c                  if ((i1.le.norb(i)).and.(j1.le.norb(i+1))) then
c                     hamie(pre+i1,pre1+j1)=superdiag(i+4,i1,j1)
c                     hamie(pre1+j1,pre+i1)=subdiag(i+4,j1,i1)
c                  else
c                     hamie(pre+i1,pre1+j1)=0.0d0
c                     hamie(pre1+j1,pre+i1)=0.0d0
c                  end if
c               end do
c            end do
c         end do
c         pre=preorblae(4)
c         pre1=preorblae(1)
c         do i1=1,nalamax*10
c            do j1=1,nalamax*10
c               if ((i1.le.norb(4)).and.(j1.le.norb(1))) then
c                  hamie(pre+i1,pre1+j1)=superdiag(8,i1,j1)*exp(ui*k*a)
c                  hamie(pre1+j1,pre+i1)=subdiag(4,j1,i1)*exp(-ui*k*a)
c               else
c                  hamie(pre+i1,pre1+j1)=0.0d0
c                  hamie(pre1+j1,pre+i1)=0.0d0
c               end if
c            end do
c         end do
cc
cc        check if H is Hermitian
cc
c         open(unit=104,file='test4.dat',status='unknown',
c     #                                                  access='append')
c         do i=1,ordhe
c            do j=1,ordhe
c               if (hamie(i,j).ne.dconjg(hamie(j,i))) then
c                  write(104,*) ll,i,j
c                  write(104,*) 'i,j',dreal(hamie(i,j)),dimag(hamie(i,j))
c                  write(104,*) 'j,i',dreal(hamie(j,i)),dimag(hamie(j,i))
c                  write(104,*) 'diff',
c     #                             dreal(hamie(i,j)-dconjg(hamie(j,i))),
c     #                             dimag(hamie(i,j)-dconjg(hamie(j,i)))
c                  write(104,*) ''
c               end if
c            end do
c         end do
c         close(104)
c
c         call zheev('N','U',ordhe,hamie,4*natol*10,en,work,2*ordhe,
c     #               rwork,info)
cc        pay attention to put the correct leading dimension
cc        (i.e. the dimension of the variable declaration)!
c
c         do i=1,ordhe
c            if (i.lt.10) then
c               str1ae=char(ichar('0')+i)
c               nc=1
c            else
c               na1=i
c               na2=na1/10
c               rem=mod(na1,10)
c               str1ae=char(ichar('0')+rem)
c               nc=1
c               na1=na2
c 64            if (na2.lt.10) goto 74
c               na2=na1/10
c               rem=mod(na1,10)
c               str1ae=char(ichar('0')+rem)//str1ae(1:nc)
c               nc=nc+1
c               na1=na2
c               goto 64
c 74            str1ae=char(ichar('0')+na1)//str1ae(1:nc)
c               nc=nc+1
c            end if 
c            if (31.lt.(nc+9)) then
c                write(*,*) 'change 31: 31<(nc+9) !!'
c                goto 10
c            end if
c            str2ae='band'//str1ae(1:nc)//'e.dat'
c            open(unit=104,file=str2ae,status='unknown',access='append')
cc            write(104,*) k/(pi/a),en(i)/q
c            write(104,*) k/(pi/a),en(i)
cc           k in units of pi/a; energy in eV
c            close(104)
c         end do
c         if (ll.eq.1) then
c            open(unit=105,file='enk0e.dat',status='unknown')
c            do i=1,ordhe
c               write(105,*) en(i)
c            end do
c            close(105)
c         end if
cc
cc
c      end do
cc
cc
c      open(unit=106,file='bandse.dat',status='unknown')
c      do i=1,ordhe
c         if (i.lt.10) then
c            str1ae=char(ichar('0')+i)
c            nc=1
c         else
c            na1=i
c            na2=na1/10
c            rem=mod(na1,10)
c            str1ae=char(ichar('0')+rem)
c            nc=1
c            na1=na2
c 84         if (na2.lt.10) goto 94
c            na2=na1/10
c            rem=mod(na1,10)
c            str1ae=char(ichar('0')+rem)//str1ae(1:nc)
c            nc=nc+1
c            na1=na2
c            goto 84
c 94         str1ae=char(ichar('0')+na1)//str1ae(1:nc)
c            nc=nc+1
c         end if 
c         if (31.lt.(nc+9)) then
c            write(*,*) 'change 31: 31<(nc+9) !!'
c            goto 10
c         end if
c         str2ae='band'//str1ae(1:nc)//'e.dat'
c         open(unit=105,file=str2ae,status='unknown')
c         do ll=1,snumk+1
c            read(105,*) val1,val2
c            write(106,*) val1,val2
c         end do
c         close(105)
c         write(106,*) ' '
c      end do
c      close(106)
c
c/////////////////////////////////////////////////
c
 90   format(i3,a1,i2,a1,3x,i3,a1,i2,a1,3x,f18.13)
 91   format(a1)
 92   format(a1,i5,a2,i5)
 93   format(a2,f18.13)
 95   format(i5,a2,i5)
c
c
      return

 10   end
      
c      return
      
c      end

      
      subroutine matrix(*,nla,nalayer,ty,up1,tup1,up2,tup2,
     #do1,tdo1,do2,tdo2,es,ep,ew,ed,ss1,ww1,sw1,sp1,
     #wp1,sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,deltae,
     #pot,prevorb,norb,diag,
c     #subdiag,
     #superdiag,nlay,natol)
c
      implicit none
      integer nlay,natol
c     maximum number of layers, maximum number of atoms per layer 
      integer i,j,i1,j1,nla,nalayer(nlay),nl,natom,ty(nlay,natol),
     #up1(nlay,natol),up2(nlay,natol),tup1(nlay,natol),
     #tup2(nlay,natol),do1(nlay,natol),do2(nlay,natol),
     #tdo1(nlay,natol),tdo2(nlay,natol),
     #prevorb(nlay,natol),norb(nlay),pr,pr1,
     #upp1,upp2,doo1,doo2
      double precision es,ep,ew,ed,ss1,ww1,sw1,sp1,wp1,sd1,wd1,
     #pp1,pp2,pd1,pd2,dd1,dd2,dd3,l,m,n,invsq3,
     #h1(10,10),pot(nlay,natol),
     #diag(nlay,natol*10,natol*10),
c     #subdiag(nlay-1,natol*10,natol*10),
     #superdiag(nlay-1,natol*10,natol*10),
     #deltae,hpass(4,4)
c
c     in the Hamiltonian we put in the order the orbitals:
c     1:s, 2:px, 3:py, 4:pz, 5:dxy, 6:dyz, 7:dzx, 
c     8:d(x**2-y**2), 9:d(3z**2-r**2), 10:w
c
      invsq3=1.0d0/sqrt(3.0d0)
c
      do i=1,nla
         do i1=1,norb(i)
            do j1=1,norb(i)
               diag(i,i1,j1)=0.0d0
            end do
         end do
      end do
      do i=1,nla-1
         do i1=1,norb(i)
            do j1=1,norb(i+1)
               superdiag(i,i1,j1)=0.0d0
            end do
         end do
      end do
c      do i=1,nla-1
c         do i1=1,norb(i+1)
c            do j1=1,norb(i)
c               subdiag(i,i1,j1)=0.0d0
c            end do
c         end do
c      end do
c
c
      do i=1,nla
         do j=1,nalayer(i)
            pr=prevorb(i,j)
c
c
            upp1=up1(i,j)
            upp2=up2(i,j)
            doo1=do1(i,j)
            doo2=do2(i,j)
            diag(i,pr+1,pr+1)=es+pot(i,j)
            diag(i,pr+2,pr+2)=ep+pot(i,j)
            diag(i,pr+3,pr+3)=ep+pot(i,j)
            diag(i,pr+4,pr+4)=ep+pot(i,j)
            if ((upp1.lt.0).or.(upp2.lt.0).or.(doo1.lt.0)
     #         .or.(doo2.lt.0)) then
               call passivation(upp1,upp2,doo1,doo2,deltae,hpass)
               do i1=1,4
                  do j1=1,4
                     diag(i,pr+i1,pr+j1)=
     #                  diag(i,pr+i1,pr+j1)+hpass(i1,j1)
                  end do
               end do
            end if
            diag(i,pr+5,pr+5)=ed+pot(i,j)
            diag(i,pr+6,pr+6)=ed+pot(i,j)
            diag(i,pr+7,pr+7)=ed+pot(i,j)
            diag(i,pr+8,pr+8)=ed+pot(i,j)
            diag(i,pr+9,pr+9)=ed+pot(i,j)
            diag(i,pr+10,pr+10)=ew+pot(i,j)
c
c
            if (up1(i,j).gt.0) then
               pr1=prevorb(i+1,up1(i,j))
               if (ty(i,j).eq.1) then
c                 nn-type=1 with respect to an anion
                  l=invsq3
                  m=invsq3
                  n=invsq3
               else
c                 nn-type=1 with respect to a cation
                  l=-invsq3
                  m=invsq3
                  n=invsq3
               end if
               call elementsitosi(l,m,n,ss1,ww1,sw1,sp1,wp1,
     #              sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1)
               do i1=1,10
                  do j1=1,10
                     superdiag(i,pr+i1,pr1+j1)=h1(i1,j1)
                  end do
               end do
            end if
c
c
            if (up2(i,j).gt.0) then
               pr1=prevorb(i+1,up2(i,j))
               if (ty(i,j).eq.1) then
c                 nn-type=2 with respect to an anion
                  l=-invsq3
                  m=-invsq3
                  n=invsq3
               else
c                 nn-type=2 with respect to a cation
                  l=invsq3
                  m=-invsq3
                  n=invsq3
               end if
               call elementsitosi(l,m,n,ss1,ww1,sw1,sp1,wp1,
     #              sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1)
               do i1=1,10
                  do j1=1,10
                     superdiag(i,pr+i1,pr1+j1)=h1(i1,j1)
                  end do
               end do
            end if
cc
cc
c            if (do1(i,j).gt.0) then
c               pr1=prevorb(i-1,do1(i,j))
c               if (ty(i,j).eq.1) then
cc                 nn-type=3 with respect to an anion
c                  l=invsq3
c                  m=-invsq3
c                  n=-invsq3
c               else
cc                 nn-type=3 with respect to a cation
c                  l=-invsq3
c                  m=-invsq3
c                  n=-invsq3
c               end if
c               call elementsitosi(l,m,n,ss1,ww1,sw1,sp1,wp1,
c     #              sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1)
c               do i1=1,10
c                  do j1=1,10
c                     subdiag(i-1,pr+i1,pr1+j1)=h1(i1,j1)
c                  end do
c               end do
c            end if
cc
cc
c            if (do2(i,j).gt.0) then
c               pr1=prevorb(i-1,do2(i,j))
c               if (ty(i,j).eq.1) then
cc                 nn-type=4 with respect to an anion
c                  l=-invsq3
c                  m=invsq3
c                  n=-invsq3
c               else
cc                 nn-type=4 with respect to a cation
c                  l=invsq3
c                  m=invsq3
c                  n=-invsq3
c               end if
c               call elementsitosi(l,m,n,ss1,ww1,sw1,sp1,wp1,
c     #              sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1)
c               do i1=1,10
c                  do j1=1,10
c                     subdiag(i-1,pr+i1,pr1+j1)=h1(i1,j1)
c                  end do
c               end do
c            end if
cc
cc
         end do
      end do
cc
cc     check if H is symmetric
cc
c      open(unit=14,file='test.dat',status='unknown')
c      do i=1,nla
c         do i1=1,norb(i)
c            do j1=1,norb(i)
c               if (diag(i,i1,j1).ne.diag(i,j1,i1)) then
c                  write(14,*) 'dr',i,i1,j1,
c     #                                diag(i,i1,j1)-diag(i,j1,i1)
c               end if
c            end do
c         end do
c      end do
c      do i=1,nla-1
c         do i1=1,norb(i)
c            do j1=1,norb(i+1)
c               if (superdiag(i,i1,j1).ne.subdiag(i,j1,i1)) then
c                  write(14,*) 's',i,i1,j1,
c     #                            superdiag(i,i1,j1)-subdiag(i,j1,i1)
c               end if
c            end do
c         end do
c      end do
c      close(14)
cc
cc
      end


      subroutine passivation(upp1,upp2,doo1,doo2,deltae,hpass)
c
      implicit none
      integer i,j,upp1,upp2,doo1,doo2,ancat,ind,l1,l2,l3,l4
      double precision deltae,hpass(4,4),const
c
      l1=0
      l2=0
      l3=0
      l4=0
c
      if (upp1.eq.-2) then
c        dangling bond in the direction [1 1 1] with respect to an anion
         ancat=1
         l1=1
      end if
      if (upp2.eq.-2) then
c        dangling bond in the direction [-1 -1 1] with respect to an anion
         ancat=1
         l2=1
      end if
      if (doo1.eq.-2) then
c        dangling bond in the direction [1 -1 -1] with respect to an anion
         ancat=1
         l3=1
      end if
      if (doo2.eq.-2) then
c        dangling bond in the direction [-1 1 -1] with respect to an anion
         ancat=1
         l4=1
      end if
c
      if (upp1.eq.-1) then
c        dangling bond in the direction [-1 1 1] with respect to a cation
         ancat=2
         l1=1
      end if
      if (upp2.eq.-1) then
c        dangling bond in the direction [1 -1 1] with respect to a cation
         ancat=2
         l2=1
      end if
      if (doo1.eq.-1) then
c        dangling bond in the direction [-1 -1 -1] with respect to a cation
         ancat=2
         l3=1
      end if
      if (doo2.eq.-1) then
c        dangling bond in the direction [1 1 -1] with respect to a cation
         ancat=2
         l4=1
      end if
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


c
c     1 is s,
c     2 is px,
c     3 is py,
c     4 is pz,
c     5 is dxy,
c     6 is dyz,
c     7 is dzx,
c     8 is d(x**2-y**2),
c     9 is d(3z**2-r**2),
c     10 is s* (i.e. w)
c


      subroutine elementsitosi(l,m,n,ss1,ww1,sw1,sp1,wp1,
     #sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1)
c     Hamiltonian elements from a Si atom to a Si atom
      implicit none
      double precision l,m,n,ss1,ww1,sw1,sp1,wp1,
     #sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1(10,10),sq3,ha,qu
c
      sq3=sqrt(3.0d0)
      ha=1/2.0d0
      qu=1/4.0d0
c
      h1(1,1)=ss1
      h1(1,2)=l*sp1
      h1(1,3)=m*sp1
      h1(1,4)=n*sp1
      h1(1,5)=sq3*l*m*sd1
      h1(1,6)=sq3*m*n*sd1
      h1(1,7)=sq3*n*l*sd1
      h1(1,8)=ha*sq3*(l*l-m*m)*sd1
      h1(1,9)=(n*n-ha*(l*l+m*m))*sd1
      h1(1,10)=sw1
c
      h1(2,1)=-l*sp1
      h1(2,2)=l*l*pp1+(1-l*l)*pp2
      h1(2,3)=l*m*(pp1-pp2)
      h1(2,4)=n*l*(pp1-pp2)
      h1(2,5)=sq3*l*l*m*pd1+m*(1-2*l*l)*pd2
      h1(2,6)=l*m*n*(sq3*pd1-2*pd2)
      h1(2,7)=sq3*n*l*l*pd1+n*(1-2*l*l)*pd2
      h1(2,8)=ha*sq3*l*(l*l-m*m)*pd1+l*(1-l*l+m*m)*pd2
      h1(2,9)=l*(n*n-ha*(l*l+m*m))*pd1-sq3*l*n*n*pd2
      h1(2,10)=-l*wp1
c
      h1(3,1)=-m*sp1
      h1(3,2)=m*l*(pp1-pp2)
      h1(3,3)=m*m*pp1+(1-m*m)*pp2
      h1(3,4)=m*n*(pp1-pp2)
      h1(3,5)=sq3*l*m*m*pd1+l*(1-2*m*m)*pd2
      h1(3,6)=sq3*n*m*m*pd1+n*(1-2*m*m)*pd2
      h1(3,7)=l*m*n*(sq3*pd1-2*pd2)
      h1(3,8)=ha*sq3*m*(l*l-m*m)*pd1-m*(1+l*l-m*m)*pd2
      h1(3,9)=m*(n*n-ha*(l*l+m*m))*pd1-sq3*m*n*n*pd2
      h1(3,10)=-m*wp1
c
      h1(4,1)=-n*sp1
      h1(4,2)=n*l*(pp1-pp2)
      h1(4,3)=n*m*(pp1-pp2)
      h1(4,4)=n*n*pp1+(1-n*n)*pp2
      h1(4,5)=l*m*n*(sq3*pd1-2*pd2)
      h1(4,6)=sq3*m*n*n*pd1+m*(1-2*n*n)*pd2
      h1(4,7)=sq3*l*n*n*pd1+l*(1-2*n*n)*pd2
      h1(4,8)=ha*sq3*n*(l*l-m*m)*pd1-n*(l*l-m*m)*pd2
      h1(4,9)=n*(n*n-ha*(l*l+m*m))*pd1+sq3*n*(l*l+m*m)*pd2
      h1(4,10)=-n*wp1
c
      h1(5,1)=sq3*l*m*sd1
      h1(5,2)=-sq3*l*l*m*pd1-m*(1-2*l*l)*pd2
      h1(5,3)=-sq3*l*m*m*pd1-l*(1-2*m*m)*pd2
      h1(5,4)=-l*m*n*(sq3*pd1-2*pd2)
      h1(5,5)=3*l*l*m*m*dd1+(l*l+m*m-4*l*l*m*m)*dd2+(n*n+l*l*m*m)*dd3
      h1(5,6)=3*l*m*m*n*dd1+n*l*(1-4*m*m)*dd2+n*l*(m*m-1)*dd3
      h1(5,7)=3*l*l*m*n*dd1+m*n*(1-4*l*l)*dd2+m*n*(l*l-1)*dd3
      h1(5,8)=3*ha*l*m*(l*l-m*m)*dd1+2*l*m*(m*m-l*l)*dd2+
     #                                         ha*l*m*(l*l-m*m)*dd3
      h1(5,9)=sq3*l*m*(n*n-ha*(l*l+m*m))*dd1-2*sq3*l*m*n*n*dd2+
     #                                       ha*sq3*l*m*(1+n*n)*dd3
      h1(5,10)=sq3*l*m*wd1
c
      h1(6,1)=sq3*m*n*sd1
      h1(6,2)=-l*m*n*(sq3*pd1-2*pd2)
      h1(6,3)=-sq3*n*m*m*pd1-n*(1-2*m*m)*pd2
      h1(6,4)=-sq3*m*n*n*pd1-m*(1-2*n*n)*pd2
      h1(6,5)=3*l*m*m*n*dd1+n*l*(1-4*m*m)*dd2+n*l*(m*m-1)*dd3
      h1(6,6)=3*m*m*n*n*dd1+(m*m+n*n-4*m*m*n*n)*dd2+(l*l+m*m*n*n)*dd3
      h1(6,7)=3*l*m*n*n*dd1+l*m*(1-4*n*n)*dd2+l*m*(n*n-1)*dd3
      h1(6,8)=3*ha*m*n*(l*l-m*m)*dd1-m*n*(1+2*(l*l-m*m))*dd2+
     #                                     m*n*(1+ha*(l*l-m*m))*dd3
      h1(6,9)=sq3*m*n*(n*n-ha*(l*l+m*m))*dd1+sq3*m*n*(l*l+m*m-n*n)*dd2-
     #                                     ha*sq3*m*n*(l*l+m*m)*dd3
      h1(6,10)=sq3*m*n*wd1
c
      h1(7,1)=sq3*n*l*sd1
      h1(7,2)=-sq3*l*l*n*pd1-n*(1-2*l*l)*pd2
      h1(7,3)=-l*m*n*(sq3*pd1-2*pd2)
      h1(7,4)=-sq3*l*n*n*pd1-l*(1-2*n*n)*pd2
      h1(7,5)=3*l*l*m*n*dd1+m*n*(1-4*l*l)*dd2+m*n*(l*l-1)*dd3
      h1(7,6)=3*l*m*n*n*dd1+l*m*(1-4*n*n)*dd2+l*m*(n*n-1)*dd3
      h1(7,7)=3*l*l*n*n*dd1+(l*l+n*n-4*l*l*n*n)*dd2+(m*m+l*l*n*n)*dd3
      h1(7,8)=3*ha*n*l*(l*l-m*m)*dd1+n*l*(1-2*(l*l-m*m))*dd2-
     #                                     n*l*(1-ha*(l*l-m*m))*dd3
      h1(7,9)=sq3*n*l*(n*n-ha*(l*l+m*m))*dd1+sq3*n*l*(l*l+m*m-n*n)*dd2-
     #                                     ha*sq3*n*l*(l*l+m*m)*dd3
      h1(7,10)=sq3*n*l*wd1
c
      h1(8,1)=ha*sq3*(l*l-m*m)*sd1
      h1(8,2)=-ha*sq3*l*(l*l-m*m)*pd1-l*(1-l*l+m*m)*pd2
      h1(8,3)=-ha*sq3*m*(l*l-m*m)*pd1+m*(1+l*l-m*m)*pd2
      h1(8,4)=-ha*sq3*n*(l*l-m*m)*pd1+n*(l*l-m*m)*pd2
      h1(8,5)=3*ha*l*m*(l*l-m*m)*dd1+2*l*m*(m*m-l*l)*dd2+
     #                                         ha*l*m*(l*l-m*m)*dd3
      h1(8,6)=3*ha*m*n*(l*l-m*m)*dd1-m*n*(1+2*(l*l-m*m))*dd2+
     #                                     m*n*(1+ha*(l*l-m*m))*dd3
      h1(8,7)=3*ha*n*l*(l*l-m*m)*dd1+n*l*(1-2*(l*l-m*m))*dd2-
     #                                     n*l*(1-ha*(l*l-m*m))*dd3
      h1(8,8)=3*qu*((l*l-m*m)**2)*dd1+(l*l+m*m-((l*l-m*m)**2))*dd2+
     #                                  (n*n+qu*((l*l-m*m)**2))*dd3
      h1(8,9)=ha*sq3*(l*l-m*m)*(n*n-ha*(l*l+m*m))*dd1+
     #           sq3*n*n*(m*m-l*l)*dd2+qu*sq3*(1+n*n)*(l*l-m*m)*dd3
      h1(8,10)=ha*sq3*(l*l-m*m)*wd1
c
      h1(9,1)=(n*n-ha*(l*l+m*m))*sd1
      h1(9,2)=-l*(n*n-ha*(l*l+m*m))*pd1+sq3*l*n*n*pd2
      h1(9,3)=-m*(n*n-ha*(l*l+m*m))*pd1+sq3*m*n*n*pd2
      h1(9,4)=-n*(n*n-ha*(l*l+m*m))*pd1-sq3*n*(l*l+m*m)*pd2
      h1(9,5)=sq3*l*m*(n*n-ha*(l*l+m*m))*dd1-2*sq3*l*m*n*n*dd2+
     #                                       ha*sq3*l*m*(1+n*n)*dd3
      h1(9,6)=sq3*m*n*(n*n-ha*(l*l+m*m))*dd1+sq3*m*n*(l*l+m*m-n*n)*dd2-
     #                                     ha*sq3*m*n*(l*l+m*m)*dd3
      h1(9,7)=sq3*n*l*(n*n-ha*(l*l+m*m))*dd1+sq3*n*l*(l*l+m*m-n*n)*dd2-
     #                                     ha*sq3*n*l*(l*l+m*m)*dd3
      h1(9,8)=ha*sq3*(l*l-m*m)*(n*n-ha*(l*l+m*m))*dd1+
     #           sq3*n*n*(m*m-l*l)*dd2+qu*sq3*(1+n*n)*(l*l-m*m)*dd3
      h1(9,9)=((n*n-ha*(l*l+m*m))**2)*dd1+3*n*n*(l*l+m*m)*dd2+
     #                                      3*qu*((l*l+m*m)**2)*dd3
      h1(9,10)=(n*n-ha*(l*l+m*m))*wd1
c
      h1(10,1)=sw1
      h1(10,2)=l*wp1
      h1(10,3)=m*wp1
      h1(10,4)=n*wp1
      h1(10,5)=sq3*l*m*wd1
      h1(10,6)=sq3*m*n*wd1
      h1(10,7)=sq3*n*l*wd1
      h1(10,8)=ha*sq3*(l*l-m*m)*wd1
      h1(10,9)=(n*n-ha*(l*l+m*m))*wd1
      h1(10,10)=ww1 
c
      end
