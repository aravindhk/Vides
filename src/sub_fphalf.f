      double precision function sub_fphalf(etaf)

c       This program calculates the definite integral F_(1/2)(etaf).
c
c       etaf = -energy relative to the Fermi level, in units of eV*q/kT

      implicit double precision (a-h,o-z)
      
      dimension a1(12), b1(12), a2(12), b2(12)
      parameter(expmax=69)
      expcap(x)=exp(max(min(x,expmax),-expmax))
        
      data m1,k1,m2,k2/7,7,10,11/
      data (a1(i),i=1,8) /5.75834152995465d+06,
     &                      1.30964880355883d+07,
     &                      1.07608632249013d+07,
     &                      3.93536421893014d+06,
     &                      6.42493233715640d+05,
     &                      4.16031909245777d+04,
     &                      7.77238678539648d+02,
     &                      1.00000000000000d+00/
      data (b1(i),i=1,8) /6.49759261942269d+06,
     &                      1.70750501625775d+07,
     &                      1.69288134856160d+07,
     &                      7.95192647756086d+06,
     &                      1.83167424554505d+06,
     &                      1.95155948326832d+05,
     &                      8.17922106644547d+03,
     &                      9.02129136642157d+01/
      data (a2(i),i=1,11)/4.85378381173415d-14,
     &                      1.64429113030738d-11,
     &                      3.76794942277806d-09,
     &                      4.69233883900644d-07,
     &                      3.40679845803144d-05,
     &                      1.32212995937796d-03,
     &                      2.60768398973913d-02,
     &                      2.48653216266227d-01,
     &                      1.08037861921488d+00,
     &                      1.91247528779676d+00,
     &                      1.00000000000000d+00/
      data (b2(i),i=1,12)/7.28067571760518d-14,
     &                      2.45745452167585d-11,
     &                      5.62152894375277d-09,
     &                      6.96888634549649d-07,
     &                      5.02360015186394d-05,
     &                      1.92040136756592d-03,
     &                      3.66887808002874d-02,
     &                      3.24095226486468d-01,
     &                      1.16434871200131d+00,
     &                      1.34981244060549d+00,
     &                      2.01311836975930d-01,
     &                     -2.14562434782759d-02/

      if (etaf .lt. 2) then
         x=expcap(etaf)
         rn=x+a1(m1)
         do i=m1-1,1,-1
            rn=rn*x+a1(i)
         end do
         den=b1(k1+1)
         do i=k1,1,-1
            den=den*x+b1(i)
         end do
         sub_fphalf=x*rn/den
      else
         x=1/(etaf*etaf)
         rn=x+a2(m2)
         do i=m2-1,1,-1
            rn=rn*x+a2(i)
         end do
         den=b2(k2+1)
         do i=k2,1,-1
            den=den*x+b2(i)
         end do
         sub_fphalf=etaf*sqrt(etaf)*rn/den
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
