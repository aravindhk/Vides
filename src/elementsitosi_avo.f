      subroutine elementsitosiavo(l,m,n,ss1,ww1,sw1,sp1,wp1,
     #sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h)
c     Hamiltonian elements from a Si atom to a Si atom
      implicit none
      double precision l,m,n,ss1,ww1,sw1,sp1,wp1,
     #sd1,wd1,pp1,pp2,pd1,pd2,dd1,dd2,dd3,h1(10,10),sq3,ha,qu,h(100)
      integer i,j,ix
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

      do i=1,10
         do j=1,10
            ix=j+(i-1)*10
            h(ix)=h1(i,j)
         enddo
      enddo
c     
      end
