q=1.60219e-19;
t=2.7*q
acc=0.144
a=sqrt(3)*acc
#V=(0.2023-0.1317)*q;
V=0.6*q;
tp=0.35*q
#set xrange[-2e10:-1.5e10]
set xrange[0:12.6]
ky=pi/a

#q(x,y)=t**2*(4*(cos(a*y/2))**2+4*cos(sqrt(3)*acc*y/2)*cos(sqrt(3)*a*x/2) +1)
#Egraphene1(x)=sqrt(q(x,y))/q
#E1(x,y)=(V/2+sqrt(q(x,y)+V**2/4+tp**2/2+0.5*sqrt(4*(V**2+tp**2)*q(x,y)+tp**4)))/q
#E2(x,y)=(V/2-sqrt(q(x,y)+V**2/4+tp**2/2+0.5*sqrt(4*(V**2+tp**2)*q(x,y)+tp**4)))/q
#E3(x,y)=(V/2-sqrt(q(x,y)+V**2/4+tp**2/2-0.5*sqrt(4*(V**2+tp**2)*q(x,y)+tp**4)))/q
#E4(x,y)=(V/2+sqrt(q(x,y)+V**2/4+tp**2/2-0.5*sqrt(4*(V**2+tp**2)*q(x,y)+tp**4)))/q

#q(x)=t**2*(4*(cos(a*x/2))**2+4*cos(sqrt(3)*acc*x/2) +1)
q(x)=t**2*(2*(cos(a*x/2))+1)**2

E1(x)=(+sqrt(q(x)+V**2/4+tp**2/2+0.5*sqrt(4*(V**2+tp**2)*q(x)+tp**4)))/q
E2(x)=(-sqrt(q(x)+V**2/4+tp**2/2+0.5*sqrt(4*(V**2+tp**2)*q(x)+tp**4)))/q
E3(x)=(-sqrt(q(x)+V**2/4+tp**2/2-0.5*sqrt(4*(V**2+tp**2)*q(x)+tp**4)))/q
E4(x)=(+sqrt(q(x)+V**2/4+tp**2/2-0.5*sqrt(4*(V**2+tp**2)*q(x)+tp**4)))/q
plot E4(x) lt 1 w l,E1(x) lt 1 w l,E4(-x-2*ky) lt 2 w l,E1(-x-2*ky) lt 2 w l
