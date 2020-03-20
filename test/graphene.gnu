acc=0.144;
a=acc*sqrt(3)
t=2.7

y=30
set xrange[-20:20]
f(x,y)=t*sqrt(1+4*cos(sqrt(3)*x*a/2)*cos(y*a/2)+4*(cos(y*a/2))**2)
#g(x)=t*sqrt(1+4*cos(sqrt(3)*x*a/2)+4)
plot f(x,y),-f(x,y)
