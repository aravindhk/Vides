#include "nanotube.h"
void nanotube(int i,int j,double *x,double *y,double *z,int n,double acc)
{
  int indz;
  double rt,angle;

  angle=2*pi/n;

  rt=sqrt(3.00)*acc*n/pi*0.5;

  if (((j%4)==0)||((j%4)==1))
    {
      *x=rt*cos(i*angle);
      *y=rt*sin(i*angle);
      if (fabs(*x)<1e-8)
	*x=0;
      if (fabs(*y)<1e-8)
	*y=0;
    }
  else
    {
      *x=rt*cos((i+0.5)*angle);
      *y=rt*sin((i+0.5)*angle);
      if (fabs(*x)<1e-8)
	*x=0;
      if (fabs(*y)<1e-8)
	*y=0;
    }
  
  indz=floor(j/4);
  if ((j%4)==0)
    *z=indz*3*acc;
  else if ((j%4)==1)
    *z=indz*3*acc+acc;
  else if ((j%4)==2)
    *z=indz*3*acc+1.5*acc;
  else if ((j%4)==3)
    *z=indz*3*acc+2.5*acc;
  return;
}
