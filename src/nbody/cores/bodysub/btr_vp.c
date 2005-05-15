#include <bodytrans.h>

real btr_vp(Body *b,real t,int  i)
{
    real r2 = x*x + y*y + z*z;
    return sqrt(
	         ((vx*vx + vy*vy + vz*vz) - sqr(x*vx + y*vy + z*vz) / r2)/r2
               );
}
