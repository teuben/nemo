#include <bodytrans.h>

real btr_vt(Body *b,real t,int  i)
{
    return (sqrt((vx*vx + vy*vy + vz*vz) -
		   sqr(x*vx + y*vy + z*vz) / (x*x + y*y + z*z)));
}
