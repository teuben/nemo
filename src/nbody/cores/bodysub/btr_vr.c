#include <bodytrans.h>

real btr_vr(Body *b,real t,int  i)
{
    return ((x*vx + y*vy + z*vz) / sqrt(x*x + y*y + z*z));
}
