#include <bodytrans.h>

real btr_vr(b,t,i)
Body *b;
real t;
int  i;
{
    return ((x*vx + y*vy + z*vz) / sqrt(x*x + y*y + z*z));
}
