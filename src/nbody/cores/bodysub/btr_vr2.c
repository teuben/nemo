#include <bodytrans.h>

real btr_vr2(Body *b,real t,int  i)
{
    return ((x*vx + y*vy) / sqrt(x*x + y*y));
}
