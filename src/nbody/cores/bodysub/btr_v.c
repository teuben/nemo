#include <bodytrans.h>

real btr_v(Body *b,real t,int  i)
{
    return (sqrt(vx*vx + vy*vy + vz*vz));
}
