#include <bodytrans.h>

real btr_v(b,t,i)
Body *b;
real t;
int  i;
{
    return (sqrt(vx*vx + vy*vy + vz*vz));
}
