#include <bodytrans.h>

real btr_jtot(b,t,i)
Body *b;
real t;
int  i;
{
    return (sqrt(sqr(x*vy - y*vx) + sqr(y*vz - z*vy) + sqr(z*vx - x*vz)));
}
