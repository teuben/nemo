#include <bodytrans.h>

real btr_etot(b,t,i)
Body *b;
real t;
int  i;
{
    return (phi + 0.5*(vx*vx + vy*vy + vz*vz));
}
