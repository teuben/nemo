#include <bodytrans.h>

real btr_etot(Body *b,real t,int  i)
{
    return (phi + 0.5*(vx*vx + vy*vy + vz*vz));
}
