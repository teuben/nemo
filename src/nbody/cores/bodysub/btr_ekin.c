#include <bodytrans.h>

real btr_ekin(Body *b,real t,int  i)
{
    return 0.5*(vx*vx + vy*vy + vz*vz);
}
