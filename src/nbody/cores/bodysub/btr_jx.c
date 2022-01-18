#include <bodytrans.h>

real btr_jx(Body *b,real t,int  i)
{
  return y*vz - z*vy;
}
