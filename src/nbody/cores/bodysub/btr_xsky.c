#include <bodytrans.h>

real btr_xsky(Body *b,real t,int  i)
{
  return atan2(x,z)*180.0/PI;    /* astronomical coord system that observers use */
}
