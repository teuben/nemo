#include <bodytrans.h>

real btr_ar(Body *b,real t,int  i)
{
    return ((x*ax + y*ay + z*az) / sqrt(x*x + y*y + z*z));
}
