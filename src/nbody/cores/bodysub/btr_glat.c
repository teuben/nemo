#include <bodytrans.h>

real btr_glat(Body *b,real t,int  i)
{
    return atan2(z,sqrt(x*x+y*y))*180.0/PI;
}
