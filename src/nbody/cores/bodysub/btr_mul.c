#include <bodytrans.h>

real btr_mul(b,t,i)
Body *b;
real t;
int  i;
{
    real r2 = x*x+y*y;
    real r = sqrt(r2+z*z);
    r2 = sqrt(r2);
    
    if (r2 > 0)
        return  (x*vy-y*vx)/r2;
    else
        return 0.0;

}
