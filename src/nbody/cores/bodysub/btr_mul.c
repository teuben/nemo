#include <bodytrans.h>

real btr_mul(Body *b,real t,int  i)
{
    real r2 = x*x+y*y;
    real r = sqrt(r2+z*z);
    r2 = sqrt(r2);
    
    if (r2 > 0 && r > 0)
        return  (x*vy-y*vx)/(r2*r);   /*  /4.74 ??? */
    else
        return 0.0;

}
