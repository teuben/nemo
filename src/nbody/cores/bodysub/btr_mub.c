#include <bodytrans.h>

real btr_mub(Body *b,real t,int  i)
{
    real r2 = x*x+y*y;
    real r = sqrt(r2+z*z);
    r2 = sqrt(r2);
    
    if (r2 > 0 && r > 0)
        return  (vz*r2-z*(y*vy+x*vx)/r2)/(r*r);	/*  /4.74 ??? */
    else
        return 0.0;

}
