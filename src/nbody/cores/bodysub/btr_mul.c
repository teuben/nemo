#include <bodytrans.h>

real btr_mul(Body *b,real t,int  i)
{
    real wr2 = x*x+y*y;
    real wr = sqrt(wr2+z*z);
    wr2 = sqrt(wr2);
    
    if (wr2 > 0 && wr > 0)
        return  (x*vy-y*vx)/(wr2*wr);   /*  /4.74 ??? */
    else
        return 0.0;

}
