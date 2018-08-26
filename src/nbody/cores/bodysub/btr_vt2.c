#include <bodytrans.h>

real btr_vt2(Body *b,real t,int  i)
{
    return (sqrt((vx*vx + vy*vy) -
		   sqr(x*vx + y*vy) / (x*x + y*y)));
}
