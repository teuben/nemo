#include <bodytrans.h>

real btr_r(b,t,i)
Body *b;
real t;
int  i;
{
    return (sqrt(x*x + y*y + z*z));
}
