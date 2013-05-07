#include <bodytrans.h>

real btr_ysky(Body *b,real t,int  i)
{
    return atan2(y,z)*180.0/PI;    /* astronomical coord system that observers use */
}
