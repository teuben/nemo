#include <bodytrans.h>

real btr_glon(Body *b,real t,int  i)
{
    return atan2(y,x)*180.0/PI;    /* GLON is an astronomical coord system that observers use */
}
