/* return angles, and keep track of the wrap. Only works for slowly
 * varying angles, i.e. < 45 degrees per call
 */

#include <bodytrans.h>

static int wrap = 0;
static int first = 1;
static real last;

real btr_ang(b,t,i)
Body *b;
real t;
int  i;
{
    real now = atan2(y,x)*57.2958;

    dprintf(2,"last=%g now=%g\n",last,now);

    if (first) {
        first = 0;
        last = now;
        return now;
    }

    if (now*last < 0 && ABS(now) > 90) {
        dprintf(1,"wrap=%d last=%g now=%g\n",wrap,last,now);
        if (now > 0)
            wrap--;
        else
            wrap++;
            
    }
    last = now;	
    return now + wrap*360.0;
}
