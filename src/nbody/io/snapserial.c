/*
 *  snapserial()
 *
 *  each time after a snapshot will be processed, call snapserial()
 *  if you want to test for serial processing of large datasets.
 *  some programs that need all the data, cannot deal with split
 *  datasets.
 */

#include <stdinc.h>

int snapserial(real tsnap, bool serial)
{
    permanent real last_tsnap;
    permanent bool reset=TRUE;

    if (reset) {
        last_tsnap = tsnap;
        reset = FALSE;
    } else {
        if (last_tsnap < tsnap)
            last_tsnap = tsnap;
        else if (last_tsnap > tsnap)
            reset = TRUE;
        else if (!serial)
            warning("time %g should not be treated as a serial process",tsnap);

    }
    return 0;
}
