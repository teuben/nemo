/* 
 * awk '{print $7}' $NEMO/adm/Usage | ustat - | tabhist -
 * 
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in=???\n       Input file in format hh:mm:ss",
    "VERSION=0\n    18-nov-90 PJT",
    NULL,
};

nemo_main()
{
    stream instr;
    char line[80];
    real t, t0;
    int  hh, mm, ss;
    bool first = TRUE;


    instr = stropen(getparam("in"), "r");
    while (fgets(line,80,instr) != NULL) {
        hh = atoi(&line[0]);
        mm = atoi(&line[3]);
        ss = atoi(&line[6]);
        t = (hh*60 + mm) + ss/60.0;
        if (!first) {
            if (t-t0 < 0)
               printf("%6.2f\n",t-t0+24*60);
            else
               printf("%6.2f\n",t-t0);
            t0 = t;
        } else
            first = FALSE;
    }
}    
