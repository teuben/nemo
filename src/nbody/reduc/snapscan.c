/*
 * SNAPSCAN.C:   scan snapshots and report times
 *		merely created for benchmarking times to scan nbody files
 *
 *	7-jun-94  1.0	created				pjt
 *     12-jan-99  1.1   make sure embedded history done pjt
 *                      and add option to look at all data  pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n              Input snapshot",
    "realloc=f\n           Realloc snapshot each read?",
    "read=t\n              Force reading all data ",
    "VERSION=1.1\n         12-jan-98 PJT",
    NULL,
};

string usage = "Scan a snapshot";


nemo_main()
{
    stream instr;
    Body   *btab = NULL;
    int    i, nbody, bits;
    bool Qrealloc;
    real tsnap;

    Qrealloc = getbparam("realloc");
    instr = stropen(getparam("in"), "r");           /* open input file */
    get_history(instr);                         /* get history */
    printf("# SnapShot Nbody Time Bits\n");
    for (i=0;;i++) {                           /* loop through snapshots */
        get_history(instr);                         /* get history */
        if (!get_tag_ok(instr, SnapShotTag))
                break;                           /* until done */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        printf("%d %d %g 0x%x\n", i+1, nbody, tsnap, bits);
        if (Qrealloc) {
            free(btab);
            btab = NULL;
        }
    }
}

