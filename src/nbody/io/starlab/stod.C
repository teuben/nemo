/*
 *  STOD: convert snapshot to dyn
 *
 *  16-jul-93   1.0     written
 *  31-jul-96   1.1     for starlab 2.1
 *  18-jul-01   1.2     for starlab 3.6:  allow to read more times
 *  14-oct-03   1.3     for starlab 4: we never kept up with the changes
 *  
 *  This has to be a C++ main, otherwise won't link ???
 *      - nemomain.c ported to nemomain.C (very small mod)
 *	- stod_subs.c is the accompanying snapshot I/O routines in C
 *
 *  15-oct-03   1.3    fixed for starlab4             PJT
 *   2-jan-03   1.4    attempt to integrated getting UBVRI for nvodemo2004
 *  30-jul-05   2.0    adapted for new StarLab V4.x with new dir.str.
 */


#include <nemo_stdinc.h>                 /* NEMO */
#include <getparam.h>
#include <history.h>
#include <extstring.h>

#include "stod_subs.h"              /* NEMO-STARLAB interface */

#include "dyn.h"                    /* STARLAB */


nemo_string defv[] = {
    "in=???\n           Input snapshot file (dyn to stdout)",
    "headline=\n        Additional comment line",
    "label=f\n          Add labels to stars?",
    "VERSION=2.0\n	30-jul-05 PJT",
    NULL,
};

nemo_string usage = "convert NEMO snapshot to STARLAB dyn";

void nemo_main(void)
{
    int i, nbody, hisc;
    double *mass, *pos, *vel, *mptr, *pptr, *vptr, tsnap;
    dyn *b = new dyn();
    dyn *bo = new dyn();
    dyn *by, *bi;
    bool i_flag = getbparam("label");
    nemo_string *hisv, headline = getparam("headline");

    check_real(sizeof(real));

    nbody = get_snap_c(getparam("in"), &tsnap, &mass, &pos, &vel);

    if (i_flag) bo->set_label(1);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    for (i = 1; i < nbody; i++) {
        by = new dyn();
        if (i_flag) by->set_label(i+1);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
        bo = by;
    }

    if (*headline == 0) 
        headline = ask_headline();
    if (headline && *headline)
        b->log_comment(headline);
    b->log_comment(" ### This is still an experimental conversion ### ");
    hisv = ask_history();
    hisc = xstrlen(hisv,sizeof(nemo_string)) - 1;
    b->log_history(hisc,hisv);


// Do not store top level info - BUG? does seem to be created anyways.
//
//    b->set_mass(1); 

    mptr = mass;    // set pointers 
    pptr = pos;
    vptr = vel;
    for (bi = b->get_oldest_daughter(); bi != NULL;     // loop over stars
         bi = bi->get_younger_sister()) {
	bi->set_mass(*mptr);	
        bi->set_pos(vec(pptr[0],pptr[1],pptr[2]));
        bi->set_vel(vec(vptr[0],vptr[1],vptr[2]));
        mptr++;
        pptr += 3;  // NDIM really
        vptr += 3;
    }
    put_dyn(b);
}
