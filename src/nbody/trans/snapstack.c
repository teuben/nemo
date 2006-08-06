/*
 * SNAPSTACK: stack two N-body systems on top of each other, optionally
 *            shifting them in phase-space w.r.t. each other
 *
 *	22-jan-89  1.1 JEB
 *	 1-jul-90  1.1a added helpvec PJT
 *	 4-feb-93  1.1b nemo_main (also to return 0 exit)  -- PJT
 *			malloc -> allocate
 *	mar94 - ansi
 *	 6-aug-96  1.1d printf -> dprintf
 *	30-dec-97  1.1e ansi 
 *      20-jun-03  1.2  using modern get_snap and put_snap
 *      13-mar-05  1.2a  fix writing the time, free unused
 *       5-aug-06  2.0  change default of zerocm=, added shift=
 *
 */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>


string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in1=???\n			  input snapshot 1",
    "in2=???\n			  input snapshot 2",
    "out=???\n			  output snapshot",
    "deltar=0.0,0.0,0.0\n	  position of in1 w.r.t. in2",
    "deltav=0.0,0.0,0.0\n	  velocity of in1 w.r.t. in2",
    "zerocm=false\n		  zero the center of mass after stacking?",
    "shift=1\n                    Shift over 1st or 2nd one?",
    "headline=\n		  random verbiage for output",
    "VERSION=2.0\n		  5-aug-06 PJT",
    NULL,
};

string usage="stack two N-body systems on top of each other";

string cvsid="$Id$";

extern string *burststring(string,string);

nemo_main()
{
    readdata();
    snapstack();
    writedata();
    freedata();
}

local int nbody, nbody1, nbody2;
local int bits,  bits1,  bits2, shift;
local Body *btab, *btab1, *btab2;
local real tsnap, tsnap1, tsnap2;

readdata()
{
    stream instr1, instr2;

    instr1 = stropen(getparam("in1"), "r");
    get_history(instr1);
    instr2 = stropen(getparam("in2"), "r");
    get_history(instr2);

    shift = getiparam("shift");

    get_snap(instr1, &btab1, &nbody1, &tsnap1, &bits1);
    get_snap(instr2, &btab2, &nbody2, &tsnap2, &bits2);

    dprintf(1,"nbody1 = %d    nbody2 = %d\n", nbody1, nbody2);
    dprintf(1,"tsnap1 = %g    tsnap2 = %g\n", tsnap1, tsnap2);
}

snapstack()
{
    vector deltar, deltav;
    Body *bp;

    setvect(deltar, getparam("deltar"));
    setvect(deltav, getparam("deltav"));

    nbody = nbody1 + nbody2;
    btab1 = (body *) reallocate(btab1, sizeof(Body)*nbody);
    btab = btab1;
    tsnap = tsnap1;
    if (tsnap2 != tsnap)
      warning("tsnap2=%g not used, since tsnap1=%g",tsnap2,tsnap1);
    bits = (bits1 & bits2);
    memcpy(&btab[nbody1],btab2,sizeof(Body)*nbody2);

    if (shift==1) {      /* shift over 1st one */
      dprintf(1,"Shifting over 1st one\n");
      for (bp=btab; bp<btab+nbody1; bp++) {
	ADDV(Pos(bp),Pos(bp),deltar);
	ADDV(Vel(bp),Vel(bp),deltav);
      }
    } else if (shift==2) {      /* shift over 2nd one */
      dprintf(1,"Shifting over 2nd one\n");
      for (bp=&btab[nbody1]; bp<btab+nbody; bp++) {
	ADDV(Pos(bp),Pos(bp),deltar);
	ADDV(Vel(bp),Vel(bp),deltav);
      }
    }
    if (getbparam("zerocm")) snapcenter();
    else warning("zerocm=false is now the default!!");
}

setvect(vector vec, string str)
{
    string *vcp;
    int i;

    vcp = burststring(str, ", ");
    for (i = 0; i < NDIM; i++)
	vec[i] = (*vcp != NULL ? atof(*vcp++) : 0.0);
}

writedata()
{
    stream outstr;

    if (! streq(getparam("headline"), ""))
	set_headline(getparam("headline"));
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
}

freedata()
{
  free(btab1);
  free(btab2);
}

snapcenter() 
{
    real w_i, w_sum;
    vector tmpv, w_pos, w_vel;
    Body *bp;

    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    for (bp = btab; bp < btab+nbody;  bp++) {
      w_i = Mass(bp);
      w_sum += w_i;
      MULVS(tmpv, Pos(bp), w_i);
      ADDV(w_pos, w_pos, tmpv);
      MULVS(tmpv, Vel(bp), w_i);
      ADDV(w_vel, w_vel, tmpv);
    }
    SDIVVS(w_pos, w_sum);
    SDIVVS(w_vel, w_sum);
 
    for (bp = btab; bp < btab+nbody; bp++) {
        SSUBV(Pos(bp), w_pos);
        SSUBV(Vel(bp), w_vel);
    }
}
