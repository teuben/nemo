/*
 * SNAPCENTER: transform snapshot file to coordinates centered
 * on some subset of particles.
 *	
 *      22-jan-89   something                   JEB
 *         mar-89   added output                PJT
 *	21-nov-90   helpvec			PJT
 *	28-apr-92   small cleanup 
 *	19-nov-92   1.3 write pos/vel into history  PJT
 *	 9-nov-93   1.4 optional non-reporting mode PJT
 *	 7-jan-96   1.5 optional output of COM system instead	PJT
 *	26-feb-97   1.6 made report=f the default               pjt
 *      31-dec-02   1.7 fixed for gccd3/SINGLEPREC              pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <history.h>
#include <filestruct.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <bodytransc.h>

string defv[] = {	
    "in=???\n       input file name ",
    "out=???\n	    output file name ",
    "weight=1\n	    factor used finding center ",
    "times=all\n    range of times to process ",
    "report=f\n	    report the c.o.m shift",
    "one=f\n        Only output COM as a snapshot?",
    "VERSION=1.7\n  31-dec-02 PJT",
    NULL,
};

string usage="Center a snapshot based on a weighed subset of particles";


void snapcenter(Body*, int, real, rproc_body, vector, vector, bool);

void nemo_main()
{
    stream instr, outstr;
    string times;
    rproc_body weight;
    Body *btab = NULL, *b;
    int i, nbody, bits;
    real tsnap, mass;
    bool Qreport, Qone;
    vector c_pos, c_vel;
    char line[256];

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");
    weight = btrtrans(getparam("weight"));
    times = getparam("times");
    Qreport = getbparam("report");
    Qone = getbparam("one");
    if (Qreport) dprintf(1,"pos vel of center(s) will be:\n");

    get_history(instr);
    put_history(outstr);

    do {
	get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
	if (bits & PhaseSpaceBit) {
	    snapcenter(btab, nbody, tsnap, weight, c_pos, c_vel, Qreport);
#if 0
            /* A better mechanism is needed for this */
            /* put_string() is not a robust in tsf/rsf conversion */
            sprintf(line,"snapcenter: time= %f pos= %f %f %f vel= %f %f %f\n",
                    tsnap, c_pos[0], c_pos[1], c_pos[2],
                           c_vel[0], c_vel[1], c_vel[2]);
            put_string(outstr,HistoryTag,line);
#endif
            if (Qone) {         /* one particle pos-vel output */
                bits = TimeBit | PhaseSpaceBit | MassBit;
                SETV(Pos(btab),c_pos);
                SETV(Vel(btab),c_vel);
                for (i = 0, b = btab; i < nbody; i++, b++) 
                    mass += Mass(b);
                Mass(btab) = mass;
                nbody = 1;
            }
	    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	}
    } while (bits != 0);
    strclose(outstr);
}

void snapcenter(
		Body *btab,
		int nbody,
		real tsnap,
		rproc_body weight,
		vector w_pos, 
		vector w_vel,
		bool Qreport)
{
    int i;
    Body *b;
    real w_i, w_sum;
    vector tmpv;

    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_i = (weight)(b, tsnap, i);
	if (w_i < 0.0)
	    warning("weight[%d] = %g < 0\n", i, w_i);
	w_sum += w_i;
	MULVS(tmpv, Pos(b), w_i);
	ADDV(w_pos, w_pos, tmpv);
	MULVS(tmpv, Vel(b), w_i);
	ADDV(w_vel, w_vel, tmpv);
    }
    if (w_sum == 0.0)
	error("total weight is zero");
    SDIVVS(w_pos, w_sum);
    SDIVVS(w_vel, w_sum);

    if (Qreport) {
      for (i=0; i<NDIM; i++)
        printf("%f ",w_pos[i]);
      for (i=0; i<NDIM; i++)
        printf("%f ",w_vel[i]);
      printf("\n");
    } else {
      for (i=0; i<NDIM; i++)
        dprintf(1,"%f ",w_pos[i]);
      for (i=0; i<NDIM; i++)
        dprintf(1,"%f ",w_vel[i]);
      dprintf(1,"\n");
    }

    for (i = 0, b = btab; i < nbody; i++, b++) {
	SSUBV(Pos(b), w_pos);
	SSUBV(Vel(b), w_vel);
    }
}
