/*
 * SNAPINERT: obtain inertia tensor & its eigenvectors, eigenvalues
 * on some subset of particles. See also: snaprect
 *
 *    Updates:
 *	??	    V1.0  Created - Jun Makino ?
 *	7-jul-89    V1.2  updated for get_snap()	PJT
 *      9-dec-90    V1.3  helpvec - new numrec macros via <nrutil>  PJT
 *	1-jun-92    V1.4  usage
 *     12-apr-97       a  fixed bug (free_fmatrix -> free_convert_matrix)
 *     17-jul-02    V2.0  convert to use with the nemofied NumRec2 routines
 *     31-dec-02    V2.1  gcc3/SINGLEPREC 
 *	
 *  Note: this program only works for NDIM=3
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <bodytransc.h>

#include <nrutil.h>

string defv[] = {	
    "in=???\n			  Input snapshot file name",
    "out=???\n			  Output file name (ascii table)",
    "times=all\n		  range of times to transform",
    "weight=m\n		          factor to use in computing center/inertia  ",
    "per_weight=t\n		  flag to give inertia per weight basis  ",
    "tab=f\n			  flag to produce one-line output",
    "VERSION=2.1\n		  31-dec-02 PJT",
    NULL,
};

string usage="get inertia tensor & its eigenvectors, eigenvalues";

#define TIMEFUZZ 0.001          /* slop tolerated in time comparisons */

void snapinert(Body *, int, real, rproc_body, float i[3][3]);


nemo_main()
{
    stream instr, outstr;
    string times;
    bool oneline;
    rproc_body weight;
    Body *btab=NULL;
    int nbody, bits, nrot, i, j;
    real tsnap;
    float inert[3][3];    
    float inert0[3][3];
    float **convi, *eigens, **eigenvs;
    
    instr = stropen(getparam("in"), "r");
    times = getparam("times");
    oneline=getbparam("tab");
    weight = btrtrans(getparam("weight"));
    outstr = stropen(getparam("out"), "w");

    eigens=fvector(1,3);
    eigenvs=fmatrix(1,3,1,3);

    do {
        get_history(instr);
	get_snap(instr, &btab, &nbody, &tsnap, &bits);

	if ((bits & PhaseSpaceBit) != 0 &&
	    (streq(times, "all") ||
	     ((bits & TimeBit) != 0 && within(tsnap, times, TIMEFUZZ)))) {

	    snapinert(btab, nbody, tsnap, weight, inert);
	    for(i=0;i<3;i++)for(j=0;j<3;j++) inert0[i][j]=inert[i][j];
	    convi = convert_matrix(&inert[0][0],1,3,1,3);
	    jacobi(convi, 3, eigens, eigenvs, &nrot);
	    eigsrt(eigens, eigenvs, 3);
	    if(!oneline){
		fprintf(outstr, "Time: %f\n", tsnap);
		for(i=0; i<3; i++){
		    for(j=0; j<3; j++)fprintf(outstr, "  %12.6f", inert0[i][j]);
		    fprintf(outstr,"\n");
		}
		fprintf(outstr,"Eigenvalues:\n");
		for(i=1; i<=3; i++){
		    fprintf(outstr, "  %12.6f  :  ", eigens[i]);
		    for(j=1; j<=3; j++)fprintf(outstr, "  %12.6f", eigenvs[j][i]);
		    fprintf(outstr,"\n");
		}
	    }else{
		fprintf(outstr,"%12.6f  ", tsnap);

		for(i=0; i<3; i++){
		    for(j=0; j<3; j++){
			fprintf(outstr, "  %12.6f", inert0[i][j]);}}
		fprintf(outstr, "   ");
		for(i=1; i<=3; i++)fprintf(outstr, " %12.6f", eigens[i]);
		for(i=1; i<=3; i++){
		    for(j=1; j<=3; j++){
			fprintf(outstr, "  %12.6f", eigenvs[j][i]);}}
		fprintf(outstr,"\n");
	    }
            free_convert_matrix(convi,1,3,1,3);
	}
    } while (bits != 0);
}

void snapinert(
	       Body *btab,
	       int nbody,
	       real tsnap,
	       rproc_body weight,
	       float inert[3][3])
  {
    int i,j,k;
    Body *b;
    real w_i, w_sum;
    vector tmpv, w_pos, w_vel;
    bool per_weight;

    per_weight=getbparam("per_weight");
    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    for (i=0; i<3; i++)CLRV(inert[i]);
    for (i = 0, b = btab; i < nbody; i++, b++) {        /* get C.O.M. */
	w_i = (weight)(b, tsnap, i);
	if (w_i < 0.0)
	    warning("weight[%d] = %g < 0", i, w_i);
	w_sum += w_i;
	MULVS(tmpv, Pos(b), w_i);
	ADDV(w_pos, w_pos, tmpv);
	MULVS(tmpv, Vel(b), w_i);
	ADDV(w_vel, w_vel, tmpv);
    }
    if (w_sum == 0.0) error("total weight is zero");
    DIVVS(w_pos, w_pos, w_sum);
    DIVVS(w_vel, w_vel, w_sum);
    for (i = 0, b = btab; i < nbody; i++, b++) {        /* get inertia */
	SUBV(Pos(b), Pos(b), w_pos);
	SUBV(Vel(b), Vel(b), w_vel);
	w_i = (weight)(b, tsnap, i);
	for(k=0; k<3; k++) for(j=0; j<3; j++){
	    inert[k][j]+= Pos(b)[k]*Pos(b)[j]*w_i;
	}
    }
    if(per_weight) for(k=0; k<3; k++) for(j=0; j<3; j++){
	inert[k][j] = inert[k][j] /w_sum;
    }
}
