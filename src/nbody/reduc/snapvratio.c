/*
 * SNAPVRATIO.C: computes (clausius) virial from present forces and potential
 *
 *	05-Nov-91  V0.0  calculate virial radius of a cluster	SMF
 *	26-mar-92  V0.1  added clausius	and full <body.h>	PJT
 *	 3-arp-92  V0.2  added keyword mode=			PJT
 *       6-apr-92  V0.3  an extra check if acc or phi present   PJT
 *	21-nov-96  V0.3a ** total mass=1 assumed ???            PJT
 *	13-mar-97  V0.4  extra column w/ mass			pjt
 *	30-jul-97  V0.5  added eps= for softening		pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>


string defv[] = {
    "in=???\n         Input file name (snapshot)",
    "times=all\n      Times to select snapshot",
    "wmode=acc\n      Use for W (acc|phi|exact)",
    "newton=f\n       Do an exact N^2 newtonian calculation too?",
    "eps=0.05\n	      Gravitational softening, if used",
    "VERSION=0.5\n    30-jul-97 PJT",
    NULL,
};

string usage = "compute various global virials (clausius, newton)";

#define ACC_MODE   0x01
#define PHI_MODE   0x02
#define EXACT_MODE 0x04


nemo_main()
{
    stream instr;
    real   tsnap,T2,v2,rij,s, eps, eps2;
    real   w_acc, w_phi, w_exact, tmass;
    vector tmp;
    int    i, nbody, bits, wmode, m, count, match();
    Body   *btab = NULL, *bp1,*bp2;
    bool   Qnewton=getbparam("newton");

    if ((m=match(getparam("wmode"),"acc phi exact",&wmode)) != 1)
        error("match=%d Bad wmode=%s\n",m,getparam("wmode"));

    eps = getdparam("eps");
    eps2 = sqr(eps);

    dprintf(1,"wmode=0x%x\n",wmode);
    printf("# T  2T/W  T+W   T  W_acc W_phi W_exact\n");
    
    instr = stropen(getparam("in"), "r");           /* open input file */
    get_history(instr);			    /* accumulate data history */
    count = 0;
    for(;;) {				 	 /* loop for all times */
        get_history(instr);                         /* for paranoidici */
        if (!get_tag_ok(instr, SnapShotTag))          /* check if done */
	    break;
        get_snap(instr, &btab, &nbody, &tsnap, &bits);      /* get one */
        if ((bits & PhaseSpaceBit) == 0) continue;
        if (wmode&ACC_MODE && (bits&AccelerationBit)==0) continue;
        if (wmode&PHI_MODE && (bits&PotentialBit)==0) continue;

        count++;
	T2 = w_acc = w_phi = w_exact = tmass = 0.0;
	for(bp1=btab;bp1<btab+nbody;bp1++) {
	    tmass += Mass(bp1);
            DOTVP(v2,Vel(bp1),Vel(bp1));
            T2 += Mass(bp1)*v2;
	    DOTVP(s,Pos(bp1),Acc(bp1));
	    w_acc += Mass(bp1)*s;
	    w_phi += Mass(bp1)*Phi(bp1);
	    if (Qnewton) 
		for(bp2=bp1+1;bp2<btab+nbody;bp2++) {
	           SUBV(tmp,Pos(bp1),Pos(bp2));
		   ABSV(rij,tmp);
		   if (eps2 > 0.0)
    		      w_exact -= Mass(bp1)*Mass(bp2)/sqrt(rij*rij+eps2);
    		   else
    		      w_exact -= Mass(bp1)*Mass(bp2)/rij;
	        }
	}
	w_phi *= 0.5;
        if (wmode==ACC_MODE)
	    printf("%f %f %f %f %f %f %f %f\n",
		tsnap,-T2/w_acc,T2/2+w_acc,T2/2,w_acc,w_phi,w_exact,tmass);
        else if (wmode==PHI_MODE)
	    printf("%f %f %f %f %f %f %f %f\n",
		tsnap,-T2/w_phi,T2/2+w_phi,T2/2,w_phi,w_acc,w_exact,tmass);
        else if (wmode==EXACT_MODE)
	    printf("%f %f %f %f %f %f %f %f\n",
		tsnap,-T2/w_exact,T2/2+w_exact,T2/2,w_exact,w_acc,w_phi,tmass);

#if 1
	free(btab);
	btab = NULL;
#endif		
    }   /* for(;;) */

    if (count==0) warning("No work done, use another wmode?");
} /* nemo_main() */
