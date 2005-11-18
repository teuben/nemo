/*
 * SNAPROTATE.C: read a snapshot file and rotate particle coordinates.
 *
 *	3-oct-88:	adapted to allow eulerian angles (order=zyz)	PJT
 *	17-feb-89:	new get+snap macros redef's			PJT
 *	27-nov-89	cosmetic defv[]                                 PJT
 *      14-nov-90       theta= now one keyword for theta1..theta3       PJT
 *	11-jun-92   V4.0 rotation is now mathematically positive, or
 *			 counterclockwise				PJT
 *	23-sep-95   V4.0a usage, RAD2DEG more accurate			pjt
 *	10-dec-95   V4.0b theta can be more than NDIM numbers		pjt
 *      14-may-96       c fixed frequent histout history                pjt
 *      16-feb-97       d SINGLEPREC fix                                pjt
 *	 8-jul-98	e fixed copying other input bodyparts		pjt
 *      21-nov-98   V5.0  Added tscale= parameter to scale with time    Pjt
 *      24-dec-04    5.0a fixed problem with uninitialized tscale       pjt
 *      18-nov-05    5.1  added option to select which vectors to rot   pjt
 *                        for Rachel's trials
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/barebody.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {	
    "in=???\n	           Input snapshot file",
    "out=???\n             Output snapshot file",
    "theta=\n              Angles (degrees) to rotate around axes",
    "order=\n              Order to do rotations (zyz=eulerian)",
    "invert=false\n        Invert transformation?",
    "tscale=\n             If used, scalefactor for angles with time of snapshot",
    "select=pos,vel,acc\n  Select the vectors for rotation (acc not implemented)",  
    "VERSION=5.1\n         18-nov-05 PJT",
    NULL,
};

string usage = "rotate a snapshot";

#define DEG2RAD PI/180.0



#define MAXANG  32

#define M_POS (1<<0)
#define M_VEL (1<<1)
#define M_ACC (1<<2)

void rotate(char axis,real angle,Body *btab,int nbody,bool Qscale,real scale, int vecmask);
void xrotate(Body *btab, int nbody, real theta, int vecmask);
void yrotate(Body *btab, int nbody, real theta, int vecmask);
void zrotate(Body *btab, int nbody, real theta, int vecmask);
void rotatevec(vector vec, matrix mat);

nemo_main()
{
    stream instr, outstr;
    char *op, *order;
    Body *btab = NULL;
    int i, nop, nopt, nbody, bits;
    real tsnap, tscale, *ang, theta[MAXANG];
    bool   invert, need_hist = TRUE;
    bool   Qtscale = hasvalue("tscale");
    string rotvects = getparam("select");
    bool Qpos, Qvel, Qacc;
    int vecmask = 0;
    
    if (strstr(rotvects,"pos")) vecmask |= M_POS;
    if (strstr(rotvects,"vel")) vecmask |= M_VEL;
    if (strstr(rotvects,"acc")) vecmask |= M_ACC;

    dprintf(0,"Warning: new option: select=%s  => mask=%d\n",rotvects,vecmask);

    instr = stropen(getparam("in"), "r");           /* open input file */
    nopt = nemoinpr(getparam("theta"),theta,MAXANG);  /* get angles */
    if (nopt<0)
        error("Parsing theta=%s returns %d",getparam("theta"),nopt);
    order = getparam("order");
    nop = strlen(order);        /* number of rotations to be applied */
    if (nop<=0) error("No order= specified");
    if (nop>MAXANG) {
        warning("Too many rotations specified in order=, only %d processed",
                                MAXANG);
        nop=MAXANG;
    }
    if (nopt>nop)
        error("Did not specify enough rotation axes: order=");
    else if (nopt<nop)
        error("Did not specify enough rotation angles: theta=");
    tscale = Qtscale ? getdparam("tscale") : 1.0;

    invert = getbparam("invert");
    outstr = stropen(getparam("out"), "w");

    get_history(instr);
    while (get_tag_ok(instr, SnapShotTag)) {
	get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (bits & PhaseSpaceBit) {
            dprintf(1,"Processing time %g bits=0x%x\n",tsnap,bits);
	    if (! invert) {
            	op = order;         /* reset order pointer */
            	ang = theta;        /* reset angle pointer */
                for (i=0; i<nop; i++) {
                    dprintf(2,"Rotating %g around %c\n",*ang,*op);
                    rotate(*op++, *ang++, btab, nbody, Qtscale, tscale*tsnap, vecmask);  /* VALGRIND error */
                }
	    } else {
                op = &order[nop-1];	/* reset: start at end */
                ang = &theta[nop-1];
                for (i=0; i<nop; i++) {
                    dprintf(2,"Rotating %g around %c\n",-*ang,*op);
		    rotate(*op--, -1.0*(*ang--), btab, nbody, Qtscale, tscale*tsnap, vecmask);
                }
	    }
	    if (need_hist) {
	        put_history(outstr);
	        need_hist = FALSE;
	    }
	    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	}
    }
    strclose(instr);
    strclose(outstr);
}

void rotate(char axis,real angle,Body *btab,int nbody,bool Qscale,real scale, int vecmask)
{
    if (Qscale) angle *= scale;

    dprintf(1,"angle=%g\n",angle);

    switch (axis) {
      case 'x':
      case 'X':
	xrotate(btab, nbody, angle, vecmask);
	break;
      case 'y':
      case 'Y':
	yrotate(btab, nbody, angle, vecmask);
	break;
      case 'z':
      case 'Z':
	zrotate(btab, nbody, angle, vecmask);
	break;
      default:
	warning("unknown axis %c: no rotation\n", axis);
    }
}

void xrotate(Body *btab, int nbody, real theta, int vecmask)
{
    matrix rmat;
    int i;

    SETMI(rmat);
    rmat[1][1] =    rmat[2][2] = cos(DEG2RAD * theta);
    rmat[1][2] =  -(rmat[2][1] = sin(DEG2RAD * theta));	/* PJT */
/*    rmat[2][1] =  -(rmat[1][2] = sin(DEG2RAD * theta));	/* JEB */
    for (i = 0; i < nbody; i++) {
      if (vecmask&M_POS) rotatevec(Pos(&btab[i]), rmat);
      if (vecmask&M_VEL) rotatevec(Vel(&btab[i]), rmat);
    }
}

void yrotate(Body *btab, int nbody, real theta, int vecmask)
{
    matrix rmat;
    int i;

    SETMI(rmat);
    rmat[2][2] =    rmat[0][0] = cos(DEG2RAD * theta);
    rmat[2][0] =  -(rmat[0][2] = sin(DEG2RAD * theta));	/* PJT */
/*    rmat[0][2] =  -(rmat[2][0] = sin(DEG2RAD * theta));		/* JEB */
    for (i = 0; i < nbody; i++) {
        if (vecmask&M_POS) rotatevec(Pos(&btab[i]), rmat);
	if (vecmask&M_VEL) rotatevec(Vel(&btab[i]), rmat);
    }
}

void zrotate(Body *btab, int nbody, real theta, int vecmask)
{
    matrix rmat;
    int i;

    SETMI(rmat);
    rmat[0][0] =    rmat[1][1] = cos(DEG2RAD * theta);
    rmat[0][1] =  -(rmat[1][0] = sin(DEG2RAD * theta));		/* PJT */
/*    rmat[1][0] =  -(rmat[0][1] = sin(DEG2RAD * theta));		/* JEB */
    for (i = 0; i < nbody; i++) {
        if (vecmask&M_POS) rotatevec(Pos(&btab[i]), rmat);
	if (vecmask&M_VEL) rotatevec(Vel(&btab[i]), rmat);
    }
}

void rotatevec(vector vec, matrix mat)
{
    vector tmp;

    MULMV(tmp, mat, vec);
    SETV(vec, tmp);
}
