/*
 * SNAPXYZ: convert snapshots to xyzc data for viewing.
 *	V1.2    1-may-90 Created by Josh Barnes  
 *	V2.0   21-jan-93 Based on bodytrans()			PJT
 *	V2.0a  29-mar-94 ansi
 *	V2.1   30-mar-97 exported in nbody/xyz, added vel=t|f
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <vectmath.h>
#include <loadobj.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

string defv[] = {	
    "in=???\n			Input file of snapshot data",
    "out=???\n			Output file of xyzc data",
    "times=all\n		Range of times to convert",
    "xvar=x\n			Expression for x-axis variable",
    "yvar=y\n			Expression for y-axis variable",
    "zvar=z\n			Expression for z-axis variable",
#ifdef XYZVEL
    "vxvar=vx\n			Expression for vx-axis variable",
    "vyvar=vy\n			Expression for vy-axis variable",
    "vzvar=vz\n			Expression for vz-axis variable",
    "vel=t\n                    Should velocitied be output at all?",
#endif
    "color=1\n			Expression for point color",
    "visib=1\n			Expression for point visibility",
    "VERSION=2.1b\n		30-may-04 PJT",
    NULL,
};
string usage = 	"Convert snapshot to xyzc data";

local rproc xfunc, yfunc, zfunc;
#ifdef XYZVEL
local rproc vxfunc, vyfunc, vzfunc;
local bool  Qvelout;
#endif

local iproc vfunc, cfunc;

local string times;

void setparams(void);
bool get_frame(stream instr);
void put_points(stream outstr);
void convert(void);

extern rproc btrtrans();    /* ??? */
extern iproc btitrans();

nemo_main()
{
    stream instr, outstr;

    setparams();
    
    instr = stropen(getparam("in"), "r");
    get_history(instr);

    outstr = stropen(getparam("out"), "w");
    put_history(outstr);

    while (get_frame(instr)) {
	convert();
	put_points(outstr);
    }
}

void setparams(void)
{
    xfunc = btrtrans(getparam("xvar"));
    yfunc = btrtrans(getparam("yvar"));
    zfunc = btrtrans(getparam("zvar"));
#ifdef XYZVEL    
    vxfunc = btrtrans(getparam("vxvar"));
    vyfunc = btrtrans(getparam("vyvar"));
    vzfunc = btrtrans(getparam("vzvar"));
    Qvelout = getbparam("vel");
#endif    

    vfunc = btitrans(getparam("visib"));
    cfunc = btitrans(getparam("color"));

    times = getparam("times");
}

local Body *btab = NULL;
local int nbody, nmax = 0;
local real tsnap = 0.0;

bool get_frame(stream instr)
{
    int bits;

    do {
	if (nmax > 0) nbody = nmax;
	get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
	nmax = MAX(nmax, nbody);
	if (bits & PhaseSpaceBit) return TRUE;
    } while (bits != 0);
    return FALSE;
}

local vector *ptab = NULL;
local vector *vtab = NULL;
local int *ctab = NULL;

local int npoint;
local bool cnew;

void put_points(stream outstr)
{
    put_set(outstr, "PointData");
    put_data(outstr, "Tpoint", RealType, &tsnap, 0);
    put_data(outstr, "Npoint", IntType, &npoint, 0);
    put_data(outstr, "Position", RealType, ptab, npoint, 3, 0);
#ifdef XYZVEL
    if (Qvelout)
        put_data(outstr, "Velocity", RealType, vtab, npoint, 3, 0);
#endif
    if (cnew)
	put_data(outstr, "Color", IntType, ctab, npoint, 0);
    put_tes(outstr, "PointData");
}


void convert(void)
{
    int visnow, vismax, new_color, i, ip, vis;
    Body b;

    if (ptab == NULL) ptab = (vector *) allocate(nbody * sizeof(vector));
#ifdef XYZVEL
    if (Qvelout)    
        if (vtab == NULL) vtab = (vector *) allocate(nbody * sizeof(vector));
#endif
    if (ctab == NULL) {
	ctab = (int *) allocate(nbody * sizeof(int));
	for (i = 0; i < nbody; i++) ctab[i] = 1;
	cnew = TRUE;
    } else
	cnew = FALSE;
    npoint = 0;
    visnow = vismax = 0;
    do {                                    /* loop over various layers */
	visnow++;                           /* activate next layer */
	for (i = 0; i < nbody; i++) {       /* loop over all bodies */
            b = btab[i];
	    vis = (*vfunc)(&b, tsnap, i);
	    vismax = MAX(vismax, vis);      /* remember highest value */
	    if (vis == visnow) {            /* if body is visible */
		ptab[npoint][0] = (*xfunc)(&b,tsnap,i);     /* get xyz */
		ptab[npoint][1] = (*yfunc)(&b,tsnap,i);
		ptab[npoint][2] = (*zfunc)(&b,tsnap,i);
#ifdef XYZVEL
                if (Qvelout) {
		   vtab[npoint][0] = (*vxfunc)(&b,tsnap,i);    /* get vel's */
		   vtab[npoint][1] = (*vyfunc)(&b,tsnap,i);
		   vtab[npoint][2] = (*vzfunc)(&b,tsnap,i);
                }
#endif
		new_color = (*cfunc)(&b,tsnap,i);           /* get color */
		if (ctab[npoint] != new_color) {
		    ctab[npoint] = new_color;
		    cnew = TRUE;
		}
		npoint++;
	    }
	}
    } while (visnow < vismax);
    dprintf(0, "[%s: npoint = %d  time = %.2f]\n",
	    getargv0(), npoint, tsnap);
}


