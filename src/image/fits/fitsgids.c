/*
 *   FITSGIDS:   dump fits file onto the GIDS server
 *
 *      14-may-93   written  (not completed yet)
 *	20-jul-94   implemented the *range= keywords: recording planes
 *      21-jul-94   cashing in the look_ routines for efficiency
 *		    (improved a 512*512 image loading time from 1'43" to 3" !)
 *                  renamed xrange= -> x= etc. and id->record->plane=...    PJT
 *	13-apr-95   add filter (transfer function) before output
 *       1-feb-98   also display the 'gid' error message from xlibc's   PJT
 *                  (see look_nemo.c)
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <image.h>
#include <fitsio.h>

string defv[] = {
    "in=???\n		Input fits file",
    "datamin=\n         Datamin, if to override",
    "datamax=\n         Datamax, if to override",
    "x=\n               Range in X to display [all]",
    "y=\n               Range in Y to display [all]",
    "z=\n               Range in Y to display [first]",
    "plane=1\n          Plane to start recording in (0=no record)",
    "filter=linear\n    Filter applied to output",
    "VERSION=1.1a\n	1-feb-98 PJT",
    NULL,
};

string usage = "display fits images on GIDS display server";

#ifndef MAXSIZE
#define MAXSIZE  2048	/* max size in each dimension of an image */
#endif

extern string *burststring(string, string);

local int make_idx(string, int *, int, int);
local real filter_lin(double), filter_sqrt(double), filter_sqr(double),
           filter_log(double), filter_exp(double);
local int  match1(string, string);
local int  isblank(FLOAT *);

void nemo_main()
{
    FITS *fitsfile;
    string fitsname, filter;
    int ndim=3, naxis[3], nx, ny, nz, i, j, k, p,
	idx[MAXSIZE], idy[MAXSIZE], idz[MAXSIZE];
    int nbval=0, plane, bitpix, nblank;
    real bval_out, rmin, rmax, tmp;
    FLOAT *bufin;
    float datamin, datamax, datarange, bscale, bzero, *bufout;
    bool Qshort;
    rproc fie;
    string valid_filters = "linear,sqrt,sqr,log,exp";

    fitsname = getparam("in");
    fitsfile = fitopen(fitsname,"old",ndim,naxis);
    if (fitsfile==NULL) error("Could not open FITS file %s",getparam("in"));
    dprintf(0,"%s: %d * %d * %d\n",
    		getparam("in"),naxis[0],naxis[1],naxis[2]);

    plane = getiparam("plane");
    filter = getparam("filter");
    switch (match1(filter,valid_filters)) {
    case 0:
        fie = filter_lin;       break;
    case 1:
        fie = filter_sqrt;      break;
    case 2:
        fie = filter_sqr;       break;
    case 3:
        fie = filter_log;       break;
    case 4:
        fie = filter_exp;       break;
    default:
        error("Bad filter=%s; valid: %s",filter,valid_filters);
    }
    
    nx = make_idx(getparam("x"),idx,MAXSIZE,naxis[0]);
    ny = make_idx(getparam("y"),idy,MAXSIZE,naxis[1]);
    nz = make_idx(getparam("z"),idz,MAXSIZE,naxis[2]);

    if (nx < 1 || ny < 1 || nz < 1 )
	error("Cannot deal with (%d,%d,%d)",nx,ny,nz);

    if (hasvalue("datamin"))
        datamin=getdparam("datamin");
    else
        if (fitexhd(fitsfile,"DATAMIN"))
            fitrdhdr(fitsfile,"DATAMIN",&datamin,0.0); 
        else
            datamin=getdparam("datamin");
    if (hasvalue("datamax"))
        datamax=getdparam("datamax");
    else
        if (fitexhd(fitsfile,"DATAMAX"))
            fitrdhdr(fitsfile,"DATAMAX",&datamax,0.0); 
        else
            datamax=getdparam("datamax");
#if 0
    datamin = (*fie)((double)datamin);
    datamax = (*fie)((double)datamax);
#endif
    dprintf(0,"Displaying (%d,%d,%d) with %s filter datamin=%g datamax=%g\n",
            nx, ny, nz, filter, datamin, datamax);
    if (datamin >= datamax) error("Cannot display with datamin/max=%g %g",
    	    datamin, datamax);

    /* buffers to read and write */
    bufin  = (FLOAT *) allocate(naxis[0]*sizeof(FLOAT));
    bufout = (float *) allocate(nx*sizeof(float));
    rmin = HUGE;
    rmax = -HUGE;
    nblank = 0;
    look_init(fitsname,nx,ny,datamin,datamax,plane);
    for (k=0; k<nz; k++) {     
        p = idz[k]-1;
        fitsetpl(fitsfile,1,&p);
        for (j=0; j<ny; j++) {      /* loop over all rows */
            fitread(fitsfile,idy[j]-1,bufin);     /* read it from fits file */
            for (i=0; i<nx; i++) {                /* stuff it in memory */
                if (isblank(&bufin[idx[i]-1])) {
                    bufout[i] = bufin[i];
                    nblank++;
                } else {
#if 1
                    if (bufin[idx[i]-1] < datamin) bufin[idx[i]-1] = datamin;
                    if (bufin[idx[i]-1] > datamax) bufin[idx[i]-1] = datamax;
#endif
                    bufout[i] = (*fie)((double)bufin[idx[i]-1]);
                    rmin=MIN(rmin,bufout[i]);
                    rmax=MAX(rmax,bufout[i]);
                }
            }
            nbval += look_loop(j, bufout, nx);
        }
        look_flush();
        if (plane>0) {
	    if (look_record(plane)) {
                dprintf(0,"[Read plane %d and recorded in %d\n",p+1,plane);
                plane++;            
            } else
                dprintf(0,"[Read plane %d but not recorded\n",p+1);
	} else
            dprintf(0,"[Read plane %d]\n",p+1);
    }
    fitclose(fitsfile);

    dprintf(0,"Clipped %d pixels; Found datamin=%g datamax=%g Last REC=%d\n",
            nbval,rmin,rmax,plane-1);
    look_finis();
}

local int make_idx(string axis, int *idx, int maxlen, int axlen)
{
    int i, j, n;
    if (axlen>maxlen)
        warning("Axis %s is too long (%d) to hold %d pixels",axis,axlen,maxlen);
    n = nemoinpi(axis,idx,maxlen);
    if (n>0) {
    	for (i=0, j=0; i<n; i++) {
            if (idx[i]<1 || idx[i]>axlen) continue; /* don't copy: */
            idx[j++] = idx[i];
    	}
    	if (j!=n)
            warning("Only %d/%d elements out of range %s used",j,n,axis);
        return j;
    }
    if (n<0) error("Parsing error %s",axis);
    for (i=0; i<axlen; i++) idx[i] = i+1;
    return axlen;
}

local int match1(string name, string options)
{
    string *sp;
    int i, nsp, retval=-1;
    int len = strlen(name);
    
    sp = burststring(options,",");
    nsp = xstrlen(sp,sizeof(string))-1;
    for (i=0; i<nsp; i++) {
        if (strncmp(sp[i],name,len) == 0) {
            if (retval>=0) {
                warning("%s: matches %s and %s",name,sp[retval],sp[i]);
                return -1;
            } else
                retval = i;
        }
    }
    freestrings(sp);
    return retval;
}

local int blank_val = 0xFFFF;		/* yuck */

local int isblank(FLOAT *bufin)
{
    return memcmp(bufin,&blank_val,sizeof(int))==0;    
}

real filter_lin(double x)  { return x; }
real filter_sqrt(double x) { return sqrt(x); }
real filter_sqr(double x)  { return x*8; }
real filter_exp(double x)  { return exp(x); }
real filter_log(double x)  { return log(x); }


