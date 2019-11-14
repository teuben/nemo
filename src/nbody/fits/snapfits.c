/* 
 *  SNAPFITS:  write a snapshot into fits format, 
 *              in random groups format that is
 *             
 *      19-mar-90  V1.0 -- created as a toy model     PJT
 *			Must be changed to 3DTABLE to be standard FITS
 *       9-aug-04  1.0a - bring code to 2000+
 */

#include <stdinc.h>             /* also gets <stdio.h>  */
#include <getparam.h>
#include <vectmath.h>           /* otherwise NDIM undefined */
#include <history.h>
#include <filestruct.h>

#include <snapshot/body.h>      /* snapshot's */
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

#include <fits.h>

string defv[] = {
        "in=???\n               Input filename (a snapshot)",
        "out=???\n              Output filename (a fits file)",
        "options=mass,phase\n   Options for output",
        "comment=\n             Additional comments for FITS file",
        "VERSION=1.0a\n         9-aug-04 PJT",            
        NULL,
};

string usage = "write a snapshot into fits format";

static stream  instr, outstr;                          /* file streams */

                /* SNAPSHOT INTERFACE */
static int    nbody, bits;
static real   tnow;
static Body   *btab = NULL;
                /* FITS interface */
struct fits_header fh;
static char *comments[] = {NULL, NULL};        /* place to save ONE comment */

#define MAXNEEDS 10
static int needs[MAXNEEDS];     /* 0=mass, 1=pos, 2=vel 3=phi, 4=acc */


nemo_main(void)
{
    instr = stropen (getparam ("in"), "r");     /* open snapshot */
    get_history(instr);
    get_snap(instr,&btab,&nbody,&tnow,&bits);
    strclose(instr);

    outstr = stropen(getparam("out"), "w");
    w_header(outstr,&fh,nbody,getparam("options"),getparam("comment"));
    w_data(outstr,&fh,nbody,btab);
}

need_pars(string options)
{
    int i, ntot;
    bool scanopt();
    for (i=0; i<MAXNEEDS; i++)
	needs[i] = 0;
    needs[0] = (scanopt(options,"mass") ? 1 : 0);
    if (scanopt(options,"phase"))
        needs[1] = needs[2] = NDIM;
    else
        needs[1] = needs[2] = 0;
    if (scanopt(options,"pos"))         /* catch pos,vel <-> phase */
        needs[1] = NDIM;
    if (scanopt(options,"vel"))
        needs[2] = NDIM;
    needs[3] = (scanopt(options,"phi") ? 1 : 0);
    needs[4] = (scanopt(options,"acc") ? NDIM : 0);

    for (ntot=0, i=0; i<MAXNEEDS; i++)
        ntot += needs[i];

    if (ntot==0)
        error("options musyt contain at least one of: mass,phase,pos,vel\n");
    else
        dprintf(0,"Options=%s outputs %d items per body\n",options,ntot);

    return(ntot);
}

w_header(stream ostr,struct fits_header *fh, int nbody, char *options,char *comment)
{   
    static int naxis[1] = { 0 };    /* random groups ! */
    int i, n, ni;
    string *ptype, copxstr();
    char pname[9];

    fts_zero(fh);          /* clear primary header */

    fh->simple = 1;          /* set various primary header stuff */
    fh->groups = 1;
    fh->bitpix = -64;        /* IEEE double precision FOR NOW */
    fh->naxis = 1;
    fh->naxisn = naxis;      /* point to the array */
    fh->pcount = n = need_pars(options);
    fh->gcount = nbody;
    fh->history = ask_history();    /* fill in some history */
    if (*comment != 0) {         /* were there any comments ? */
        comments[0] = comment;        
        fh->comment = comments;
    }
    
    ptype = (char **) malloc(n*sizeof(char *));
    ni = 0;
    if (needs[0]>0)
        ptype[ni++] = copxstr("MASS",sizeof(char));
    if (needs[1]>0)
        for (i=0; i<NDIM; i++) {
            sprintf(pname,"POS%1d",i+1);
            ptype[ni++] = copxstr(pname,sizeof(char));
        }
    if (needs[2]>0)
        for (i=0; i<NDIM; i++) {
            sprintf(pname,"VEL%1d",i+1);
            ptype[ni++] = copxstr(pname,sizeof(char));
        }
    if (needs[3]>0)
        ptype[ni++] = copxstr("PHI",sizeof(char));
    if (needs[4]>0)
        for (i=0; i<NDIM; i++) {
            sprintf(pname,"ACC%1d",i+1);
            ptype[ni++] = copxstr(pname,sizeof(char));
        }
    for (ni=0; ni<n; ni++)
        dprintf(0,"PTYPE%d = %s\n",ni+1,ptype[ni]);
    fh->ptypen = ptype;
    
    
    fts_whead(fh,ostr);      /* write the header */
}

w_data(stream ostr,struct fits_header *fh,int nbody,Body *btab)
{
    int i, nw, need_mass,need_pos,need_vel, ntot, isize;
    Body *bp;
    double out[32];   /* scratch array for buffering up a body's parameters */
    char null = 0;
    
    ntot = 0;
    isize = sizeof(double);
    for (bp=btab;bp<btab+nbody;bp++) {
        nw = 0;
        if (needs[0]>0)
            out[nw++] = Mass(bp);
        if (needs[1]>0)
            for (i=0; i<NDIM; i++)
                out[nw++] = Pos(bp)[i];
        if (needs[2]>0)
            for (i=0; i<NDIM; i++)
                out[nw++] = Vel(bp)[i];
        if (needs[3]>0)
            out[nw++] = Phi(bp);
        if (needs[4]>0)
            for (i=0; i<NDIM; i++)
                out[nw++] = Acc(bp)[i];
            
        fts_wdata(fh,ostr,nw*isize,(char *)out);
        ntot += nw*isize;
    }
    /* make sure tail end filled with zero's */
    i = FTSBLKSIZ - (ntot % FTSBLKSIZ);
    if (i==FTSBLKSIZ) i=0;
    dprintf(0,"Wrote %d bytes, flushing %d zero's\n",ntot,i);
    while (i-- > 0)
       fts_wdata(fh,ostr,1,&null);
    
}
