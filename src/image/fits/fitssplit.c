/*
 *   FITSSPLIT:   split a fits dataset (often a cube) into smaller ones
 *		  (often planes)
 *
 *	14-aug-90       V1.0 created under NEMO		PJT
 *	 7-mar-92	V1.1 new nemo, gcc2 happy	PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input fits file",
    "out=???\n	           Basename for output fits file",
    "file=0\n              Select which planes(s)? [0=all]",
    "blocking=1,1\n        Blocking factors for I/O",
    "VERSION=1.2a\n        25-mar-94 PJT",
    NULL,
};

string usage = "split fits image dataset into subsets";

#include "fits.h"

void nemo_main()
{
    stream instr, outstr;
    int    i,n,nfile, blocking[2], planesize, totsize;
    int    zbufsiz;
    struct fits_header fh;
    char  *basename;
    char  *buffer, *zbuf, fname[128];

    instr = stropen(getparam("in"),"r");    /* open input */
    basename = getparam("out");          /* basename for output */

    nfile = getiparam("file");

    n = nemoinpi(getparam("blocking"),blocking,2);
    if (n<1 || (outstr!=NULL && n<2))
        error("Not enough blocking factors supplied");

    fts_setiblk(blocking[0]);    /* set input blocking factor */
    fts_setoblk(blocking[1]);    /* set input blocking factor */
	
    fts_zero(&fh);			             /* clean out header */
    n = fts_rhead(&fh,instr);	                     /* read header */
    if (n<0) {				              /* if no data (EOF) .. */
        warning("No data???");
        return ;
    }
    totsize = ABS(fh.bitpix) * fh.dlen / 8;
    planesize = ABS(fh.bitpix) * fh.naxisn[0] * fh.naxisn[1] / 8;
    if (totsize % planesize != 0) {
        dprintf(0,"Total size = %d bytes\n",totsize);
        dprintf(0,"Plane size = %d bytes\n",planesize);
        dprintf(0,"*** Not commensurate\n");
    }
    if (totsize/planesize != fh.naxisn[2]) {
        dprintf(0,"Number of size/  planes = %d \n",totsize/planesize);
        dprintf(0,"Number of naxis3 planes = %d \n",fh.naxisn[2]);
        dprintf(0,"*** Not the same - last one taken\n");
    }
    buffer = allocate(planesize);
    zbufsiz = planesize % 2880;
    if (zbufsiz>0) {
        zbufsiz = 2880 - zbufsiz;
        zbuf = allocate(zbufsiz);        
    }
    fh.naxis = 2;
    for (i=0; i<fh.naxisn[2]; i++) {
        sprintf(fname,"%s.%d",basename,i);
        dprintf(0,"[Working on %s]\n",fname);
        outstr = stropen(fname,"w");
        fts_whead(&fh,outstr);                            /* copy header */
        fread(buffer,1,planesize,instr);
        fwrite(buffer,1,planesize,outstr);
        if (zbufsiz>0) fwrite(zbuf,1,zbufsiz,outstr);    
        strclose(outstr);
    }
    strclose(instr);                    /* close files */
}
