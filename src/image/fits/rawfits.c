/*
 *   RAWFITS:   raw image formats to fits
 *
 *	11-feb-99  V1.0	Created for A310 to replace obscure J routine   PJT
 *                      does not pad the data yet ... 
 *			Usage: rawfits $f.img $f.fits 16 384,576 160
 *      12-feb-99  1.1  Added swap=
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input raw image bytes",
    "out=???\n	           Output fits file",
    "bitpix=8\n            Bitpix of the raw image",
    "naxis=\n              Size of image (naxis1,naxis2)",
    "skip=0\n              Number of bytes to skip before data",
    "swap=f\n              Additional swap on data?",
    "VERSION=1.1\n         12-feb-99 PJT",
    NULL,
};

string usage = "raw image to fits format conversion";

#include "fits.h"

local string minmax[] = { "DATAMIN", "DATAMAX", NULL };

void nemo_main()
{
    stream instr, outstr;
    int    i,j, n, blocking[2], planesize, totsize;
    int    bitpix = getiparam("bitpix");
    int    zbufsiz, ninfiles, naxis, naxis1, naxis2, naxis3, naxisn[3], ntmp;
    struct fits_header fh, fh_out;
    string outname, inname;
    char   *buffer, *zbuf, fname[128];
    int    nret, skip = getiparam("skip");
    double datamin, datamax;
    bool   Qswap = getbparam("swap");

    nret = nemoinpi(getparam("naxis"),naxisn,3);
    if (nret==1) {
        naxisn[1] = naxisn[0];
        naxisn[2] = 1;
    } else if (nret == 2) {
	naxisn[2] = 1;
    } else if (nret != 3)
	error("Need 1,2 or 3 values; Parsing naxis=%s",getparam("naxis"));

    inname = getparam("in");
    outname = getparam("out"); 

    outstr = stropen(outname,"w");
    instr = stropen(inname,"r");

    
    fts_zero(&fh);                                    /* make blank header */
    if (Qswap) fh.flip = 1;
    fts_zero(&fh_out);
    n = fts_xhead(&fh,instr,bitpix,nret,naxisn,skip); /* fake reading header */
    if (n<0) error("xhead bad header");
    memcpy((char *)&fh_out,(char *)&fh,sizeof(fh));   /* copy header */
	/* should delete DATAMIN/DATAMAX */
    fts_whead(&fh_out,outstr);                        /* write header */
#if 1
	/* this appears to have a bug for rawfits; some part of
           header not well initialized */
    fts_cdata(&fh,instr,outstr,FALSE,TRUE);           /* copy data + pad */
#else    
    fts_cdata(&fh,instr,outstr,FALSE);                /* copy data */
#endif    
    strclose(instr);
    strclose(outstr);
}
