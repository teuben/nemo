/*
 *   FITSGLUE:   catentate fits planes into a cube
 *
 *	9-aug-93  V1.0	Created				pjt
 *	7-oct-84  V1.0a fixed bug
 *     28-aug-86  V1.1  added compact=                  pjt
 *     24-feb-98   1.2  added inlist= to allow large in= lists  PJT
 *     25             a fixed datamin/max in the header         PJT
 *     12-feb-99      b changed for new fts_cdata 		PJT
 *	5-apr-01      c increased default MAXIN			pjt
 *
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input fits files or template for list",
    "out=???\n	           Output fits file",
    "multiple=t\n          Check for multiple fits files on input?",
    "naxis3=0\n            Preset NAXIS3 here, or let program scan first",
    "blocking=1,1\n        Blocking factors for I/O",
    "compact=f\n           Compact (move) dummy axes to the end",
    "inlist=\n             optional nemoinp(1) list expression for in=",
    "VERSION=1.2c\n        5-apr-01 PJT",
    NULL,
};

string usage = "catenate fits images into a fitscube";

#include "fits.h"

extern string *burststring(string,string);

#ifndef MAXIN
#define MAXIN 8096
#endif

local string minmax[] = { "DATAMIN", "DATAMAX", NULL };

void nemo_main()
{
    stream instr, outstr;
    int    i,j, n, blocking[2], planesize, totsize;
    int    zbufsiz, ninfiles, naxis, naxis1, naxis2, naxis3, naxisn[3], ntmp;
    bool   Qmultiple, Qcompact = getbparam("compact");
    struct fits_header fh, fh_out;
    string outname, *innames;
    char   *buffer, *zbuf, fname[128];
    int    inlist[MAXIN];
    double datamin, datamax;

    innames = burststring(getparam("in"),", ");
    ninfiles = xstrlen(innames,sizeof(string)) - 1;
    if (hasvalue("inlist")) {
        if (ninfiles!=1)
            error("inlist= requires one template in=%s\n",getparam("in"));
        if (strchr(getparam("in"),'%')==0)
            error("input template file needs printf-type formatting directive");
        ninfiles = nemoinpi(getparam("inlist"),inlist,MAXIN);
        if (ninfiles < 1)
            error("Need a list of numbers, or Syntax error inlist=");
        innames = (string *) allocate(ninfiles*sizeof(string));
        for (j=0; j<ninfiles; j++) {
            sprintf(fname,getparam("in"),inlist[j]);
            dprintf(1,"%d: Using file %s\n",fname);
            innames[j] = strdup(fname);
        }
    }
    outname = getparam("out"); 

    n = nemoinpi(getparam("blocking"),blocking,2);
    if (n<1 || (outstr!=NULL && n<2))
        error("Not enough blocking factors supplied");
    fts_setiblk(blocking[0]);    /* set input blocking factor */
    fts_setoblk(blocking[1]);    /* set input blocking factor */

    Qmultiple = getbparam("multiple");
    naxis3 = getiparam("naxis3");

    if (naxis3 < 1) {
      for (j=0; j<ninfiles; j++) {        /* scan all files for compatibility */
        instr = stropen(innames[j],"r");
        for(;;) {
            fts_zero(&fh);
            n = fts_rhead(&fh,instr);       /* read header */
            if (n<0) break;                 /* done reading */
            if (naxis3 < 1) {        /* first time around */
                naxis  = fh.naxis;
                naxis1 = fh.naxisn[0];
                naxis2 = fh.naxisn[1];
                naxis3 = naxis > 2 ? fh.naxisn[2] : 1;
                datamin = fh.datamin;
                datamax = fh.datamax;
            } else {
                if (naxis1 != fh.naxisn[0])
                    error("FITS file %s has incompatible X size",innames[j]);
                if (naxis2 != fh.naxisn[1])
                    error("FITS file %s has incompatible Y size",innames[j]);
                naxis3 += fh.naxis > 2 ? fh.naxisn[2] : 1;
                datamin = MIN(fh.datamin,datamin);
                datamax = MAX(fh.datamax,datamax);
            }
            for (i=3; i<fh.naxis; i++)
                if (fh.naxisn[i] > 1)
                    error("[%d]=>%d: Cannot handle 4+D FITS file %s",
				i+1, fh.naxisn[i], innames[j]);
            if (!Qmultiple) break;          /* only need first HDU */
            fts_sdata(&fh,instr);           /* skip data */
        }
        strclose(instr);
      } /* for(j) */
    } else {
        naxis = 0;                          /* signal don't have a size yet */
    }

    printf("Output NAXIS3 = %d\n",naxis3);
    outstr = stropen(outname,"w");

    fts_zero(&fh_out);
    for (j=0; j<ninfiles; j++) {
        instr = stropen(innames[j],"r");
        for (;;) {
            fts_zero(&fh);
            n = fts_rhead(&fh,instr);
            if (n<0) break;
            if (naxis==0) {
                naxis = fh.naxis;
                naxis1 = fh.naxisn[0];
                naxis2 = fh.naxisn[1];
            }
            if (fh_out.naxis <= 0) {
            	dprintf(1,"Copying Header: %d %d %d\n",naxis1,naxis2,naxis3);
                memcpy((char *)&fh_out,(char *)&fh,sizeof(fh));
                fh_out.naxis = 3;   /* naxis */
                naxisn[0] = naxis1;
                naxisn[1] = naxis2;
                naxisn[2] = naxis3;
                fh_out.naxisn = naxisn;
                fh_out.datamin = datamin;
                fh_out.datamax = datamax;
                fts_dhead(&fh_out,minmax);
                if (Qcompact) {
                    for (i=1; i<fh_out.naxis; i++) {
                        if (fh_out.naxisn[i-1] == 1) {
                            dprintf(0,"Moving NAXIS%d\n",i);
                            ntmp = fh_out.naxisn[i];
                            fh_out.naxisn[i] = fh_out.naxisn[i-1];
                            fh_out.naxisn[i-1] = ntmp;
                        }
                    }
                }
                fts_whead(&fh_out,outstr);
            }
            fts_cdata(&fh,instr,outstr,FALSE,FALSE);
            fts_sdata(&fh,instr);
            dprintf(1,"Copying %d-D Data: %d\n",fh.naxis, fh_out.naxis);
            if (!Qmultiple) break;
        }
        strclose(instr);
    } /* for (j) */
    strclose(outstr);
}
