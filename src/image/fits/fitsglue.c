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
 *     20-jun-01   1.3  removed the blocking= keyword           pjt
 *     10-aug-02   1.4  copy over WCS of the first file         pjt
 *      3-oct-02   1.4a flush 0s to fill to 2880 
 *                        fts_cdata_(....,TRUE) does not work   PJT
 *     14-jan-03   1.5  fix when CROTA is present               PJT
 *
 * TODO:
 *   if no WCS, wcs=t coredumps the program
 *   CROTA, various CD matrix
 */

#include <stdinc.h>
#include <getparam.h>
#include <fits.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input fits files or template for list",
    "out=???\n	           Output fits file",
    "multiple=t\n          Check for multiple fits files on input?",
    "naxis3=0\n            Preset NAXIS3 here, or let program scan first",
    "compact=f\n           Compact (move) dummy axes to the end",
    "inlist=\n             optional nemoinp(1) list expression for in=",
    "wcs=f\n               try and copy a reasonably WCS from input to output",
    "VERSION=1.5a\n        4-may-04 PJT",
    NULL,
};

string usage = "catenate fits images into a fitscube";

extern string *burststring(string,string);

#ifndef MAXIN
#define MAXIN 8096
#endif

local string minmax[] = { "DATAMIN", "DATAMAX", NULL };

typedef struct fts_wcs {
  int naxis;
  double crval;
  double crpix;
  double cdelt;
  double crota;    /* ieck */
  char ctype[32];
  char cunit[32];  /* optional */
} fts_wcs;

void fts_iniwcs(int, fts_wcs *);
void fts_addwcs(fits_header *,int ,fts_wcs *);
void fts_setwcs(fits_header *,int ,fts_wcs *);

void nemo_main()
{
    stream instr, outstr;
    int    i,j, n, planesize, totsize;
    int    zbufsiz, ninfiles, naxis, naxis1, naxis2, naxis3, naxisn[3], ntmp;
    bool   Qmultiple, Qcompact = getbparam("compact"), Qwcs = getbparam("wcs");
    struct fits_header fh, fh_out;
    string outname, *innames;
    char   *buffer, *zbuf, fname[128], zero = 0;
    int    inlist[MAXIN];
    double datamin, datamax;
    fts_wcs wcs[3];

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

    Qmultiple = getbparam("multiple");
    naxis3 = getiparam("naxis3");
    fts_iniwcs(3,wcs);

    if (naxis3 < 1) {
      for (j=0; j<ninfiles; j++) {        /* scan all files for compatibility */
        instr = stropen(innames[j],"r");
        for(;;) {
            fts_zero(&fh);
            n = fts_rhead(&fh,instr);       /* read header */
            if (n<0) break;                 /* done reading */
	    if (Qwcs) fts_addwcs(&fh,3,wcs);
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
    for (j=0; j<ninfiles; j++) {              /* now accumulate data */
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
		if (Qwcs) fts_setwcs(&fh_out,3,wcs);
                fts_whead(&fh_out,outstr);
            }
	    fts_cdata(&fh,instr,outstr,FALSE,FALSE);
            fts_sdata(&fh,instr);
            dprintf(1,"Copying %d-D Data: %d\n",fh.naxis, fh_out.naxis);
            if (!Qmultiple) break;
        }
        strclose(instr);
    } /* for (j) */
    n = fts_tsize(&fh_out);
    if (n > 0)
      fwrite(&zero,1,n,outstr);
    strclose(outstr);
}

/* 
 * this WCS stuff should be absorbed into fits.c
 */

void fts_iniwcs(int ndim, fts_wcs *wcs)
{
  int i;

  for (i=0; i<ndim; i++)
    wcs[i].naxis = 0;
}

void fts_addwcs(fits_header *fh, int ndim, fts_wcs *wcs)
{
  int i;

  if (wcs[0].naxis == 0) {  /* first time around */
    if (fh->naxis > 3) error("Cannot handle > 3D data yet");
    for (i=0; i < fh->naxis; i++) {
      wcs[i].naxis = fh->naxisn[i];
      if (fh->crpixn)
	wcs[i].crpix = fh->crpixn[i];   /* need to check if they exist */
      else 
	wcs[i].crpix = 1.0;
      if (fh->crvaln)
	wcs[i].crval = fh->crvaln[i];
      else
	wcs[i].crval = 0.0;
      if (fh->cdeltn)
	wcs[i].cdelt = fh->cdeltn[i];
      else
	wcs[i].cdelt = 1.0;
      if (fh->crotan)
	wcs[i].crota = fh->crotan[i];
      else
	wcs[i].crota = 0.0;
      if (fh->ctypen)
	strcpy(wcs[i].ctype, fh->ctypen[i]);
      else {
	if (i==0)
	  strcpy(wcs[i].ctype,"RA");
	else if (i==1)
	  strcpy(wcs[i].ctype,"DEC");
	else if (i==2)
	  strcpy(wcs[i].ctype,"VLSR");
      }
    }
    if (wcs[2].naxis == 0) {
      wcs[2].naxis = 1;
      wcs[2].crval = 0;
      wcs[2].cdelt = 1;
      wcs[2].crpix = 1;
      wcs[2].crota = 0;
      strcpy(wcs[2].ctype,"VLSR");
    }
    dprintf(0,"fts_addwcs: Starting at %d x %d x %d [%s %s]\n",
	    wcs[0].naxis, wcs[1].naxis, wcs[2].naxis,
	    wcs[0].ctype, wcs[1].ctype);
  } else {                  /* accumulate */
    if (fh->naxis > 2) {
      wcs[2].naxis += fh->naxisn[2];      
    } else {
      wcs[2].naxis++;
    }
    dprintf(0,"fts_addwcs: Continue to %d x %d x %d\n",
	    wcs[0].naxis, wcs[1].naxis, wcs[2].naxis);
  }
}

void fts_setwcs(fits_header *fh, int ndim, fts_wcs *wcs)
{
  int i,n;

  n = fh->naxis;
  dprintf(0,"fts_setwcs: naxis=%d\n",n);
  if (n>3) error("Cannot do > 3D yet");
  fh->crvaln = (double *) allocate(n*sizeof(double));
  fh->crpixn = (double *) allocate(n*sizeof(double));
  fh->cdeltn = (double *) allocate(n*sizeof(double));
  fh->crotan = (double *) allocate(n*sizeof(double));
  fh->ctypen = (char **) allocate(n*sizeof(char *));

  for (i=0; i<n; i++) {
    if (fh->naxisn[i] != wcs[i].naxis) 
      warning("fts_setwcs: mismatch at axis %d: %d != %d",
	      i+1,fh->naxisn[i],wcs[i].naxis);
    fh->crvaln[i] = wcs[i].crval;
    fh->crpixn[i] = wcs[i].crpix;
    fh->cdeltn[i] = wcs[i].cdelt;
    fh->ctypen[i] = wcs[i].ctype;
    fh->crotan[i] = wcs[i].crota;
  }
}
