/*
 *  TSD: show what's in a HDF Scientific Data Set (SD),
 *		modeled after 'tsf'
 *
 *      24-dec-94	V1.0  toy model		Peter Teuben
 *                      the MULTI_FILE doesn't work yet
 *	12-feb-95       added verbatim type info	PJT
 *	19-may-95	local symbols - 'fmt' clashed with something on linux
 *      21-may-95       V1.2: added out= format=, coordinate output
 *      27-aug-96       V1.3: header output now using dprintf()
 *       6-dec-04       V1.4: warn and disable SDS maps that differ in size
 *                            allow user to select different SDS
 *      10-dec-04       V1.5: dummy= introduced
 *      14-dec-05       V1.6  changed order of keywords, select= now 3rd as with other hdf tools
 *
 *  TODO: fix rank=3 with coordinates and select=
 *     
 */

 
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <history.h>

//#if defined(HAVE_HDF)
#include <hdf.h> 	/* some conflicts with nemo include files */
//#endif

string defv[] = {
    "in=???\n			Input file (HDF SD)",
    "out=\n                     ascii dump of the data to this file?",
    "select=all\n               Select which SDS# for display? (all|1..)",
    "format=%g\n                Format used in dump",
    "coord=f\n                  Add coordinates?",
    "dummy=t\n                  Also print out dummy axis (axis with length 1)",
    "VERSION=1.6\n		14-dec-04 PJT",
    NULL,
};

/* #define MULTI_FILE 1 == does not work yet */

string usage="Scan and optionally ascii dump of an HDF SDS";

string cvsid="$Id$";
 

#define MAXRANK 10
#define MAXSDS  32

local int rank, shape[MAXRANK], run[MAXRANK], old_shape[MAXRANK];
local char label[256], unit[256], fmt[256], coordsys[256];

local void scan_sd(string infile);

extern void hdf_type_info(int type, char *msg);     // see hdf_info.c

  
void nemo_main()
{
    scan_sd(getparam("in"));
}

#ifndef MULTI_FILE
void scan_sd(string infile)		/* this is the fancy new version */
{
    float **dump, **coord;
    int i, j, k, ret, old_size, size, type, nsds, old_rank;
    char ntype[32];
    bool *visib, axis_visib[MAXRANK];
    int nselect, select[MAXSDS];
    string format, sselect;
    stream outstr;
    bool Qdummy = getbparam("dummy");
    bool Qout = hasvalue("out");

    if (!Qdummy) warning("dummy may not be working so well");

    nsds = DFSDndatasets(infile);
    if (nsds<0) 
        error("%s is probably not an HDF scientific dataset",infile);
    dprintf(0,"Found %d scientific data set%s in %s\n",
            nsds,(nsds > 1 ? "s" : ""),infile);
    visib = (bool *) allocate(sizeof(bool)*nsds);
    if (hasvalue("out")) {
        dump = (float **) allocate(nsds*sizeof(float *));   /* all data !! */
        format = getparam("format");
        outstr = stropen(getparam("out"),"w");
    } else
        outstr = NULL;

    for (k=0; k<nsds; k++)    /* first flag all SDS to be shown */
        visib[k] = TRUE;
    for (k=0; k<MAXRANK; k++)  /* first flag all axes to be shown (in case coord=t) */
        axis_visib[k] = TRUE;

    sselect = getparam("select");
    if (!streq(sselect,"all")) {
      nselect = nemoinpi(sselect,select,MAXSDS);
      if (nselect < 0) error("%d error parsing %s",nselect,sselect);
      if (nselect > nsds) error("%s: too many specified, nsds=%d",nselect,nsds);
      for (k=0; k<nsds; k++)
        visib[k] = FALSE;
      for (k=0; k<nselect; k++) {
	if (select[k] < 1 || select[k] > nsds) error("%d: bad SDS selection, nsds=%d",select[k],nsds);
	visib[select[k]-1] = TRUE;
      }
    }
    for (k=0; k<nsds; k++) 
      dprintf(1,"%d: %s\n", k+1, visib[k] ? "OK" : "hidden");

    old_size = -1;
    for (k=0; k<nsds; k++) {               /* loop over all SDS to get the rank, shape, size and info */
    	ret = DFSDgetdims(infile,&rank, shape, MAXRANK);
    	if (ret < 0) error("Problem getting rank/shape at SDS #%d",k+1);

    	label[0] = unit[0] = fmt[0] = coordsys[0] = 0;
        ret = DFSDgetdatastrs(label, unit, fmt, coordsys);
	if (ret < 0) error("Problem getting labels at SDS #%d",k+1);

        ret = DFSDgetNT(&type);
	if (ret < 0) error("Problem getting data type at SDS #%d",k+1);

        if (! visib[k]) continue;          /* don't count SDS# that were not selected */

        if (old_size < 0) {                /* first time around allocate coordinates */
            if (getbparam("coord")) {
                coord = (float **) allocate(rank*sizeof(float *));
            } else {
                coord = NULL;
            }
        }

    	dprintf(0,"%d: %s(",k+1,label);
    	for (i=0, size=1; i<rank; i++) {
    	    if (i==rank-1)
                dprintf(0,"%d)",shape[i]);
            else
                dprintf(0,"%d,",shape[i]);
    	    size *= shape[i];

            if (old_size < 0 && coord) {
                coord[i] = (float *) allocate(shape[i] * sizeof(float));
                ret = DFSDgetdimscale(i+1, shape[i],coord[i]);
                if (ret<0) error("getting shape[%d]",i+1);
            }
    	}
    	hdf_type_info(type,ntype);
    	dprintf(0," %s ",unit);
    	dprintf(0," -> [%d elements of type: %d (%s)]\n", size, type, ntype);
        if (outstr) {
            dump[k] = (float *) allocate(size * sizeof(float));
            ret = DFSDgetdata(infile,rank,shape,dump[k]);
        }

	if (old_size < 0) {       /* first time around */
	  old_size = size;
	  old_rank = rank;
	  for (i=0; i<rank; i++) old_shape[i] = shape[i];
	} else {                   /* make sure subsequent ones have the same size for display */
	  if (old_size != size) {
	    if (Qout) warning("bad shape for SDS #%d, removing from selection list",k+1);
	    visib[k] = FALSE;
	    size = old_size;
	  } else if (old_rank != rank) {
	    if (Qout) warning("bad rank for SDS #%d, removing from selection list",k+1);
	    visib[k] = FALSE;
	    rank = old_rank;
	  }
	}
    } /* k */
    rank = old_rank;   /* restore it, in case the last one was a bad one */
    for (i=0; i<rank; i++) shape[i] = old_shape[i] ;

#if 0
    for (k=0; k<nsds; k++)
      if (!visib[k]) warning("SDS #%d cannot be displayed, it has a different size",k+1);
#endif

    if (outstr) {
        for (i=0; i<rank; i++) {
	  run[i] = 0;      /* reset run (coordinate index) array */
	  axis_visib[i] = Qdummy ? TRUE : shape[i] > 1;
	  dprintf(1,"axis %d=%d => %d\n",i+1,shape[i],axis_visib[i]);
	}

        for (i=0; i<size; i++) {                /* loop over all data */
            if (coord) {                        /* print coord system ? */
                for (j=rank-1; j>=0; j--) {
		  if (Qdummy || axis_visib[j]) {
                    fprintf(outstr,format,coord[j][run[j]]);  
                    fprintf(outstr," ");
		  }
                }
                run[rank-1]++;
                for (j=rank-1; j>=0; j--) {     /* check if axis reset needed */
                    if (run[j] >= shape[j]) {
                        run[j] = 0;
                        if (j>0) run[j-1]++;
                    } else
                        break;
                }
            } /* coord */
	
            for (k=0; k<nsds; k++) {            /* loop over all columns */
	      if (visib[k]) {
                fprintf(outstr,format,dump[k][i]);
                fprintf(outstr," ");
	      }
            }        
            fprintf(outstr,"\n");
        }
    }
}

#else


     /* NOTE: the above code was substantially changed, the code below not */

#define DFACC_RDONLY 1

void scan_sd(string infile)
{
    int id, sds, n, i, j, k, ret, size, type, na, nd;
    char name[256];

    id = SDstart(infile, DFACC_RDONLY);
    if (id<0) error("%s is probably not an HDF SD (scientific dataset)",
                      infile);

    ret = SDfileinfo(id, &nd, &na);
    
    dprintf(0,"Found %d scientific data set%s in %s with %d attributes\n",
            nd,(nd > 1 ? "s" : ""),infile, na);
    for (k=0;;k++) {
        sds = SDselect(id, nd);
        ret = SDgetinfo(sds, name, &rank, shape, &type, &na);
    	label[0] = unit[0] = fmt[0] = coordsys[0] = 0;
        ret =  SDgetdatastrs(sds, label, unit, fmt, coordsys, 256);
    	dprintf(0,"%d: %s(",k,label);
    	for (i=0, size=1; i<rank; i++) {
    	    if (i==rank-1)
                dprintf(0,"%d)",shape[i]);
            else
                dprintf(0,"%d,",shape[i]);
    	    size *= shape[i];
    	}
    	dprintf(0," %s ",unit);
    	dprintf(0," -> [%d elements of type: %d]\n", size, type);
        /* ret = SDendaccess(sds); */
    }   
    SDend(id);
}

#endif
/*

   32-bit float DFNT_FLOAT32         5
   64-bit float DFNT_FLOAT64         6
    8-bit signed int DFNT_INT8      20
    8-bit unsigned int DFNT_UINT8   21
   16-bit signed int DFNT_INT16     22
   16-bit unsigned int DFNT_UINT16  23
   32-bit signed int DFNT_INT32     24
   32-bit unsigned int DFNT_UINT32  25

*/                        
