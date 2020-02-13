/*
 *   TABFITS:   turn a table into a fits file 
 *
 *	4-dec-98	coded, quick and dirty as usual		PJT
 *		  ** oops, still ccd, not fits yet **
 *      7-dec-98    read commented header instead               pjt
 *     17-dec-98    V2.0:  nxy= is now nx=, ny=                 pjt
 *      5-jan-98    V2.1 fixed header item parsing              pjt
 *     30-nov-99    V3.0: allow nz=                             pjt
 *      9-jul-00    V3.1: allow dcol=0 for reading all columns  pjt
 *     22-jul-00    V3.2: fix for new style fitsio              pjt
 *     23-oct-03    V3.2b: check for MAXDAT (and made it bigger) pjt
 *      1-jan-04        c: get_line check return value
 *     18-may-06        d: fix for too long comment lines        pjt
 *      7-nov-2013  V3.3  added scale=                           pjt
 *     13-feb-2020  V3.4  promote COMMENT and HISTORY for Wolfire pjt
 *
 * todo:  read one more line, check if file is done!!
 *        if header contains NAXIS, NAXIS1,...., nx=,ny=,nz= is still needed
 *        handle COMMENT and HISTORY "keywords" like keys, not comments
 */

#include <stdinc.h>
#include <getparam.h>
#include <fitsio_nemo.h>
#include <ctype.h>

string defv[] = {
    "in=???\n		Input sig (table) file",
    "out=???\n		Output image file",
    "dcol=1\n           Column for data (use 0 for continues read)",
    "nx=\n              X-Size of data, or else specify NAXIS1 in header",
    "ny=\n              Y-Size of data, or else specify NAXIS2 in header",
    "nz=\n              Z-Size of data, or else specify NAXIS3 in header",
    "nmax=100000\n      Allocation space for piped I/O",
    "scale=1.0\n        Scale factor to multiply data by",
    "VERSION=3.4\n	13-feb-2020 PJT",
    NULL,
};

string usage = "convert tables into FITS images";

string cvsid="$Id$";

#define MAXCOL 2048
#define MAXDAT 8192

#ifndef MAX_LINELEN
#define MAX_LINELEN  16384
#endif

int get_key(char *line, char *key, char *value, int *his);

void init_data(int);
double get_next_data(int);

static  char line[MAX_LINELEN];
static  double data[MAXCOL];
static  stream instr;
static  double scale;

void nemo_main()
{
    stream outstr;
    int nmax, nx, ny, nz, i, j, k, ierr, naxis[3], nheader=0, his;
    real rmin, rmax, dval;
    string name = getparam("in");
    int dcol = getiparam("dcol");
    char namex[32], namey[32], namez[32], key[80], value[80];
    FITS *fitsfile = NULL;
    float rvalue, rdata[MAXDAT];

    if (dcol < 0 || dcol > MAXCOL) error("Bad dcol=%d",dcol);
    dcol--;
    
    nmax = nemo_file_lines(name, getiparam("nmax"));

    instr = stropen(name,"r");

    nx = naxis[0] = hasvalue("nx") ? getiparam("nx") : 0;
    ny = naxis[1] = hasvalue("ny") ? getiparam("ny") : 0;
    nz = naxis[2] = hasvalue("nz") ? getiparam("nz") : 0;

    scale = getdparam("scale");


    while (1) {             /* read commented header values */

        if (fitsfile==NULL && nx*ny*nz > 0) {
            dprintf(0,"Creating initial fits primary HDU %d*%d*%d\n",nx,ny,nz);
            fitsfile = fitopen(getparam("out"),"new",3,naxis);
            if (fitsfile==NULL)
                error("Could not open fitsfile for writing");
        }
    
        if (get_line(instr, line) < 0)     /* EOF */
	  break;

        if (!get_key(line,key,value,&his)) {                /*   check if comment is # key = value */
            dprintf(0,"Read %d header lines\n",nheader);
	    break;
	}
        nheader++;
        
        dprintf(1,"HEADER: %-8s = %s\n",key,value);
        if (streq(key,"NAXIS1")) {
            nx = naxis[0] = atoi(value);
        } else if (streq(key,"NAXIS2")) {
            ny = naxis[1] = atoi(value);
        } else if (streq(key,"NAXIS3")) {
            nz = naxis[2] = atoi(value);
        } else if (strncmp(key,"CRVAL",5)==0 ||
		   strncmp(key,"CDELT",5)==0 ||        
		   strncmp(key,"CRPIX",5)==0) {
            /* write real fits header */
	    rvalue = atof(value);
            fitwrhdr(fitsfile,key,rvalue);		/* casting problem ??*/
        } else if (strncmp(key,"CTYPE",5)==0) {
            fitwrhd(fitsfile,key,value);
	} else if (his) {
	  // HISTORY or COMMENT
	  fitwra(fitsfile,key,value);
        } else if (*key) {
            /* write keyword verbose */
            /* fitwrhda(fitsfile,key,value); */
	  if (strlen(key) > 8) error("can't write key=%s",key);
	  if (strlen(value) > 60) value[60] = 0;
	  fitwrhd(fitsfile,key,value);
        } else {
	  if (strlen(line) > 72) line[72]=0;
	  /*               12345678        */
	  fitwra(fitsfile,"        ",line);
      
        }
    }

    if (nx > MAXDAT) error("MAXDAT=%d not large enough for nx=%d",MAXDAT,nx);

    if (nx*ny*nz <= 0 || fitsfile == NULL) {         /* figure out sizes */
    	if (nx == 0 && ny == 0 && nz == 0) {
    	    error("Need to specify at least nx= or ny= or nz=");
        } else if (nx == 0 && nz == 0) {
            nx = nmax/ny;
            nz = 1;
            if (nmax % ny) warning("%d odd trailing lines in file",nmax%ny);
        } else if (ny == 0 && nz == 0) {
            ny = nmax/nx;
            nz = 1;
            if (nmax % nx) warning("%d odd trailing lines in file",nmax%nx);
        } 

        if (nx*ny*nz <= 0)
            error("[%d,%d,%d] I give up, cannot determine size of array",
                        nx,ny,nz);

        naxis[0] = nx;
        naxis[1] = ny;
        naxis[2] = nz;

        dprintf(0,"Creating fits primary HDU %d*%d*%d\n",nx,ny,nz);
        fitsfile = fitopen(getparam("out"),"new",3,naxis);
        if (fitsfile==NULL)
            error("Could not open fitsfile for writing");
                
    }

    init_data(dcol);
    
    for (k=0; k<nz; k++) {
      fitsetpl(fitsfile,1,&k);
      for (j=0; j<ny; j++) {                       /* loop over whole array */
        for (i=0; i<nx; i++) {
            dval = get_next_data(dcol);
            rdata[i] = (float) dval * scale;
            dprintf(1,"%d %d %d = %g\n",i,j,k,dval);
            if (i || j || k) {
                rmin = MIN(rmin, dval);
                rmax = MAX(rmax, dval);
            } else
                rmin = rmax = dval;
        }
        fitwrite(fitsfile,j,rdata);
      }
    }
    dprintf(0,"MinMax in map: %g - %g\n",rmin,rmax);    
    fitclose(fitsfile);
}

int get_key(char *line, char *key, char *value, int *his)
{
    char *cp = line, *lkey=key, *lvalue=value;

    *key = 0;                               /* initialize no keyword */
    *value = 0;                             /* and no value part */
    *his = 0;

    if (cp==0) return 1;
    while (isspace(*cp))                    /* skip leading whitespace */
        cp++;
    if (cp==0) return 1;                    /* if only whitespace, done */
    
    if (*cp == '#')			    /* expect a comment for header lines */
	cp++;    
    else
    	return 0;                           /* if not, return bad line */
    	
    while (isspace(*cp))                    /* skip leading whitespace before key */
        cp++;  

    while (!isspace(*cp) && *cp != '=')     /* grab the keyword */
        *lkey++ = *cp++;
    
    while (isspace(*cp))                    /* skip over more potential whitespace */
        cp++;
    
    if (*cp == '=')			
        *lkey = 0;			    /* terminate keyword */
    else {
        *lkey = 0;			    /* terminate keyword */	  
        dprintf(2,"PJT [%s]\n",key);
	if (streq(key,"COMMENT") || streq(key,"HISTORY")) {
	  *his = 1;
	} else {
	  *key = 0;			    /* no keyword */
	  return 1;
	}
    }
    if (*his == 0)
      while (isspace(*cp) || *cp == '=')      /* skip whitespace around '=' */
        cp++;
    else
      while (isspace(*cp))                    /* skip whitespace */
        cp++;
    while (*cp)                             /* grab the value */
        *value++ = *cp++;
    *value = 0;

    if (strlen(key) > 8) warning("key %s too long for FITS",key);

    return 1;
  
}


static int active_row = 0;
static int active_col = 0;
static int icol = 0;

void init_data(int dcol)
{
    int ierr;

    ierr=nemoinpd(line,data,MAXCOL);
    if (ierr < 1) error("badly formatted line: \"%s\" err=%d\n", line,ierr);
    if (dcol >= ierr) error("Not enough columns in data: %d < %d",ierr,dcol);
    
    active_col = ierr;
    active_row = 1;
    icol = 0;

}

double get_next_data(int dcol)
{
    if (!active_row && !active_col) {
        while(1) {
            if (get_line(instr, line) <= 0)
            	error("No more data available");
            if (line[0]!='#') break;
        }
        init_data(dcol);
    }

    if (dcol>=0) {
        active_row = active_col = 0;
        return data[dcol];
    } else {
	active_col--;
	active_row = 0;
	return data[icol++];
    }
    
}

