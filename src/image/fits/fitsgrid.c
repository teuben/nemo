/*
 *   FITSGRID:   grid a binary table into an image,
 *	 	 quick and dirty for the DIRBE data
 * 
 *		(Note: code still contains lacks proper error checking
 *		       against the incoming FITS file !!! )
 *
 *
 *     13-aug-93    V1.0  Written - for DIRBE only      PJT
 *     16-aug-93    V1.1  added ncell=
 *      2-sep-93    V1.2  added log= for table of (mean,sigma)
 *                        (now called moment=)
 *
 *  Deficiencies:
 *      - no wrapping around 0..360 if ncell > 0
 *	- edge pixels are not correct if ncell > 0
 *      - force -180..180 range if GC in the middle, ie cannot map -10..190
 *      - no weighing of data according to DIRBE detector rules (yet)
 *        i.e. entry.w is always 1.0
 */

#include <nemo.h>
#include <history.h>
#include <nemo_fitsio.h> /* for processing output file */
#include <fits.h>   /* for processing input file */

string defv[] = {
    "in=???\n              Input fits BINTABLE file (HDU=2)",
    "out=???\n             Output fits image",
    "long=180,-180\n       Range in (galactic) longitudes",
    "lat=-15,15\n          Range in (galactic) latitudes",
    "nlong=360\n           Number of pixels in longitude",
    "nlat=30\n             Number of pixels in latitude",
    "band=1\n              DIRBE Band to process (1..10)",
    "coord=galactic\n	   Coordinate system to grid in (gal|ecl)",
    "ncell=0\n             Number of neighbor cells to use",
    "sigma=0\n             Core weighting scale (if >0) for neighbor cells",
    "moment=1\n            Moment 1 (mean in cell) or 2 (dispersion) map? ",
    "VERSION=1.3c\n        12-apr-04 PJT",
    NULL,
};

string usage = "convert and grid fits table to regular fits image";


#define SYM_ANGLE(x) ((x) - 360.0 * floor(((x)+180.0)/360.0 ))
#define MAXCNT  100

typedef struct my_table_header {        /* only for DIRBE */
    int     pixel_no;           /* 'J' unscaled */
    float   photomet[10];
    short   photoqual[10];      /* 'I' unscaled */
    double  time;
    float   eclon;              /* 'I' scaled */
    float   eclat;              /* 'I' scaled */
    float   galon;              /* 'I' scaled */
    float   galat;              /* 'I' scaled */
} my_table_header;

typedef struct entry {          /* each pixel has a linked list of these: */
    float x, y;                 /* position */
    float f, w;                 /* flux and weight */
    struct entry *next, *last;  /* point to next and last (for fast append) */
} entry, *entryptr;

entryptr *grid;                 /* Nlon * Nlat grid of pointers to entries */
int band;                       /* band to use */
float lonmin, lonmax;           /* grid edges  in longitude */
float latmin, latmax;           /* grid edges  in latitude */
float dlat, dlon;               /* pixel sizes */
int nlat, nlon;                 /* number of pixels */
int nused=0;                    /* counter */
int ncell;
float sigma2;
bool gc_middle;	                /* see of Gal.Center needs to be in middle */
bool gal_coord;			/* use galactic (T) or ecliptic (F) */

/* local's */
void fitwral();     /* should go to fitsio */





void nemo_main()
{
    stream instr, outstr;
    int    i, n, naxis1, naxis2, naxis[2], moment;
    double edges[2];
    struct fits_header fh;
    struct my_table_header r;
    char   *record, *cp, mesg[80];
    string *hitem;
    FITS   *map;

/* Setup */
    
    instr = stropen(getparam("in"),"r");    /* open input */
    moment = getiparam("moment");
    if (moment < 1 || moment > 2) error("moment must be 1 or 2");

    
    band = getiparam("band");
    if (band < 1 || band > 10) {
        band = band_id(getdparam("band"));
        if (band < 1) error("Invalid DIRBE band");
    }

    naxis[0] = nlon = getiparam("nlong");
    naxis[1] = nlat = getiparam("nlat");
    grid = (entryptr *) allocate(nlat*nlon*sizeof(entryptr));
    for (i=0; i<nlat*nlon; i++)
        grid[i] = NULL;
    cp = getparam("coord");
    switch (*cp) {
      case 'g':  gal_coord = TRUE; break;
      case 'e':  gal_coord = FALSE; break;
      default: error("Bad coordinate system choosen; try gal or ecl");
    }

    if (nemoinpd(getparam("long"),edges,2) != 2) error("long= needs 2 values");
    if (edges[0] <= edges[1]) error("long= needs left edge to be largest");
    lonmin = edges[0];
    lonmax = edges[1];
    dlon = (lonmax-lonmin)/(float)nlon;
    if (nemoinpd(getparam("lat"),edges,2) != 2) error("lat= needs 2 values");
    if (edges[0] >= edges[1]) error("lat= needs right edge to be largest");
    latmin = edges[0];
    latmax = edges[1];
    dlat = (latmax-latmin)/(float)nlat;
    dprintf(1,"GridSize: %d * %d Pixels: %g * %g\n",nlon,nlat,dlon,dlat);
    gc_middle = (lonmax < 0.0 && lonmin > 0.0); /* see if to use SYM_ANGLE */
    ncell = getiparam("ncell");
    sigma2 = 2*sqr(getdparam("sigma"));

/* Open output FITS file, and write a small yet descriptive enough header */
    
    map = fitopen(getparam("out"),"new",2,naxis);
    fitwrhda(map,"CTYPE1", gal_coord ? "GLON" : "ELON");
    fitwrhdr(map,"CRPIX1",(float) 1.0);     /* should use center */
    fitwrhdr(map,"CRVAL1",(float) (lonmin + 0.5 * dlon));
    fitwrhdr(map,"CDELT1",(float) dlon);
    
    fitwrhda(map,"CTYPE2", gal_coord ? "GLAT" : "ELAT");
    fitwrhdr(map,"CRPIX2",(float) 1.0);     /* should use center */
    fitwrhdr(map,"CRVAL2",(float) (latmin + 0.5 * dlat));
    fitwrhdr(map,"CDELT2",(float) dlat);

    fitwrhda(map,"TELESCOP","COBE");
    fitwrhda(map,"INSTRUME","DIRBE");
    fitwrhda(map,"ORIGIN","NEMO processing on CDAC data");
    fitwrhda(map,"BUNIT","MJy/sr");


    sprintf(mesg,"NEMO: %s VERSION=%s",getargv0(),getparam("VERSION"));
    fitwra(map,"HISTORY", mesg);
    hitem = ask_history();
    fitwra(map,"HISTORY","NEMO: History in reversed order");
    for (i=0, cp=hitem[0]; cp != NULL; i++) {
        fitwral(map,"HISTORY",cp);
        cp = hitem[i+1];		/* point to next history item */
    }
        

/* Open input file, and process all rows */
    
    fts_zero(&fh);		               /* clean out header */
    n = fts_rhead(&fh,instr);	               /* read primary header */
    if (n<0) error("Error reading primary HDU");
    fts_sdata(&fh,instr);                      /* and skip data .. */

    fts_zero(&fh);                             /* clean out header */
    n = fts_rhead(&fh,instr);	               /* read primary header */
    if (n<0) error("Error reading second HDU");
    naxis1 = fh.naxisn[0];                      /* size of one row */
    naxis2 = fh.naxisn[1];                      /* number of rows */
    record = allocate(naxis1);
    for (i=0; i<naxis2; i++) {                  /* loop reading rows */
        n = fread(record,1,naxis1,instr);
        if (n != naxis1) error("Early EOF on record %d",i+1);
        stuffit(&fh,record,&r);
    }
    printf("Used %d/%d points in gridding\n",nused,naxis2);

/* map the data on a grid */

    mapit(map,moment);

/* finish off */

    fitclose(map);
    strclose(instr);
}


stuffit(fh, record, th)
struct fits_header *fh;
char *record;
my_table_header *th;
{
    short tmp;
    int ilat, ilon, igrid;
    float x, y;
    entry *ep;
    
    memcpy(&th->pixel_no,   &record[0], 4);          /* col 1 */
    memcpy(th->photomet,    &record[4], 40);         /* col 2 */
    memcpy(th->photoqual,   &record[44],20);         /* col 3 */
    memcpy(&th->time,       &record[64], 8);         /* col 4 */
    memcpy(&tmp,            &record[72], 2);         /* col 5 */
    th->eclon = tmp * fh->tscaln[4] + fh->tzeron[4];
    memcpy(&tmp,            &record[74], 2);         /* col 6 */
    th->eclat = tmp * fh->tscaln[5] + fh->tzeron[5];
    memcpy(&tmp,            &record[76], 2);         /* col 7 */
    th->galon = tmp * fh->tscaln[6] + fh->tzeron[6];
    memcpy(&tmp,            &record[78], 2);         /* col 8 */
    th->galat = tmp * fh->tscaln[7] + fh->tzeron[7];

    x = gal_coord ? th->galon : th->eclon;
    y = gal_coord ? th->galat : th->eclat;
    if (gc_middle) x = SYM_ANGLE(x);

    ilon = (int) floor( (x - lonmin)/dlon);
    if (ilon < 0 || ilon >= nlon) return;       /* outside grid */

    ilat = (int) floor( (y - latmin)/dlat);
    if (ilat < 0 || ilat >= nlat) return;       /* outside grid */
    
    igrid = ilon + ilat*nlon;       /* 1d-index into grid */
    ep = grid[igrid];
    if (ep==NULL) {         /* first entry in cell : initialize */
        ep = (entry *) allocate(sizeof(entry));
        ep->next = NULL;
        ep->last = ep;      /* point to itself ... only for first in list */
        grid[igrid] = ep;
    } else {                /* follow linked list, and append */
        ep->last->next = (entry *)allocate(sizeof(entry));
        ep->last       = ep->last->next;
        ep             = ep->last;
        ep->next = NULL;
        ep->last = NULL;    /* since this is the last */
    }
    ep->x = x;
    ep->y = y;
    ep->f = th->photomet[band-1];
    ep->w = 1.0;    /* for now ... */
    nused++;
}




mapit(map,moment)
FITS *map;
int moment;
{
    entry *ep;
    float w, sum0, sum1, sum2, *row, xcell, ycell, xcell1, ycell1;
    int i, j, ilat, ilon;
    int ihisto, histo[MAXCNT+1];
    float weight();

    row = (float *) allocate(sizeof(float)*nlon);
    for (i=0; i<=MAXCNT; i++) histo[i] = 0;

    for (j=0; j<nlat; j++) {                    /* loop over all cells */
        ycell = latmin + (0.5+j)*dlat;
        for (i=0; i<nlon; i++) {
            xcell = lonmin + (0.5+i)*dlon;
            ihisto = 0;         /* count number of entries in cell */
            sum0 = sum1 = sum2 = 0.0;

            for (ilat=j-ncell; ilat<=j+ncell; ilat++) { /* loop over neighbors */
              if (ilat < 0 || ilat >= nlat) continue;
              for (ilon=i-ncell; ilon<=i+ncell; ilon++) {
                if (ilon < 0 || ilon >= nlon) continue;
                /* could allow wrapping here if full 0--360 range is used */
                ep = grid[ilon + ilat*nlon];
                while (ep != NULL) {
                    w = weight(ep,xcell,ycell,sigma2);
                    sum0 += w;
                    sum1 += w * ep->f;
                    sum2 += w * ep->f * ep->f;
                    if (ihisto < MAXCNT) ihisto++;
                    ep = ep->next;
                }
              }
            }
            histo[ihisto]++;
            row[i] = 0.0;           /* value for no-data or bad */
            if (ihisto > 0) {
                sum1 /= sum0;
                if (moment == 1) 
            	    row[i] = sum1;
            	if (ihisto > 1 && moment==2) {
                    sum2 /= sum0;
                    sum2 -= sum1*sum1;
                    if (sum2 > 0.0) row[i] = sqrt(sum2);
            	}
            } 
        }
        fitwrite(map,j,row);
    }
    free(row);

    dprintf(0,"Histogram of cell population:\n");
    dprintf(0,"#Cells:    #Entries:\n");
    dprintf(0,"__________________\n");
    dprintf(0,"%5d      0\n",histo[0]);
    for (i=1; i<MAXCNT; i++)
        if (histo[i]>0) dprintf(0,"%5d      %d\n",histo[i],i);
    dprintf(0,"%5d    > %d\n",histo[MAXCNT],MAXCNT);

}

float weight(ep,xcell,ycell,sigma2)
entry *ep;
float xcell, ycell, sigma2;
{
    float r;
    
    if (sigma2 <= 0.0) return 1.0;      /* plain average */

    r = sqr(xcell - ep->x) + sqr(ycell - ep->y);
    r /= sigma2;
    if (r > 5) return 0.0;      /* appr. < 0.006 */
    return exp(-r);
}

static real dirbe_freq[] = 
    { 1.25, 2.2, 3.5, 4.9, 12.0, 25.0, 60.0, 100.0, 140.0, 240.0, -1.0 };

/*
 * BAND_ID: return DIRBE band  (1..10) if one if found by wavelenght
 *                 0 if an error, or not found
 */
 
int band_id(freq)
real freq;
{
    int i, imin;
    real d, dmin;

    if (freq <= 0.0) {
        warning("Cannot match wavelength %g",freq);
        return 0;
    }
    i=0;
    while (dirbe_freq[i] > 0) {                     /* find match */
        if (freq == dirbe_freq[i]) {
            dprintf(0,"Using band %d: exact match found at wavelength %g \\mu",
                    i+1, freq);
            return i+1;
        }
        i++;
    }
    i = imin = 0;
    dmin = HUGE;
    while (dirbe_freq[i] > 0) {                 /* find close match */
        d = ABS(freq-dirbe_freq[i]);
        if (d < dmin) {
            dmin = d;
            imin = i+1;
        }
        i++;
    }
    if (imin > 0 && ABS(1-d/freq) < 0.25) {  /* accept if near match < 25% */
        dprintf(0,"Using band %d: near match found at wavelength %g \\mu",
                    imin, dirbe_freq[imin-1]);
        return imin;
    }
    warning("No match found at wavelength %g, nearest was %g",freq,
        imin > 0 ? dirbe_freq[imin-1] : -1);
    return 0;
}

/*
 * output a very long string, spanning multiple lines
 * since the fits community hasn't made up it's mind
 * this routine has to be used with care.
 *
 */

void 
fitwral(map, key, value)
FITS *map;
string key, value;
{
    char tmp[80], *cp = value;
    int tmplen = 70;

    for (;;) {
      strncpy(tmp,cp,tmplen);
      tmp[tmplen] = 0;
      if (cp==value)
        fitwra(map,key,tmp);
      else
        fitwra(map," ",tmp);
      if ((int)strlen(cp) < tmplen) break;
      cp += tmplen;
    } 
}
