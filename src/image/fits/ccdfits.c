/* 
 *      CCDFITS: convert CCD image to a fits file
 *
 *	Also contains a controversial 'long strings' implementation
 *
 *      29-apr-88  V1.0 created (again) PJT
 *       2-jun-88  V1.1 renamed 'wfits' to 'ccdfits'  
 *                      new filestruct: code is almost same
 *      23-dec-88  V1.2 added velocity information to header
 *      30-jan-89  V1.3 velocity information is now in Z
 *       6-mar-89  V1.3a 2D cube written when NZ=1
 *      28-jul-90  V2.0 Bob's new fitsio.c routines, finally
 *       1-oct-90  V2.1 default scale is now 1, not 1/60   PJT
 *                      crpix is now real !!
 *	11-oct-90  V2.2 blocking factor introduced 	   PJT
 *      14-feb-91  V2.3 write axes names for AIPS          PJT
 *	 7-mar-92  V2.4 gcc2.0 happy - new style	   PJT
 *	23-feb-93  V2.5  axisnames (CTYPE's) now ok	   PJT
 *      18-may-93  V2.6 history trial following sci.astro.fits proposal  PJT
 *      21-feb-94  V2.6a ansi
 *      22-nov-94      b ???
 *	20-mar-95      c long histories now written multiline		PJT
 *	 8-oct-95      d force CTYPEs for some WCS freaks		pjt
 *      28-mar-98      e write out Headline[] as COMMENT                PJT
 *      11-mar-99      e fixed bug in the above change                  pjt
 *	21-mar-99  V2.7  fixed scaling problems in bitpix=8		pjt
 *      18-may-99  V2.8  cdmatrix optional output                       pjt
 *      15-oct-99  V2.9  force a more silly RA---SIN/DEC--SIN axis      pjt
 *      23-mar-01  V3.0  allow fits reference image to inherit WCS from PJT
 *       8-apr-01      a fixed SINGLEPREC operation
 *      23-may-01  V3.1  added dummy= to be able to not write dummy axes  PJT
 *                        NOT WORKING YET
 *	 7-aug-01  
 *      18-dec-01  V4.0  work with new fitsio_nemo.h
 *       6-may-02      b implemented dummy=f
 *       6-jul-02   5.0  changed ref*= to crval/cdelt/crpix=
 *      10-jul-02   5.1  better handling of long COMMENT fields
 *       4-feb-04   5.2  also listen to changed crval/cdelt/crpix= without refmap
 *       8-may-05   5.3  deal with the new  axis type 1 images          PJT
 *      17-aug-12   5.5  default restfreq to keep some WCS routines happy    PJT
 *       8-sep-15   5.7  add dummy 4th stokes axis to keep CASA/ADMIT happy  PJT
 *      26-may-16   5.8  added CUNITn and BUNIT, and better WCS output when radecvel=t
 *       6-apr-17   5.9  allow refmap to update only certain WCS (1,2,3)     PJT
 *       8-apr-17   6.0  new approach inheriting a new WCS
 *      14-jun-19   6.0a correct VSYS when in freq=t mode, fix cdelt1 in one common case
 *      19-jun-10   6.1  Output now in km/s
 *
 *  TODO:
 *      reference mapping has not been well tested, especially for 2D
 *
 *      WCS should follow the following strategy:
 *      - use the (usually crappy) WCS from the NEMO image (radecval=True)
 *      - inherit certain refaxis from a refmap
 *      - override any of crpix/crval/cdelt/
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <history.h>
#include <fitsio_nemo.h>

string defv[] = {
        "in=???\n        Input image filename",
        "out=???\n       Output fits filename",
	"bitpix=-32\n	 FITS bitpix value {16,32,-32}",
        "scale=1,1,1\n   Extra scalefactor for cdelt&crval",
        "iscale=1,0\n    Scale and Offset Intensity conversion factor",
        "object=\n       Object name",
        "comment=\n      Extra comments",
	"cdmatrix=f\n    Use standard CD matrix instead of CDELT?",
	"blocking=1\n	 Blocking factor for output (blocksize/2880)",
	"refmap=\n       reference map to inherit WCS from",
	"refaxis=\n      which axes from refmap to be used (1,2,3,4)",
	"crpix=\n        reference pixel, if different from default",
	"crval=\n        reference value, if different from default",
	"cdelt=\n        pixel value increment, if different from default",
	"radecvel=f\n    Enforce reasonable RA/DEC/VEL axis descriptor",
	"proj=SIN\n      Projection type if RA/DEC used (SIN,TAN)",
	"restfreq=115271204000\n   RESTFRQ (in Hz) if a doppler axis is used",  /* 1.420405751786 */
	"vsys=0\n        VSYS correction in km/s",
	"freq=f\n        Output axis in FREQ or VEL",
	"dummy=t\n       Write dummy axes also ?",
	"nfill=0\n	 Add some dummy comment cards to test fitsio",
	"ndim=\n         Testing if only that many dimensions need to be written",
	"select=1\n      Which image (if more than 1 present, 1=first) to select",
	"blank=\n        If set, use this is the BLANK value in FITS (usual NaN)",
        "VERSION=6.2\n   1-aug-2019 PJT",
        NULL,
};

string usage = "convert image to a fits file";

string cvsid = "$Id$";

stream  instr, outstr;                         /* file streams */

imageptr iptr=NULL;                     /* image, allocated dynamically */
int  isel = 0;

double scale[4];        /* scale conversion for FITS (CDELT) */
double iscale[2];	/* intensity rescale */
string object;           /* name of object in FITS header */
string comment;          /* extra comments */
string headline;         /* optional NEMO headline, added as COMMENT */
string proj;             /* projection type for WCS */
bool Qcdmatrix;          /* writing out new-style cdmatrix ? */
bool Qradecvel;          /* fake astronomy WCS header */
bool Qfreq;              /* freq or vel output ? */
bool Qrefmap;
bool Qcrval, Qcdelt, Qcrpix;
bool Qdummy;             /* write dummy axes ? */
bool Qblank;
real blankval;
int  nrefaxis, refaxis[4];
bool Qrefaxis[4];

int   nref = 0, nfill = 0;
FLOAT ref_crval[4] = {0,0,0,1},
      ref_crpix[4] = {1,1,1,1},
      ref_cdelt[4] = {1,1,1,1};
char  ref_ctype[4][80], ref_cunit[4][80];
FLOAT restfreq;
real  vsys;

void setparams(void);
void write_fits(string,imageptr);
void stuffit(FITS *, string, string);
void set_refmap(string);
void permute(int *x ,int *idx, int n);

void nemo_main()
{
  int i;

  setparams();                               /* set cmdln par's */
  instr = stropen (getparam("in"), "r");     /* open image file */
  for (i=0; i<isel; i++)
    if (read_image (instr,&iptr)==0)
      error("Cannot process image select=%d",i);
  headline = ask_headline();                 /* possible headline */
  strclose(instr);                           /* close image file */
  write_fits(getparam("out"),iptr);          /* write fits file */
  free_image(iptr);
}

void setparams(void)
{
  int i,n;
  real tmpr[4];

  isel = getiparam("select");

  switch (nemoinpd(getparam("scale"),scale,3)) {
  case 1:
    scale[1] = scale[0];
    scale[2] = 1.0;
    break;
  case 2:
    scale[2] = 1.0;
    break;
  case 3:
    break;
  case 0:
    scale[0] = scale[1] = scale[2] = 1.0;
    break;
  default:
    error("parsing error scale=%s",getparam("scale"));
  }
  if (nemoinpd(getparam("iscale"),iscale,2) != 2)
    error("parsing error scale=%s - must be 2 numbers",getparam("scale"));
  object = getparam("object");
  comment = getparam("comment");
  restfreq = getrparam("restfreq");
  Qcdmatrix = getbparam("cdmatrix");
  Qradecvel = getbparam("radecvel");
  Qfreq     = getbparam("freq");
  Qblank    = hasvalue("blank");
  if (Qblank) blankval = getrparam("blank");

  
  Qrefmap = hasvalue("refmap");
  if (Qrefmap) {
    set_refmap(getparam("refmap"));
    nrefaxis = nemoinpi(getparam("refaxis"),refaxis,4);
    for (i=0; i<4; i++) Qrefaxis[i] = FALSE;
    for (i=0; i<4; i++) Qrefaxis[refaxis[i]-1] = TRUE;
  }

  Qcrpix = hasvalue("crpix");
  if (Qcrpix) {
    nref =  nemoinpr(getparam("crpix"),tmpr,4);
    for (i=0; i<nref; i++)
      ref_crpix[i] = tmpr[i];
  }
  Qcrval = hasvalue("crval");
  if (Qcrval) {
    nref =  nemoinpr(getparam("crval"),tmpr,4);
    for (i=0; i<nref; i++)
      ref_crval[i] = tmpr[i];
  }
  Qcdelt = hasvalue("cdelt");
  if (Qcdelt) {
    nref =  nemoinpr(getparam("cdelt"),tmpr,4);
    for (i=0; i<nref; i++)
      ref_cdelt[i] = tmpr[i];
  }
  Qdummy = getbparam("dummy");
  nfill = getiparam("nfill");
  proj = getparam("proj");
  vsys = getrparam("vsys");
}

static string ctypes[4] = { "CTYPE1",   "CTYPE2",   "CTYPE3",   "CTYPE4"};
static string cdelts[4] = { "CDELT1",   "CDELT2",   "CDELT3",   "CDELT4"};
static string crvals[4] = { "CRVAL1",   "CRVAL2",   "CRVAL3",   "CRVAL4"};
static string crpixs[4] = { "CRPIX1",   "CRPIX2",   "CRPIX3",   "CRPIX4"};
static string radeves[4]= { "RA---SIN", "DEC--SIN", "VRAD",     "STOKES"};
static string radevet[4]= { "RA---TAN", "DEC--TAN", "VRAD",     "STOKES"};
//static string radeves[4]= { "RA---SIN", "DEC--SIN", "VELO-LSR", "STOKES"};
//static string radevet[4]= { "RA---TAN", "DEC--TAN", "VELO-LSR", "STOKES"};
static string radefrs[4]= { "RA---SIN", "DEC--SIN", "FREQ    ", "STOKES"};
static string radefrt[4]= { "RA---TAN", "DEC--TAN", "FREQ    ", "STOKES"};
static string xyz[4]    = { "X",        "Y",        "Z",        "S"};

void write_fits(string name,imageptr iptr)
{
    FLOAT tmpr,xmin[4],xref[4],dx[4],mapmin,mapmax;   /* fitsio FLOAT !!! */
    FLOAT bmaj,bmin,bpa;
    FLOAT fnan;
    FITS *fitsfile;
    char *cp, origin[80];
    char *ctype1_name, *ctype2_name, *ctype3_name, *ctype4_name;
    string *hitem, axname[4];
    float *buffer, *bp;
    int i, j, k, axistype, bitpix, keepaxis[4], nx[4], p[4], nx_out[4], ndim=3;
    double bscale, bzero;

    get_nanf(&fnan);    

    if (Qfreq)
      vsys = -restfreq * vsys / 300000.0;        // convert km/s to Hz
    
    if (hasvalue("ndim")) ndim = getiparam("ndim");
    nx[0] = Nx(iptr);
    nx[1] = Ny(iptr);
    nx[2] = Nz(iptr);   if (nx[2] <= 0) nx[2] = 1;
    nx[3] = 1;
    xmin[0] = Xmin(iptr)*scale[0];
    xmin[1] = Ymin(iptr)*scale[1];
    xmin[2] = Zmin(iptr)*scale[2];
    xmin[3] = 1.0;
    dx[0] = Dx(iptr)*scale[0];
    dx[1] = Dy(iptr)*scale[1];
    dx[2] = Dz(iptr)*scale[2];
    dx[3] = 1.0;
    xref[0] = Xref(iptr)+1.0;
    xref[1] = Yref(iptr)+1.0;
    xref[2] = Zref(iptr)+1.0;
    xref[3] = 1.0;
    axistype = Axis(iptr);
    axname[0] = (Namex(iptr) ? Namex(iptr) : xyz[0]);
    axname[1] = (Namey(iptr) ? Namey(iptr) : xyz[1]);
    axname[2] = (Namez(iptr) ? Namez(iptr) : xyz[2]);
    mapmin = MapMin(iptr);
    mapmax = MapMax(iptr);
    bmaj = Beamx(iptr);
    bmin = Beamy(iptr);
    bpa  = 0.0;        /* only spherical beams for now */

    if (Qdummy) 
      for (i=0; i<4; i++) p[i] = i;   /* set permute order */
    else {
      if (Qrefmap) warning("dummy=f and usage of refmap will result in bad headers");
      permute(nx,p,3);
      dprintf(0,"Reordering axes: %d %d %d\n",p[0],p[1],p[2]);
    }
#if 1
    for (i=0; i<4; i++)  nx_out[i] = nx[p[i]];
    /* fix this so CubeValue works */
    Nx(iptr) = nx_out[0];
    Ny(iptr) = nx_out[1];
    Nz(iptr) = nx_out[2];
#else
    for (i=0; i<4; i++)  nx_out[i] = nx[i];
#endif
    sprintf(origin,"NEMO ccdfits %s",getparam("VERSION"));

    dprintf(1,"NEMO Image file written to FITS disk file\n");
    dprintf(1,"%d %d %d   %f %f %f   %f %f %f  %f %f %f   %f %f \n",
	    nx[0],nx[1],nx[2],xmin[0],xmin[1],xmin[2],dx[0],dx[1],dx[2],xref[0],xref[1],xref[2],
	    mapmin,mapmax);
    dprintf(1,"keepaxis(%d,%d,%d)\n",keepaxis[0],keepaxis[1],keepaxis[2]);
    
    fit_setblocksize(2880*getiparam("blocking"));
    bitpix = getiparam("bitpix");
    fit_setbitpix(bitpix);
    if (bitpix == 16) {      /* scale from -2^(bitpix-1) .. 2^(bitpix-1)-1 */
        bscale = (mapmax - mapmin) / (2.0*32768.0 - 1.0);
        bzero = mapmax - bscale*32767.0;
        fit_setscale(bscale,bzero);
    } else if (bitpix == 32) {
        bscale = (mapmax - mapmin) / (2.0*2147483648.0 - 1.0);
        bzero = mapmax - bscale*2147483647.0;
        fit_setscale(bscale,bzero);
    } else if (bitpix == 8) {
    	bscale = (mapmax - mapmin) / (2.0*128.0 - 1.0);
    	bzero = mapmin;
    	fit_setscale(bscale,bzero);
    }
    dprintf(1,"bscale,bzero=%g %g\n",bscale,bzero);

    fitsfile = fitopen(name,"new",ndim,nx_out);
    if (fitsfile==NULL) error("Could not open fitsfile %s for writing\n",name);

    if (Qrefmap || Qcrpix) {
      dprintf(1,"Using ref_crpix\n");
      fitwrhdr(fitsfile,"CRPIX1",ref_crpix[0]);       
      fitwrhdr(fitsfile,"CRPIX2",ref_crpix[1]);       
      if (ndim>2) fitwrhdr(fitsfile,"CRPIX3",ref_crpix[2]);
      if (ndim>3) fitwrhdr(fitsfile,"CRPIX4",ref_crpix[4]);
    } else {
      if (axistype==1) {
	fitwrhdr(fitsfile,"CRPIX1",xref[0]);      
	fitwrhdr(fitsfile,"CRPIX2",xref[1]);
	if (ndim>2) fitwrhdr(fitsfile,"CRPIX3",xref[2]);
	if (ndim>3) fitwrhdr(fitsfile,"CRPIX4",xref[3]);
      } else {
	fitwrhdr(fitsfile,"CRPIX1",1.0);        /* CRPIX = 1 by NEMO definition */
	fitwrhdr(fitsfile,"CRPIX2",1.0);
	if (ndim>2) fitwrhdr(fitsfile,"CRPIX3",1.0);
	if (ndim>3) fitwrhdr(fitsfile,"CRPIX4",1.0);
      }
    }
    if (Qrefmap || Qcrval) {
      dprintf(1,"Using ref_crval\n");
      fitwrhdr(fitsfile,"CRVAL1",ref_crval[0]);
      fitwrhdr(fitsfile,"CRVAL2",ref_crval[1]);
      if (ndim>2) {
	fitwrhdr(fitsfile,"CRVAL3",ref_crval[2]+vsys);
      }
      if (ndim>3) fitwrhdr(fitsfile,"CRVAL4",ref_crval[3]);
    } else {
      fitwrhdr(fitsfile,"CRVAL1",xmin[p[0]]);
      fitwrhdr(fitsfile,"CRVAL2",xmin[p[1]]);
      if (ndim>2) fitwrhdr(fitsfile,"CRVAL3",xmin[p[2]]);
      if (ndim>3) fitwrhdr(fitsfile,"CRVAL4",xmin[p[3]]);
    }

    if (Qcdmatrix) {
      fitwrhdr(fitsfile,"CD1_1",dx[p[0]]);    
      fitwrhdr(fitsfile,"CD2_2",dx[p[1]]);    
      if (ndim>2) fitwrhdr(fitsfile,"CD3_3",dx[p[2]]);    
    } else {
      if (Qrefmap || Qcdelt) {
	dprintf(1,"Using ref_cdelt\n");
	fitwrhdr(fitsfile,"CDELT1",ref_cdelt[0]*scale[0]);
	fitwrhdr(fitsfile,"CDELT2",ref_cdelt[1]*scale[1]);
	if (ndim>2) fitwrhdr(fitsfile,"CDELT3",ref_cdelt[2]*scale[2]);
	if (ndim>3) fitwrhdr(fitsfile,"CDELT4",1.0);
      } else {
	fitwrhdr(fitsfile,"CDELT1",-dx[p[0]]);    
	fitwrhdr(fitsfile,"CDELT2",dx[p[1]]);    
	if (ndim>2) {
	  if (Qfreq)
	    fitwrhdr(fitsfile,"CDELT3",-dx[p[2]]);
	  else
	    fitwrhdr(fitsfile,"CDELT3",dx[p[2]]);
	}
	if (ndim>3) fitwrhdr(fitsfile,"CDELT4",1.0);
      }
    }

    if (Qradecvel) {
      if (streq(proj,"SIN")) {
	ctype1_name = radeves[p[0]];
	ctype2_name = radeves[p[1]];
	ctype3_name = Qfreq ? radefrs[p[2]] : radeves[p[2]];
	ctype4_name = radeves[p[3]];
      } else if (streq(proj,"TAN")) {
	ctype1_name = radevet[p[0]];
	ctype2_name = radevet[p[1]];
	ctype3_name = Qfreq ? radefrt[p[2]] : radevet[p[2]];	
	ctype4_name = radevet[p[3]];
      } else
	error("Illegal projection scheme %s",proj);

      dprintf(0,"Axes names written as %s, %s, %s\n",
	      ctype1_name,ctype2_name,ctype3_name);
      fitwrhda(fitsfile,"CTYPE1",ctype1_name);
      fitwrhda(fitsfile,"CTYPE2",ctype2_name);
      if (ndim>2) fitwrhda(fitsfile,"CTYPE3",ctype3_name);
      fitwrhdr(fitsfile,"RESTFRQ",restfreq);  

      if (ndim>3) fitwrhda(fitsfile,"CTYPE4",ctype4_name);
      fitwrhdr(fitsfile,"RESTFREQ",restfreq);
      //
      fitwrhda(fitsfile,"CUNIT1","deg");
      fitwrhda(fitsfile,"CUNIT2","deg");
      if (Qfreq)
	fitwrhda(fitsfile,"CUNIT3","Hz");
      else
	fitwrhda(fitsfile,"CUNIT3","km/s");            /* or km/s */
      fitwrhda(fitsfile,"CUNIT4","");
      if (bmaj > 0.0)
	fitwrhda(fitsfile,"BUNIT","JY/BEAM");
      else
	fitwrhda(fitsfile,"BUNIT","JY/PIXEL");        /* if we have no beam */
      fitwrhda(fitsfile,"TELESCOP","NEMO");  
      
    } else {
      if (Qrefmap) {
        fitwrhda(fitsfile,"CTYPE1",ref_ctype[0]);
        fitwrhda(fitsfile,"CTYPE2",ref_ctype[1]);
        if (ndim>2) fitwrhda(fitsfile,"CTYPE3",ref_ctype[2]);
      } else {
	fitwrhda(fitsfile,"CTYPE1",axname[p[0]]);
	fitwrhda(fitsfile,"CTYPE2",axname[p[1]]);
	if (ndim>2) fitwrhda(fitsfile,"CTYPE3",axname[p[2]]);
	if (ndim>3) fitwrhda(fitsfile,"CTYPE4",axname[p[3]]);	
      }
    }

    fitwrhdr(fitsfile,"BMAJ",bmaj*scale[0]);
    fitwrhdr(fitsfile,"BMIN",bmin*scale[1]);
    fitwrhdr(fitsfile,"BPA",bpa);

    fitwrhdr(fitsfile,"DATAMIN",mapmin);
    fitwrhdr(fitsfile,"DATAMAX",mapmax);
    fitwrhda(fitsfile,"ORIGIN",origin);
    fitwrhda(fitsfile,"SPECSYS","LSRK");     /* spectral reference frame */
    

    cp = getenv("USER");                                /* AUTHOR */
    if (cp) {
        fitwrhda(fitsfile,"AUTHOR",cp);
        fitwrhda(fitsfile,"OBSERVER",cp);
    } else {
        fitwrhda(fitsfile,"AUTHOR","NEMO");
        fitwrhda(fitsfile,"OBSERVER","NEMO");
    }
    
    // DATE-OBS= '2015-06-14T15:42:02.319999'
    
    if (object)                                        /* OBJECT */
        fitwrhda(fitsfile,"OBJECT",object);

    if (comment)                                       /* COMMENT */
        stuffit(fitsfile,"COMMENT",comment);
    if (headline)
        stuffit(fitsfile,"COMMENT",headline);

    hitem = ask_history();                              /* HISTORY */
    fitwra(fitsfile,"HISTORY","NEMO: History in reversed order");
    for (i=0, cp=hitem[0]; cp != NULL; i++) {
    	stuffit(fitsfile,"HISTORY",cp);
        cp = hitem[i+1];
    }

    for(i=0; i<nfill; i++)   /* debugging header I/O */
        fitwra(fitsfile,"COMMENT","Dummy filler space");

    buffer = (float *) allocate(nx[p[0]]*sizeof(float));

    for (k=0; k<nx_out[2]; k++) {          /* loop over all planes */
        fitsetpl(fitsfile,1,&k);
        for (j=0; j<nx_out[1]; j++) {      /* loop over all rows */
 	  for (i=0, bp=buffer; i<nx_out[0]; i++, bp++) {
	    if (Qblank && CubeValue(iptr,i,j,k) == blankval)
	      *bp = fnan;
	    else
	      *bp =  iscale[0] * CubeValue(iptr,i,j,k) + iscale[1];
	  }
            fitwrite(fitsfile,j,buffer);
        }
    }
    free(buffer);
    fitclose(fitsfile);
}


/*

 stuff a character string accross the 80-line boundary
 NOTE: the exact implementation of this routines is still controversial

 From: pence@tetra.gsfc.nasa.gov Sat May 15, 1993 10:37 "Continuation Keywords"

 Option 1.   Blank Continuation Keyword Convention

 A continued string value is recognized by a backslash (\)
 as the last character in the first keyword string, followed
 by a keyword with a blank name followed by the continuation
 of the keyword value enclosed in quotes.  Example:

 KEYNAME = 'This is a very long keyword string value which \'
         'continues over several lines\'
                 ' of the FITS header.'

 */
                 
void stuffit(FITS *fitsfile, string fkey, string cp)
{
  char line[81], *hp;
  int i;
  int maxlen = 70;	/* comment field , minus 1 */
  
  hp = cp;
  strncpy(line,hp,maxlen);
  line[maxlen] = 0;
  fitwra(fitsfile,fkey,line);
  while ((int)strlen(hp) > maxlen) {
    hp += maxlen;
    strncpy(line,hp,maxlen);
    line[maxlen] = 0;
    fitwra(fitsfile," ",line);
  }
}

/* refmap stuff */

void set_refmap(string name)
{
  FITS *fitsfile;
  FLOAT tmpr, defval;
  int ndim = 3;
  int naxis[4], tmpi;
  int wcsaxes = -1;

  fitsfile = fitopen(name,"old",ndim,naxis);
  dprintf(0,"[Reading reference map %s [%d,%d,%d]\n",name,naxis[0],naxis[1],naxis[2]);

  /* set defaults according to Greisen & Calabretta 2002 WCS paper-I */
  defval = 1.0;
  fitrdhdr(fitsfile,"CDELT1",&ref_cdelt[0], defval);
  fitrdhdr(fitsfile,"CDELT2",&ref_cdelt[1], defval);
  fitrdhdr(fitsfile,"CDELT3",&ref_cdelt[2], defval);
  /* ieck; what if no CDELT's present, but all in CD matrix */

  defval = 0.0;
  fitrdhdr(fitsfile,"CRPIX1",&ref_crpix[0], defval);
  fitrdhdr(fitsfile,"CRPIX2",&ref_crpix[1], defval);
  fitrdhdr(fitsfile,"CRPIX3",&ref_crpix[2], defval);

  defval = 0.0;
  fitrdhdr(fitsfile,"CRVAL1",&ref_crval[0], defval);
  fitrdhdr(fitsfile,"CRVAL2",&ref_crval[1], defval);
  fitrdhdr(fitsfile,"CRVAL3",&ref_crval[2], defval);

  fitrdhda(fitsfile,"CTYPE1",ref_ctype[0],"");
  fitrdhda(fitsfile,"CTYPE2",ref_ctype[1],"");
  fitrdhda(fitsfile,"CTYPE3",ref_ctype[2],"");

  fitrdhda(fitsfile,"CUNIT1",ref_cunit[0],"");
  fitrdhda(fitsfile,"CUNIT2",ref_cunit[1],"");
  fitrdhda(fitsfile,"CUNUT3",ref_cunit[2],"");

  fitrdhdi(fitsfile,"WCSAXES",&wcsaxes, -1);
  if (wcsaxes != -1) warning("WCSAXES = %d\n",wcsaxes);

  fitclose(fitsfile);


}


/* return a reverse order index array  */
/* if done to an array 1 2 3, it will result in buggy data */
/* this routine should just shift the 1's to the end , not sort */

void permute (int *x ,int *p, int n)
{
  int    i, j, tmp;

  for (i=0; i<n; i++)
    p[i]=i;               /*  one-to-one */
                
  for (j=0; j<n; j++) {
    for (i=1; i<n; i++) {
      if (x[p[i-1]]==1 && x[p[i]] > 1) {
	tmp = p[i];
	p[i]   = p[i-1];
	p[i-1] = tmp;
      }
    }
  }
}


/* right now this is a silly bubble sort !!! */

void permute_old (int *x ,int *idx, int n)
{
    int    gap, i, j;
    int    tmp;
    for (i=0; i<n; i++)
        idx[i]=i;               /*  one-to-one */
                
    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (x[idx[j]] > x[idx[j+gap]])
                    break;          /* in order */
                tmp = idx[j];
                idx[j]    = idx[j+gap];
                idx[j+gap]= tmp;
            }
}
