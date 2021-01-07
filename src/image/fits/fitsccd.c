/*
 *   FITSCCD:   read fits files - and near-fits files
 *
 *	ToDO:
 *		allow bscale/bzero/crpix etc. to be overriden
 *		even if fits file; i.e. make their
 *		defaults blanks, and use hasvalue() to detect mode.
 *              fix header if box= or planes= is used
 *             
 *
 *        ==>>  fix the NaN problem
 *
 *	 7-sep-89       started coding on 3b1
 *	28-jan-90	[32] V1.0 finished finally
 *       8-mar-90       V1.1 fixed bug on suns - added slice=
 *	16-mar-90       V2.0 New fits I/O routines, options=	PJT
 *	 1-oct-90	V3.0 Also allow Sault's fitsio routines PJT
 *	11-oct-90	V3.1 allow near-fits files PJT
 *       2-may-91       V3.2 fixed Xmin() bug       PJT
 *	10-sep-91       V3.3 allow blank value re-substitution	PJT
 *                       3.3a   also check datamin/max          PJT
 *	 7-mar-92	 3.4  gcc2.0 happy - new style		PJT
 *	25-aug-92	 3.5 fixed minmax bug if not in fits file PJT
 *      31-mar-93           c  fixed bug in fitrdhda() calling  PJT
 *	 9-sep-93           d  fixed bug - forgot bytepix in rawopen?
 *      29-jun-94           e  new fts_rdata() arg order (old code)   PJT
 *	26-jul-94	    g  allow bitpix=8 in rawopen	      PJT
 *	11-oct-95	    h  fixed bug in writing CTYPE3		pjt
 *	17-jan-99           i  nomore ny=1 or nx=1 checking needed      pjt
 *      18-may-99        3.6  read support for CD_i_j matrix            pjt
 *                            (should really also store the CD matrix)
 *	21-mar-00	    a never read the offset= properly in raw mode  pjt
 *       7-nov-00        3.7  allow relative coordinates read           pjt
 *      28-nov-00           a fixed un-important leak
 *      12-mar-01           b allow blank to be read from datamin/datamax
 *                          c 
 *                          d printf -> dprintf
 *      18-dec-01        4.0  convert to use fitsio_nemo.h and optional CFITSIO
 *      12-aug-02           a NaN problems with miriad
 *      14-apr-03           b fix NaN problems with FITS files
 *      23-nov-04        4.9  deal with axistype 1 images, but forced keyword   pjt
 *       3-dec-2013      5.0  showcs option      pjt
 *      18-feb-2015      5.1  add box=           pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <image.h>
#include <strlib.h>
#include <fitsio_nemo.h>

string defv[] = {
    "in=???\n		Input fits file",
    "out=\n		Output image file",
    "planes=\n          Planes to select from fits cube [all]", 
    "box=\n             Box (xmin,ymin,xmax,ymax) to select from fits cube [all]",
    "blocking=1\n	Extra blocking factor for input (blocksize/2880)",
    "mode=fits\n        Format mode of input file {fits,raw}",
    "blocksize=1\n      Blocksize in raw mode",
    "offset=0\n         Data offset in raw mode",
    "bitpix=32\n        Bitpix in raw mode {16,32,-32}",
    "naxis=\n           Length of the axes in raw mode",
    "cdelt=\n           Pixel separation in raw mode (2/naxis)",
    "crpix=0,0,0\n      Reference pixel in raw mode (0)",
    "crval=0,0,0\n      Coordinates of reference pixel in raw mode (0)",
    "bscale=1\n         Scale conversion factor in raw mode (1)",
    "bzero=0\n          Offset conversion factor in raw mode (0)",
    "blank=\n           Blank value re-substitution value?",
    "relcoords=f\n      Use relative (to crpix) coordinates instead abs",
    "axistype=1\n       Force axistype 0 (old, crpix==1) or 1 (new, crpix as is)",
    "VERSION=5.2\n	6-jan-2021 PJT",
    NULL,
};

string usage = "convert (near)fits files into ccd images";

string cvsid="$Id$";

#define MAXPLANES 2048

void make_fitheader(FITS *fitsfile, imageptr iptr, bool Qrel, bool Qout, int, FLOAT *, FLOAT *);
void make_rawheader(FITS *fitsfile, imageptr iptr, bool Qrel);
FITS *rawopen(string name, string status, int naxis, int *nsize);
void print_axis(int axis, int naxis, real crpix, real crval, real cdelt);
int is_feq(int *a, int *b);

void nemo_main()
{
    stream outstr;
    FITS *fitsfile;
    int ndim=3, naxis[3], nx, ny, nz, i, j, k, npl, p, planes[MAXPLANES];
    int nbox, box[4], i0, j0;
    int nbval=0;
    real bval_out, rmin, rmax, tmp, fbval;
    FLOAT *buffer, *bp, bval_in;  /* fitsio- is in FLOAT !!! */
    FLOAT fdata_min, fdata_max, fnan;
    int mir_nan = -1;    /* MIRIAD's FITS NaN */
    imageptr iptr;
    string mode, blankval;
    bool   Qblank, Qrel, Qout;
    int axistype = getiparam("axistype");

    Qout = hasvalue("out");
    if (Qout)
      outstr = stropen(getparam("out"),"w");   /* open image file for output */
    npl = nemoinpi(getparam("planes"),planes,MAXPLANES);
    if (npl<0) warning("Error parsing %d planes=",MAXPLANES);
    nbox = nemoinpi(getparam("box"),box,4);
    if (nbox<0) 
      error("Error parsing box=");
    else if (nbox == 0)
      box[0] = box[1] = box[2] = box[3] = 0;
    else if (nbox != 4) 
      error("Need 4 integers box=xmin,ymin,xmax,ymax");
    if (nbox>0) warning("New feature box= not fully tested");  /* > 2GP files wrong? */
    dprintf(1,"npl=%d nbox=%d\n",npl,nbox);
    mode = getparam("mode");
    blankval = getparam("blank");
    Qblank = (*blankval != 0);
    Qrel = getbparam("relcoords");
    get_nanf(&fnan);
#if 1
    memcpy(&fnan,&mir_nan,sizeof(int));    /* PORTABILITY ! */
#endif    
    if (streq(mode,"fits"))
      fitsfile = fitopen(getparam("in"),"old",ndim,naxis);
    else if (streq(mode,"raw"))
      fitsfile = rawopen(getparam("in"),"old",ndim,naxis);
    else
      error("Illegal mode %s: must be one of {fits,raw}",mode);

    if (fitsfile==NULL) error("Could not open file %s in %s mode\n",
               getparam("in"),mode);
    /* box=xmin,ymin,xmax,ymax */
    /*      0    1    2    3   */
    nx = (nbox>0 ? box[2]-box[0]+1 : naxis[0]);
    ny = (nbox>0 ? box[3]-box[1]+1 : naxis[1]);
    nz = (npl>0) ? npl : naxis[2];
    if (nx*ny*nz <=0) error("Bad fits image: nx*ny*nz=0");
    if (nx==1) warning("Fits image: nx=1");
    if (ny==1) warning("Fits image: ny=1");
    if (nz > 1) {
        dprintf(0,"FITS file: Image size %d %d %d\n",nx,ny,nz);
        create_cube(&iptr,nx,ny,nz);
    } else {
        dprintf(0,"FITS file: Image size %d %d\n",nx,ny);
        create_image(&iptr,nx,ny);
    }
    if (iptr==NULL) error("No memory to allocate image");

    if (streq(mode,"fits")) {
      make_fitheader(fitsfile,iptr,Qrel,Qout,axistype, &fdata_min, &fdata_max);
      dprintf(1,"Datamin/max read: %g - %g\n",fdata_min, fdata_max);
    } else if (streq(mode,"raw"))
      make_rawheader(fitsfile,iptr,Qrel);
    else
      error("Never Reached");

    if (!Qout) return;

    bval_out = 0.0;
    if (Qblank) {    
      if (streq(blankval,"datamin"))
	bval_in = fdata_min;
      else if (streq(blankval,"datamax"))
	bval_in = fdata_max;
      else if (streq(blankval,"mirnan")) {
	memcpy(&bval_in,&mir_nan,sizeof(int));    /* PORTABILITY ! */
      } else
	bval_in = (FLOAT) getdparam("blank");
      dprintf(0,"Substituting blank=%g [%s] with %g\n",
	      bval_in,blankval,bval_out);
    }


    buffer = (FLOAT *) allocate(naxis[0]*sizeof(FLOAT));

    rmin = HUGE;
    rmax = -HUGE;
    for (k=0; k<nz; k++) {          /* loop over all/selected planes */
      p = (npl>0) ? planes[k] : k;        /* select plane number */
      dprintf(2,"Reading plane %d\n",p);
      fitsetpl(fitsfile,1,&p);
      for (j=0; j<ny; j++) {      /* loop over all rows */
	j0 = (nbox == 0 ? j : j+box[1]-1);
	fitread(fitsfile,j0,buffer);     /* read it from fits file */
	i0 = (nbox == 0 ? 0 : box[0]-1);
	for (i=0, bp=&buffer[i0]; i<nx; i++, bp++) {   /* stuff it in memory */
	  if (Qblank) {
	    if (is_feq((int *)bp,(int *)&bval_in)) {
	      nbval++;
	      *bp = bval_out;
	    } else if (isnan(*bp)) {
	      nbval++;
	      *bp = bval_out;
	    }
	  } else {
	    if (is_feq((int *)bp,(int *)&fnan)) {
	      nbval++;
	      *bp = bval_out;
	    } else if (isnan(*bp)) {
	      nbval++;
	      *bp = bval_out;
	    }
	    dprintf(2,"%g %g %g %g, %d\n",*bp,fdata_min, fdata_max,fnan,nbval); 
	  }
	  tmp = CubeValue(iptr,i,j,k) = *bp;
	  rmin=MIN(rmin,tmp);
	  rmax=MAX(rmax,tmp);
	}
      }
    }
    if (rmin != MapMin(iptr)) {
      warning("Setting map minimum from %g to %g",MapMin(iptr),rmin);
      MapMin(iptr) = rmin;
    } 
    if (rmax != MapMax(iptr)) {
      warning("Setting map maximum from %g to %g",MapMax(iptr),rmax);
      MapMax(iptr) = rmax;
    }
    fitclose(fitsfile);
    write_image(outstr,iptr);
    if (nbval==0)
      dprintf(0,"There were no blank values set in the image\n");
    else {
      fbval = (1.0*nbval)/(nx*ny*nz);
      dprintf(0,"There were %d blank values in the image (%f %%)\n",nbval,fbval*100);
    }
    free(buffer);
}

#ifdef HAVE_LIBCFITSIO
FITS *rawopen(string name, string status, int naxis, int *nsize)
{
  error("raw I/O not implemented for CFITSIO");
}
void make_rawheader(FITS *fitsfile, imageptr iptr, bool Qrel)
{
  error("raw I/O not implemented for CFITSIO");
}
#else
FITS *rawopen(string name, string status, int naxis, int *nsize)
{
    FITS *f;
    int i,ndim,bitpix;
    
    f = (FITS *) allocate(sizeof(FITS));

    if (streq(status,"old")){       /*  handle an old file */
        f->fd = stropen(name,"r");
        ndim = nemoinpi(getparam("naxis"),nsize,naxis);
        if (ndim < 1) error("naxis= returns %d from parser",ndim);
        f->naxis = ndim;
        for (i=ndim; i<naxis; i++) nsize[i] = 1;
        for (i=0; i<ndim; i++)  f->axes[i] = nsize[i];
        for (i=ndim; i<MAXNAX; i++) f->axes[i] = 1;
        bitpix = getiparam("bitpix");
        fit_setbitpix(bitpix);
        if (bitpix==8)        f->type = TYPE_8INT;
        else if (bitpix==16)  f->type = TYPE_16INT;
        else if (bitpix==32)  f->type = TYPE_32INT;
        else if (bitpix==-32) f->type = TYPE_FLOAT;
        else error("Invalid bitpix: %d",bitpix);
        f->bytepix = ABS(bitpix)/8;
        f->status = STATUS_OLD;
        f->bscale = getdparam("bscale");
        f->bzero  = getdparam("bzero");
        f->offset = getiparam("offset");
        f->skip = f->offset;
        f->ncards = 0;
    } else
        error("Cannot rawopen %s as \"%s\"\n",name,status);
    fit_setblocksize( getiparam("blocksize") * getiparam("blocking") );
    return f;
}

void make_rawheader(FITS *fitsfile, imageptr iptr, bool Qrel)
{
    int ndim, n;
    double tmpr[3];
    
    ndim = (Nz(iptr) > 1) ? 3 : 2;

    n = nemoinpd(getparam("cdelt"),tmpr,3);
    if (n == 0) {
        Dx(iptr) = 2.0/Nx(iptr);
        Dy(iptr) = 2.0/Ny(iptr);
        Dz(iptr) = 2.0/Nz(iptr);
    } else if (n == 3) {
        Dx(iptr) = tmpr[0];
        Dy(iptr) = tmpr[1];
        Dz(iptr) = tmpr[2];
    } else
        error("cdelt= needs either 3 or 0 entries");

    Xmin(iptr) = 0.0;
    Ymin(iptr) = 0.0;
    Zmin(iptr) = 0.0;

}
#endif

void print_axis(int axis, int naxis, real crpix, real crval, real cdelt)
{
  real xmin = (1.5-crpix)*cdelt;
  real xmax = (naxis+0.5-crpix)*cdelt;
  printf("AXIS%d: %d %g %g %g   %g %g\n", axis,naxis,crpix,crval,cdelt,xmin,xmax);
}


/* 
 * @todo   if nbox>0 header needs adjusting
 *         crval1 -= box[0]    *         crval2 -= box[1]
 */


void make_fitheader(FITS *fitsfile, imageptr iptr, bool Qrel, bool Qout, int axistype,
		    FLOAT *data_min, FLOAT *data_max)
{
    int nz, tmpi, i, j;
    real crpix1, crpix2, crpix3;
    FLOAT tmpr, cd[3][3];
    char cdname[10], ctype[32];

    nz = Nz(iptr);

    fitrdhdr(fitsfile,"CRVAL1",&tmpr,0.0); Xmin(iptr) = tmpr;
    fitrdhdr(fitsfile,"CRVAL2",&tmpr,0.0); Ymin(iptr) = tmpr;
    fitrdhdr(fitsfile,"CRVAL3",&tmpr,0.0); Zmin(iptr) = tmpr;
    if (Qrel) {
      Xmin(iptr) = Ymin(iptr) = Zmin(iptr) = 0.0;
    }

    fitrdhdr(fitsfile,"CDELT1",&tmpr,1.0); Dx(iptr)=tmpr;
    fitrdhdr(fitsfile,"CDELT2",&tmpr,1.0); Dy(iptr)=tmpr;
    fitrdhdr(fitsfile,"CDELT3",&tmpr,1.0); Dz(iptr)=tmpr;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	sprintf(cdname,"CD%d_%d",i,j);
	fitrdhdr(fitsfile,cdname,&tmpr, (i==j ? 1.0 : 0.0));
	dprintf(1,"%s: %g\n",cdname,tmpr); 
	cd[i-1][j-1] = tmpr;
      }
    if (fitexhd(fitsfile,"CD1_1")) {
      if (fitexhd(fitsfile,"CDELT1"))
	  warning("CDELT1 as well as CD_1_1 have been specified - using CD");
      Dx(iptr) = cd[0][0];      
      Dy(iptr) = cd[1][1];
      Dz(iptr) = cd[2][2];
    }
	
    fitrdhdr(fitsfile,"CRPIX1",&tmpr,1.0); crpix1 = tmpr;
    fitrdhdr(fitsfile,"CRPIX2",&tmpr,1.0); crpix2 = tmpr;
    if (nz>1) {
      fitrdhdr(fitsfile,"CRPIX3",&tmpr,1.0); crpix3 = tmpr;
    } else
      crpix3 = 0.0;

    if (!Qout) {
      print_axis(1, Nx(iptr), crpix1, Xmin(iptr), Dx(iptr));
      print_axis(2, Ny(iptr), crpix2, Ymin(iptr), Dy(iptr));
      print_axis(3, Nz(iptr), crpix3, Zmin(iptr), Dz(iptr));
      return;
    }

    if (axistype==0) {
      Axis(iptr) = 0;
      if (crpix1 != 1.0)
        Xmin(iptr) -= (crpix1-1.0)*Dx(iptr);
      if (crpix2 != 1.0)
        Ymin(iptr) -= (crpix2-1.0)*Dy(iptr);
      if (nz>1 && crpix3 != 1.0)
        Zmin(iptr) -= (crpix3-1.0)*Dz(iptr);
    } else if (axistype==1) {
      Axis(iptr) = 1;
      Xref(iptr) = crpix1-1;
      Yref(iptr) = crpix2-1;
      Zref(iptr) = crpix3-1;
    } else
      error("Illegal axistype=%d",axistype);

    fitrdhda(fitsfile,"CTYPE1",ctype,"");
    Namex(iptr) = scopy(ctype);

    fitrdhda(fitsfile,"CTYPE2",ctype,"");
    Namey(iptr) = scopy(ctype);

    if (nz>1) {
           fitrdhda(fitsfile,"CTYPE3",ctype,"");
	   Namez(iptr) = scopy(ctype);
    }


    fitrdhdr(fitsfile,"DATAMIN",&tmpr,0.0); 
    MapMin(iptr) = *data_min = tmpr;
    dprintf(1,"DATAMIN: %g %g\n",*data_min,MapMin(iptr));
    fitrdhdr(fitsfile,"DATAMAX",&tmpr,0.0); 
    MapMax(iptr) = *data_max = tmpr;
    dprintf(1,"DATAMAX: %g %g\n",*data_max,MapMax(iptr));
    if (fitexhd(fitsfile,"BLANK")) {
        dprintf(1,"BLANK keyword exists\n");
        fitrdhdi(fitsfile,"BITPIX",&tmpi,-1);
        if (tmpi<0) {
            warning("FITS keyword BLANK not interpreted - BITPIX<0");
        } else {
            fitrdhdi(fitsfile,"BLANK",&tmpi,0);
            dprintf(1,"BLANK = %d\n",tmpi);
        }
    }
}


/*
 *  the _nan routines should supply this
 */

is_feq(int *ia, int *ib)
{
  int r;
#if 0
  r = memcmp(&a,&b,sizeof(float));
  int ia, ib;
  memcpy(&ia,&a,sizeof(float));
  memcpy(&ib,&b,sizeof(float));
  printf("A = %d    B=%d   \n",ia,ib);
#else
  r = (*ia == *ib);
#endif
  return r;
}
