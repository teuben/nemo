/* 
 *	CCDDUMP: dump contents of an image files in some simple dump format
 *               Because item by item is output, it's bound to be slow
 *
 *	29-jul-89  V1.0 created 	Peter Teuben
 *	23-sep-91  V1.1 added double-prec output mode		PJT
 *	21-may-92   1.1a  (SGI) fixed an ANSI complaint		pjt
 *	25-may-95 
 *	27-mar-97   1.1c  remove nested declarations            pjt
 *      28-mar-07   1.2   allow SM style                        pjt
 *       9-jan-12   1.3   ds9 style scale options, add 3D       pjt
 *                        lots of unfinished options
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <image.h>

string defv[] = {	/* keywords + help string for user interface */
  "in=???\n	    Input filename (image)",
  "out=???\n	    Output filename (dump)",
  "option=float\n   Dump option [byte|short|int|float|double]",
  "range=\n	    Range in case scaling is needed (A:B)",
  "mode=row\n       Output mode [row|column]",
  "swap=false\n	    Swap bytes?[t|f]",
  "scale=linear\n   ds9 style scaling [linear|log|pow|sqrt|square|gamma|mid|asinh]",
  "gamma=2\n        gamma factor in gamma and mid scaling",
  "b=1\n            softening parameter for the asinh magnitude scale",
  "sm=f\n           use SM file_type mode CH?",
  "VERSION=1.3\n    9-jan-2012 PJT",
  NULL,
};

string usage="dump the bytes of an image, optional scaling";

static bool Qsm;
static real p_gamma;
static real p_a = 1.08574;   /* 2.5 log(e)  ; Pogson 1856 */
static real p_b;

static rproc my_scale;

/* scaling functions that should stay with [0,1] on input and output */

real my_linear(real x)  {   return x; }
real my_sqrt(real x)    {   return sqrt(x); }
real my_sqr(real x)     {   return x*x; }
real my_pow(real x)     {   return pow(x,10); }
real my_log(real x)     {   return log10(x); }
real my_gamma(real x)   {   
  dprintf(2,"my_gamma %g\n",x);
  return pow(x,p_gamma); }
real my_mid(real x)     { 
  if (x<0.5) return 0.5*(1-pow(1-2*x,p_gamma));
  return 0.5*(1+pow(2*x-1,p_gamma));
}
real my_asinh(real x)     { 
  return x;
}

/* 'asinh' see also: 1999AJ....118.1406L  */

int write_dump(stream outstr, imageptr iptr, char *mode, char *option, char *range, int swap);
int write_dump_byte(stream outstr, imageptr iptr, real omin, real omax, int outmode);
int write_dump_float(stream outstr, imageptr iptr, int outmode, int swap);
int write_dump_double(stream outstr, imageptr iptr, int outmode, int swap);
int get_range(char *s, real *a, real *b);

nemo_main()
{
  imageptr iptr=NULL;
  stream instr, outstr;
  string scale;

  Qsm = getbparam("sm");

  instr = stropen (getparam("in"),"r");	/* get stream */
  read_image (instr,&iptr);               /* read image */
  strclose(instr);                        /* close image file */

  p_gamma = getrparam("gamma");
  p_b     = getrparam("b");
  scale = getparam("scale");
  if (streq(scale,"linear")) {
    dprintf(0,"linear scaling\n");
    my_scale = my_linear;
  } else if (streq(scale,"sqrt")) {
    dprintf(0,"sqrt scaling\n");
    my_scale = my_sqrt;
  } else if (streq(scale,"square")) {
    dprintf(0,"square scaling\n");
    my_scale = my_sqr;
  } else if (streq(scale,"mid")) {
    dprintf(0,"mid scaling, gamma=%g\n",p_gamma);
    my_scale = my_mid;
  } else if (streq(scale,"gamma")) {
    dprintf(0,"gamma scaling, gamma=%g\n",p_gamma);
    my_scale = my_gamma;
  } else if (streq(scale,"asinh")) {
    dprintf(0,"asinh scaling, b=%g\n",p_b);
    my_scale = my_asinh;
  } else
    error("unknown scale=%s",scale);
  
  outstr = stropen(getparam("out"),"w");
  write_dump(outstr, iptr,
	     getparam("mode"),getparam("option"),
	     getparam("range"),getbparam("swap"));
  strclose(outstr);
}

int write_dump(stream outstr, imageptr iptr, 
	       char *mode, char *option, char *range, int swap) 
{
  int  outmode=0;
  real omin,omax;

#if 0
  if (strncmp(Storage(iptr),"CDef",4)==0) {              /* a kludge */
    outmode = 0;
    if (*mode=='c') outmode=1;
  } else if (strncmp(Storage(iptr),"ForDef",6)==0) {
    outmode = 1;
    if (*mode=='r') outmode=0;
  } else
    dprintf(0,"Warning: Unknown storage scheme in image file: %s\n",
	    Storage(iptr));
#endif

  if (*option == 'b') {                                   /* byte */
    if (*range) {
      if (!get_range(range,&omin,&omax))
	error("Error: (%s) range=min:max expected",range);
      dprintf(0,"scaling.. %g : %g\n",omin,omax);
    } else {
      omin = MapMin(iptr);
      omax = MapMax(iptr);
      dprintf(0,"autoscaling.. %g : %g\n",omin,omax);
    }
    write_dump_byte(outstr,iptr,omin,omax,outmode);
  } else if (*option == 's') {                            /* short */
    error("16 byte not implemented yet");
  } else if (*option == 'i') {                            /* int */
    error("32 byte not implemented yet");
  } else if (*option == 'f') {                            /* float */
    write_dump_float(outstr,iptr,outmode,swap);
  } else if (*option == 'd') {                            /* double */
    write_dump_double(outstr,iptr,outmode,swap);
  } else
    error ("Unknown option %s",option);
}

int write_dump_byte(stream outstr, imageptr iptr, 
		    real omin, real omax, int outmode) 
{
  real z, scale, offset;
  int  ix, iy, iz, nx, ny, nz, nlen;
  char zout;
  
  nx = Nx(iptr);
  ny = Ny(iptr);
  nz = Nz(iptr);
  dprintf(0,"write_dump_byte: range: %g:%g mode=%d map=%d * %d * %d\n",
	  omin, omax, outmode, nx, ny, nz);
  scale = 1.0/(omax-omin);
  offset = omin;
  nlen = sizeof(char);
  if (outmode) {
    /* @todo: fix for 3D */
    for (ix=0; ix<nx; ix++)
      for (iy=0; iy<ny; iy++) {
	z = MapValue(iptr,ix,iy);
	if (z<=omin)
	  zout = 0;
	else if (z>=omax)
	  zout = (char) 255;
	else {
	  z = (z - offset)*scale;
	  z = (*my_scale)(z);
	  zout = (char) (z*256); 
	}
	fwrite(&zout,1,1,outstr);
      }
  } else {
    /* this is the default way */
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++) {
	  z = CubeValue(iptr,ix,iy,iz);
	  if (z<=omin)
	    zout = 0;
	  else if (z>=omax)
	    zout = (char) 255;
	  else {
	    z = (z - offset)*scale;
	    z = (*my_scale)(z);
	    zout = (char) (z*256);  /* should never be out of 0..255 */
	  }
	  fwrite(&zout,1,1,outstr);
	}
  }
}


int write_dump_float(stream outstr, imageptr iptr, 
		     int outmode, int swap)

{
  int   ix, iy, iz, nx, ny, nz, nlen, itmp;
  float zout, ftmp;

  nx = Nx(iptr);
  ny = Ny(iptr);
  nz = Nz(iptr);
  
  if (Qsm) {
    if (nz>1) error("3D images not supported in SM mode");
    warning("testing new SM mode, not swapping bytes though");

    /* this is the CH filetype in sm */
    
    itmp=nx;              if(swap)bswapi(&itmp,1);    fwrite(&itmp,sizeof(int),1,outstr);
    itmp=ny;              if(swap)bswapi(&itmp,1);    fwrite(&itmp,sizeof(int),1,outstr);
    ftmp = Xmin(iptr);    if(swap)bswapf(&ftmp,1);    fwrite(&ftmp,sizeof(float),1,outstr);
    ftmp += nx*Dx(iptr);  if(swap)bswapf(&ftmp,1);    fwrite(&ftmp,sizeof(float),1,outstr);
    ftmp = Ymin(iptr);    if(swap)bswapf(&ftmp,1);    fwrite(&ftmp,sizeof(float),1,outstr);
    ftmp += ny*Dy(iptr);  if(swap)bswapf(&ftmp,1);    fwrite(&ftmp,sizeof(float),1,outstr);
  }
  nlen = sizeof(float);
  if (outmode) {          /* column by column; plane after plane */
    for (iz=0; iz<nz; iz++)
      for (ix=0; ix<nx; ix++)
	for (iy=0; iy<ny; iy++) {
	  zout = (float) CubeValue(iptr,ix,iy,iz);
	  if (swap) bswapf(&zout,1);
	  fwrite(&zout,nlen,1,outstr);
	}
  } else {                    /* row by row; plane after plane */
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++) {
	  zout = (float) CubeValue(iptr,ix,iy,iz);
	  if (swap) bswapf(&zout,1);
	  fwrite(&zout,nlen,1,outstr);
	}
  }
}
            
int write_dump_double(stream outstr, imageptr iptr, 
		      int outmode, int swap)
{
    int   ix, iy, iz, nx, ny, nz, nlen;
    double zout;

    nx = Nx(iptr);
    ny = Ny(iptr);
    nz = Nz(iptr);
    nlen = sizeof(double);
    if (outmode) {          /* column by column; plane after plane */
      for (iz=0; iz<nz; iz++)
        for (ix=0; ix<nx; ix++)
            for (iy=0; iy<ny; iy++) {
                zout = (double) CubeValue(iptr,ix,iy,iz);
                if (swap) bswapd(&zout,1);
                fwrite(&zout,nlen,1,outstr);
            }
    } else {                    /* row by row; plane after plane */
      for (iz=0; iz<nz; iz++)
        for (iy=0; iy<ny; iy++)
            for (ix=0; ix<nx; ix++) {
                zout = (double) CubeValue(iptr,ix,iy,iz);
                if (swap) bswapd(&zout,1);
                fwrite(&zout,nlen,1,outstr);
            }
    }
}
            
    
get_range(char *s, real *a, real *b)
{
  char *cp;

  if (*s==0) {
    dprintf(0,"error: (empty string) range 'A:B' expected\n");
    return 0;
  }
  cp = s;
  *a = atof(cp);
  cp = strchr(cp,':');
  if (*cp++==0) {
    dprintf(0,"error: (%s) range 'A:B' expected\n",s);
    return 0;
  }
  *b = atof(cp);
  return 1;
}

