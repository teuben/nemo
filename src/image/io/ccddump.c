/* 
 *	CCDDUMP: dump contents of an image files in some silly format
 *               Because item by item is output, it's bound to be a slow
 *               program - feel free to buffer...
 *		 Perhaps useful for ImageTool
 *
 *	29-jul-89  V1.0 created 	Peter Teuben
 *	23-sep-91  V1.1 added double-prec output mode		PJT
 *	21-may-92   1.1a  (SGI) fixed an ANSI complaint		pjt
 *	25-may-95 
 *	27-mar-97   1.1c  remove nested declarations            pjt
 *      28-mar-07   1.2   allow SM style                        pjt
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <image.h>

string defv[] = {	/* keywords + help string for user interface */
	"in=???\n	Input filename (image)",
	"out=???\n	Output filename (dump)",
	"option=float\n	Dump option [byte|float|double]",
	"range=\n	Range in case scaling is needed (A:B)",
        "mode=row\n	Output mode [row|column]",
        "swap=false\n	Swap bytes?[t|f]",
	"sm=f\n         use SM file_type mode CH?",
	"VERSION=1.2\n	28-mar-07 PJT",
	NULL,
};

string usage="dump the bytes of an image, optional scaling";

static bool Qsm;

nemo_main()
{
        imageptr iptr=NULL;
        stream instr, outstr;

	Qsm = getbparam("sm");

	instr = stropen (getparam("in"),"r");	/* get stream */
	read_image (instr,&iptr);               /* read image */
	strclose(instr);                        /* close image file */

	outstr = stropen(getparam("out"),"w");
	write_dump(outstr, iptr,
            getparam("mode"),getparam("option"),
            getparam("range"),getbparam("swap"));
	strclose(outstr);
}


write_dump(outstr,iptr,mode,option,range,swap)
stream   outstr;
imageptr iptr;
char     *mode, *option, *range;
bool     swap;
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
                error("Error: (%s) range A:B expected",range);
	    dprintf(0,"scaling.. %g : %g\n",omin,omax);
	} else {
	    omin = MapMin(iptr);
            omax = MapMax(iptr);
            dprintf(0,"autoscaling.. %g : %g\n",omin,omax);
	}
        write_dump_byte(outstr,iptr,omin,omax,outmode);
    } else if (*option == 'f') {                            /* float */
        write_dump_float(outstr,iptr,outmode,swap);
    } else if (*option == 'd') {                            /* double */
        write_dump_double(outstr,iptr,outmode,swap);
    } else
        error ("Unknown option %s",option);
}

write_dump_byte(outstr,iptr,omin,omax,outmode)
stream   outstr;
imageptr iptr;
real     omin, omax;
int      outmode;
{
    real z, scale, offset;
    int  ix, iy, nx, ny, nlen;
    char zout;

    nx = Nx(iptr);
    ny = Ny(iptr);
	dprintf(0,"write_dump_byte: range: %g:%g mode=%d map=%d * %d\n",
			omin, omax, outmode, nx, ny);
    scale = 255.0/(omax-omin);
    offset = omin;
    nlen = sizeof(char);
    if (outmode) {
        for (ix=0; ix<nx; ix++)
            for (iy=0; iy<ny; iy++) {
                z = MapValue(iptr,ix,iy);
		if (z<omin)
                    zout = 0;
                else if (z>omax)
                    zout = 255;
                else
                    zout = (char)  ((z - offset)*scale);
		fwrite(&zout,1,1,outstr);
            }
    } else {
        for (iy=0; iy<ny; iy++)
            for (ix=0; ix<nx; ix++) {
                z = MapValue(iptr,ix,iy);
		if (z<omin)
                    zout = 0;
                else if (z>omax)
                    zout = 255;
                else
                    zout = (char)  ((z - offset)*scale);
		fwrite(&zout,1,1,outstr);
            }
    }
}


write_dump_float(outstr,iptr,outmode,swap)
stream   outstr;
imageptr iptr;
int      outmode;
bool     swap;
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
            

write_dump_double(outstr,iptr,outmode,swap)
stream   outstr;
imageptr iptr;
int      outmode;
bool     swap;
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
    return(1);
}
    
