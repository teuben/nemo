/* 
 * CCDSMOOTH: smooth image cube in XYZ
 *
 *	25-jun-87  V1.0 adapted from old CCD, handles new image format only   PJT
 *	30-jun-87  V1.1 handles new filestructure with 'struct'	PJT
 *	 9-jul-87  V1.2 keyword order, added 'dir' keyword 	PJT
 *	 1-jun-88  V2.0 new filestruct; although source code is same	PJT
 *	 7-jan-89  V2.0a: less output 
 *	18-jan-89  V2.1: a little compatible with new 3D cubes PJT
 *	19-feb-89  V2.2a: bug removed - beam size when smoothing in Y only PJT
 *       1-nov-90  V2.3  for new NEMO       PJT
 *	 7-mar-92  V2.4  modern NEMO style - gcc2.0 happy          PJT
 *      13-feb-97  V3.0  added bad= to ignore bad values	   PJT
 *	 7-mar-98      a attempted to handle maps with < 0 Dx/Dy/Dz	PJT
 *      12-mar-98  V3.1  handle gauss=0 and added cut= keyword          PJT
 *	20-apr-01      a bigger default size for MSIZE			pjt
 *      30-jun-2016 V3.4 option to use a moffat smoothing
 *
 *	"Smoothing is art, not science"
 *				- Numerical Recipies, p495
 *
 *  todo
 *     - near an edge, normalization is done with full beam, so the signal tapers off....
 *
 *  moffat:   https://ui.adsabs.harvard.edu/abs/1969A%26A.....3..455M/abstract
 */

#include <stdinc.h>
#include <string.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n               Input filename",
	"out=???\n              Output filename",
	"gauss=\n               FWHP gaussian beamwidth, if used",
	"moffat=\n              FWHM moffat function, if used",
	"dir=xy\n               Coordinates to smooth (xyz)",			
	"smooth=0.25,0.5,0.25\n Alternate smoothing array if gauss= not used",
	"wiener=\n              Size of a wiener filter, if used (square in X,Y)",
	"nsmooth=1\n            Number of smoothings",
	"bad=\n			Optional ignoring this bad value",
	"beta=4.765\n           beta parameter for the moffat function",
	"cut=0.01\n             Cutoff value for gaussian, if used",
	"mode=0\n               Special edge smoothing modes (testing)",
	"VERSION=3.4\n          29-jun-2016 PJT",
	NULL,
};

string usage = "smooth image cube in XYZ";

#ifndef HUGE
# define HUGE 1.0e20
#endif

string	infile, outfile;			/* file names */
stream  instr, outstr;				/* file streams */

#define MSIZE  8196		      /* maximum # pixels along one dimension */
#define MSMOOTH 101 		    /* maximum full beam-size (has to be odd) */
	              /* because of symmetry, you could try and be smart here */

imageptr iptr=NULL;			/* will be allocated dynamically */
int    nx,ny,nz,nsize;			/* actual size of map */
real   xmin,ymin,zmin,dx,dy,dz;
real   size;				/* size of frame (square) */
real   cell;				/* cell or pixel size (square) */
int    nxw=0, nyw=0, nzw=0;             /* wiener filter size */

real   smooth[MSMOOTH];			/* full 1D beam */
int    lsmooth;				/* actual smoothing length */
int    nsmooth;				/* number of smoothings */
real   gauss_fwhm;			/* beam size in case gauss */
real   gauss_cut;                       /* cutoff of gaussian */
real   moffat_fwhm;
real   moffat_alpha;
real   moffat_beta;
string dir;                             /* direction to smooth 'xyz' */
int    mode;                            /* special edge smoothing modes */

bool   Qbad;                            /* ignore smoothing for */
real   bad;                             /* this value */

void setparams(), smooth_it(), wiener();
int convolve_cube (), convolve_x(), convolve_y(), convolve_z();

void make_gauss_beam(char *sdir);
void make_moffat_beam(char *sdir);


void nemo_main(void)
{
    setparams();			/* set par's in globals */

    instr = stropen(infile, "r");
    read_image(instr,&iptr);
				/* set some global paramters */
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    xmin = Xmin(iptr);
    ymin = Ymin(iptr);
    zmin = Zmin(iptr);
    dx = Dx(iptr);
    dy = Dy(iptr);
    dz = Dz(iptr);

    if(hasvalue("gauss"))
        make_gauss_beam(getparam("dir"));
    else if (hasvalue("moffat")) {  
	moffat_fwhm = getrparam("moffat");
        moffat_beta = getrparam("beta");
	moffat_alpha = moffat_fwhm / (2*sqrt(pow(2.0,1.0/moffat_beta)-1));
	dprintf(0,"Moffat: alpha=%g beta=%g\n",moffat_alpha, moffat_beta);
	make_moffat_beam(getparam("dir"));
    }

    outstr = stropen(outfile,"w");
    if (hasvalue("wiener"))
      wiener();
    else
      smooth_it();
    write_image(outstr,iptr);

    strclose(instr);
    strclose(outstr);
}

void setparams()
{
    int     nw,nws[3];
    infile = getparam ("in");
    outfile = getparam ("out");

    if(hasvalue("gauss")) {         /* gaussian beam */
        if(nemoinpr(getparam("gauss"),&gauss_fwhm,1) != 1)
            error ("gauss=%s : beam must have 1 dimension",
                    getparam("gauss"));
        gauss_cut = getrparam("cut");
    } else              		/* beam by hand */
        lsmooth = nemoinpr(getparam("smooth"),smooth,MSMOOTH);
    nsmooth = getdparam("nsmooth");
    dir = getparam("dir");
    Qbad = hasvalue("bad");
    if (Qbad) bad = getdparam("bad");
    mode = getiparam("mode");
    nw = nemoinpi("wiener",nws,3);
    if (nw>0) {
      /* for now */
      nxw = nyw = nws[0];
    }
}

void smooth_it()
{
    real m_min, m_max, brightness, total;
    int    i, ix, iy, iz, kounter, idir;
    char   *cp;

    m_min = HUGE;
    m_max = -HUGE;
    total = 0.0;

    kounter = nsmooth;
    while (kounter-- > 0) {
     	dprintf (1,"Convolving %s with %d-length beam: ",dir,lsmooth);
	for (i=0; i<lsmooth; i++)
		dprintf (1," %f ",smooth[i]);
	cp = dir;			/* point to direction again */
        while (*cp) {
            if (*cp=='x')
	        idir=1;
	    else if (*cp=='y')
                idir=2;
            else if (*cp=='z')
                idir=3;
            else
	        error("Wrong direction %c for beamsmoothing\n",*cp);
	    convolve_cube (Frame(iptr),nx,ny,nz,smooth,lsmooth,idir);
            cp++;
	}
    }

    m_max = -HUGE;                      /* determine new min/max */
    m_min =  HUGE;
    for (ix=0; ix<Nx(iptr); ix++)   	
    for (iy=0; iy<Ny(iptr); iy++)
    for (iz=0; iz<Nz(iptr); iz++) {
          brightness = CubeValue(iptr,ix,iy,iz);
	  total += brightness;
          m_max = MAX(m_max, brightness);
          m_min = MIN(m_min, brightness);
    }
    MapMin(iptr) = m_min; 		/* update map headers */
    MapMax(iptr) = m_max;
    if (hasvalue("gauss")) {
    	BeamType(iptr)=GAUSS;
	if (strchr(getparam("dir"),'x'))
    	   Beamx(iptr)=gauss_fwhm;
 	if (strchr(getparam("dir"),'y'))
    	   Beamy(iptr)=gauss_fwhm;
	if (strchr(getparam("dir"),'z'))
    	   Beamz(iptr)=gauss_fwhm;
     } else {
    	BeamType(iptr)=ANYBEAM;
	/*  factor (3+nsmooth)/2   is very roughly ok for nsmooth=1,2,3,4 */
	Beamx(iptr) = (1.5+0.5*nsmooth) * ABS(Dx(iptr));
	Beamy(iptr) = (1.5+0.5*nsmooth) * ABS(Dy(iptr));
	Beamz(iptr) = (1.5+0.5*nsmooth) * ABS(Dz(iptr));
     }
    	
    dprintf (1,"New min and max in map are: %f %f\n",m_min,m_max);
    dprintf (1,"New total brightness/mass is %f\n",total*Dx(iptr)*Dy(iptr));
}


/*  sdir:    smoothing direction: 'x', 'y' or 'z' */

void make_gauss_beam(char *sdir)

{
    real sigma, sum, x, fx;
    int    i;

    sigma = gauss_fwhm/2.355;   /* FWHM = 2 sqrt(2 ln(2)) * sigma */
    lsmooth=1;              /* count how many we need for the 'smooth' ARRAY */
    x=0;                    /* minimum is 1, at x=0 */
    sum=1.0;                /* unnormalized value at x=0 */
    if (*sdir == 'x')
        cell = ABS(dx);
    else if (*sdir == 'y')
        cell = ABS(dy);
    else if (*sdir == 'z')
        cell = ABS(dz);
    if (cell==0.0)    
        error("make_gauss_beam: cellsize zero for dir=%c",*sdir);
    dprintf(0,"dir=%c sigma=%g cell=%g\n",*sdir,sigma,cell);

    if (gauss_fwhm == 0.0) {
        warning("No smoothing done in %c",*sdir);
        lsmooth = 1;
        smooth[0] = 1.0;
        return;
    }

    do {
        x += cell;
        lsmooth += 2;
        fx = exp(-0.5*sqr(x/sigma));            
        sum += 2*fx;
    } while ((fx>gauss_cut) && (lsmooth<MSMOOTH));  /* cutoff beam */
    sum = 1.0/sum;
    dprintf (1,"Gaussean beam will contain %d points at cutoff %g \n",
                lsmooth, gauss_cut);
    if ((lsmooth==MSMOOTH) && (fx>gauss_cut)) 
        dprintf(0,"Warning: beam cutoff %g at %d points\n",
                gauss_cut, MSMOOTH);
    for (i=0; i<lsmooth; i++) {
        x = (i-(lsmooth-1)/2)*cell/sigma;
        smooth[i] = exp(-0.5*x*x) * sum;   /* normalize */
    }
}

void make_moffat_beam(char *sdir)
{
    real sigma, sum, x, fx;
    int    i;

    sigma = moffat_fwhm;    /* FWHM = 2*a*sqrt(2^(1/b)-1) */
    lsmooth=1;              /* count how many we need for the 'smooth' ARRAY */
    x=0;                    /* minimum is 1, at x=0 */
    sum=1.0;                /* unnormalized value at x=0 */
    if (*sdir == 'x')
        cell = ABS(dx);
    else if (*sdir == 'y')
        cell = ABS(dy);
    else if (*sdir == 'z')
        cell = ABS(dz);
    if (cell==0.0)    
        error("make_moffat_beam: cellsize zero for dir=%c",*sdir);
    dprintf(0,"dir=%c sigma=%g cell=%g\n",*sdir,sigma,cell);

    do {
        x += cell;
        lsmooth += 2;
        fx = 1.0 / pow( 1+ sqr(x/moffat_alpha), moffat_beta);
        sum += 2*fx;
    } while ((fx>gauss_cut) && (lsmooth<MSMOOTH));  /* cutoff beam */
    sum = 1.0/sum;
    dprintf (1,"Moffat beam will contain %d points at cutoff %g \n",
                lsmooth, gauss_cut);
    if ((lsmooth==MSMOOTH) && (fx>gauss_cut)) 
        dprintf(0,"Warning: beam cutoff %g at %d points\n",
                gauss_cut, MSMOOTH);
    for (i=0; i<lsmooth; i++) {
        x = (i-(lsmooth-1)/2)*cell;
        smooth[i] = sum / pow(1+sqr(x/moffat_alpha),moffat_beta);
    }
}
                

int convolve_cube (a, nx, ny, nz, b, nb, idir)
real *a, b[];
int  nx,ny,nz,nb,idir;
{
    int ix,iy,iz, ier;
    
    if (idir==1)
        for (iy=0; iy<ny; iy++)
        for (iz=0; iz<nz; iz++)
            ier += convolve_x (a,iy,iz,nx,ny,nz,b,nb);
    else if (idir==2)
        for (ix=0; ix<nx; ix++)
        for (iz=0; iz<nz; iz++)
            ier += convolve_y (a,ix,iz,nx,ny,nz,b,nb);
    else if (idir==3)
        for (ix=0; ix<nx; ix++)
        for (iy=0; iy<ny; iy++)
            ier += convolve_z (a,ix,iy,nx,ny,nz,b,nb);
    else
        return 0;
    return ( ier==0 ? 1 : 0 );
}


int convolve_x (a, iy, iz, nx, ny, nz, b, nb)
real *a, b[];
int    nx, ny, nz, nb, iy, iz;
{
	real c[MSIZE];
	int    ix, jx, kx, offset;
	
	if (nx>MSIZE) {
		warning("convolve_x: MSIZE=%d array too small for input %d\n",
			MSIZE,nx);
                return 0;
	}
	
	for (ix=0; ix<nx; ix++) {
		offset = iz + iy*nz + ix*nz*ny;	/* CDEF */
		c[ix] = *(a+offset);		/* copy array */
		*(a+offset) = 0.0;		/* reset for accumulation */
	}
		
	for (ix=0; ix<nx; ix++) 
	    for (jx=0; jx<nb; jx++) {
	        kx = ix + jx - (nb-1)/2;
		if (kx>=0)
		    if (kx<nx) {
                        if (Qbad && c[ix]==bad) continue;
		        *(a+iz+iy*nz+kx*ny*nz) += b[jx]*c[ix];
		    } else
			continue;
	    }
	return 1;
}

int convolve_y (a, ix, iz, nx, ny, nz, b, nb)
real *a, b[];
int    nx, ny, nz, nb, ix, iz;
{
	real c[MSIZE];
	int    iy, jy, ky, offset;
	
	if (ny>MSIZE) {
		warning("convolve_y: MSIZE=%d array too small for input %d\n",
			MSIZE,ny);
                return 0;
	}
	
	for (iy=0; iy<ny; iy++) {
	    offset = iz+iy*nz+ix*nz*ny;     /* CDEF */
	    c[iy] = *(a+offset);		/* copy array */
	    *(a+offset) = 0.0;		/* reset for accumulation */
	}
		
	for (iy=0; iy<ny; iy++) 
    	    for (jy=0; jy<nb; jy++) {
	        ky = iy + jy - (nb-1)/2;
		if (ky>=0)
		    if (ky<ny){
                        if (Qbad && c[iy]==bad) continue;
		        *(a+iz+ky*nz+ix*ny*nz) += b[jy]*c[iy];
		    } else
			continue;
	    }
	return 1;
}

int convolve_z (a, ix, iy, nx, ny, nz, b, nb)
real *a, b[];
int    nx, ny, nz, nb, ix, iy;
{
	real c[MSIZE];
	int    iz, jz, kz, offset;
	
	if (nz>MSIZE) {
		warning("convolve_z: MSIZE=%d array too small for input %d\n",
			MSIZE,nz);
		return 0;
	}
	
	for (iz=0; iz<nz; iz++) {
            offset = iz+iy*nz+ix*nz*ny;         /* CDEF */
            c[iz] = *(a+offset);		/* copy array */
            *(a+offset) = 0.0;		/* reset for accumulation */
	}
		
	for (iz=0; iz<nz; iz++) 
            for (jz=0; jz<nb; jz++) {
		kz = iz + jz - (nb-1)/2;
		if (kz>=0)
	            if (kz<nz) {
                        if (Qbad && c[iz]==bad) continue;
			*(a+kz+iy*nz+ix*nz*ny) += b[jz]*c[iz];
		    } else
			continue;
	    }
	return 1;
}



void wiener(void)
{
#if 0
  imageptr itmp = NULL;
  int ix, iy, ixd, iyd; ix1, iy1;
  
  itmp = create_image(iptr,nx,ny);

  for (iz=0; iz< nz; iz++) {
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++)   	
	MapValue(iptr,ix,iy,iz) = CubeValue(iptr,ix,iy,iz);
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++) {
	sumi = sumii = 0.0;
	for (iyd=-nxw; iyd<=nxw; iyd++) {
	  iy1 = iy + iyd;
	  if (iy1<0 || iy1>=ny) continue;
	  for (ixd=-nxw; ixd<=nxw; ixd++) {
	    ix1 = ix + ixd;
	    if (ix1<0 || ix1>=nz) continue;
	    sumi  += MapValue(iptr,ix1,iy1,iz);
	    sumii += MapValue(iptr,ix1,iy1,iz) * MapValue(iptr,ix1,iy1,iz);
	    
	    
	    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++)   	
          brightness = CubeValue(iptr,ix,iy,iz);
	  total += brightness;
          m_max = MAX(m_max, brightness);
          m_min = MIN(m_min, brightness);
	  
#else
	  error("wiener not implmented");
#endif
}
