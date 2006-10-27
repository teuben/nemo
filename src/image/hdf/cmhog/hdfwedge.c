/*
 *  HDFWEDGE: cloned off hdfgrid dec 17, 2004 to accomodate wedges that
 *            can optionally be made into a full plane using a periodic=t
 *            keyword.  Shetty, Ostriker, Teuben
 *
 *     17-dec-2004    V0.1 Initial version, only deals with select=
 *                         no copies of data made, so zvar= does not work
 *                         periodic=f does not work either
 *     20-dec-2004    V0.2 patched up too wide edges for non-integral 2pi intervals 
 *     23-dec-2004    V0.3 fixed periodic for np=odd
 *     24-dec-2004    V0.4 fixed for periodic=f
 *     26-oct-2006    V0.5 attempt to handle 3D cubes by selecting a 'Z' value
 *
 * TODO
 *  - check the interpolation on the 2nd and 3rd quadrant, this this is where
 *    hdfgrid and hdfwedge differ on the rounding level
 *    but with the jp=idx[ip+1] it's bad (rounding or higher) all over the plane
 *  - enable zvar=
 *  - enable periodic=f, it's currently broken
 *
 */

 
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <history.h>
#include <image.h>
#include <grid.h>

#ifdef INC_HDF
#include <hdf.h> 	/* some conflicts with nemo include files */
#endif

string defv[] = {
    "in=???\n			Input file (HDF SD)",
    "out=???\n                  Output image file",
    "select=1\n			Select which SDS from the file (1=first)",
    "nx=256\n			Number of pixels in X",
    "ny=256\n			Number of pixels in Y",
    "xrange=-16:16\n		Range in X",
    "yrange=-16:16\n		Range in Y",
    "zvar=\n                    Optional selections: {vr,vt,den,vx,vy}",
    "periodic=t\n               Attempt to fill the plan by cloning the wedge N times",
    "z=0\n                      Select a z-slice (0=first) from the RPZ cube if NZ > 1",
    "order=zpr\n                Order of HDF axes",
    "VERSION=0.5\n		26-oct-06 PJT",
    NULL,
};

string usage="Regrid a CMHOG polar (RPZ) HDF SDS wedge to a cartesian NEMO image";

string cvsid="$Id$";

/* 3rd dimension is ignored though.... */
#define MAXRANK 3

#define OOST    0.7071067812        /*  (float) 1/sqrt(2)  */

local int rank, shape[MAXRANK], run[MAXRANK];
local char label[256], unit[256], fmt[256], coordsys[256];
local int it0;

local void fshift(int n, float *a, int i0);
local void setrange(real *rval, string rexp, int nr);
local real string_real(string,string);
float *expand_array(float *a, int size, int n, int *idx);


void nemo_main()
{
  float **image, **coord, *buffer, *rads, *phis, *phisf, rad, phi, phi_orig, 
    dmin, dmax, phi1, phi2, rad1, rad2;
  float *buffer1, *buffer2, *buffer3, *zero;
  float **image1, **image2;
  int i0, i1, i, j, k, ret, size, type, nsds, isel, nx, ny;
  int ir, ip, jr, jp, nr, np, n1, n2, is, ix, iy, ncopy,nu,nl,n,npf,nzero;
  int zoff, z = getiparam("z");     /* the slice to be selected */
  int *idx;
  char **output, *cbuff;
  char ntype[32];
  string filter, zvar, infile = getparam("in");
  stream outstr;
  real xrange[3], yrange[3], cosp, sinp, pratio, deltap,edge;
  real x, y, a1, a2, c1,c2,c3,c4, cmin, cmax, dcon, tmp, vr, vt;
  real sds_time = -1.0;
  imageptr iptr, iptr0;
  bool mirror, first = TRUE, both = FALSE, flip=FALSE, Qrot;
  bool Qmirror, Qperiodic  = getbparam("periodic");
  Grid gx, gy;
  string axis = getparam("order");

  if (Qperiodic)
    warning("Periodic=f is currently broken");

  if (!streq(axis,"zpr"))
    error("sorry...cannot handle cubes that are not in ZPR order....");
  
  nsds = DFSDndatasets(infile);
  if (nsds<0) 
    error("%s is probably not an HDF scientific dataset",infile);
  dprintf(1,"Found %d scientific data set%s in %s\n",
	  nsds,(nsds > 1 ? "s" : ""),infile);
  if (hasvalue("zvar")) {
    zvar = getparam("zvar");
    dprintf(0,"Polar CMHOG Gridding variable: %s\n",zvar);
    if (streq(zvar,"vr"))
      isel = 1;
    else if (streq(zvar,"vt"))
      isel = 2;
    else if (streq(zvar,"den"))
      isel = 3;
    else if (streq(zvar,"vx") || streq(zvar,"vy") || streq(zvar,"vp") || streq(zvar,"vm")) {
      both = TRUE;
      flip = TRUE;
      isel = 1;   /* !! will need to read 1 + 2 !! */
      Qrot = (streq(zvar,"vp") || streq(zvar,"vm"));
    } else
      error("Unsupported gridding output variable %s",zvar);
  } else
    isel = getiparam("select");
  if (isel>nsds || isel<1) error("Illegal isel, must be 1..%d",nsds);
  
  for (k=0; k<nsds; k++) {        /* read until we have the right SDS */
    ret = DFSDgetdims(infile,&rank, shape, MAXRANK);
    if (ret < 0) error("Problem getting rank/shape at SDS #%d, MAXRANK=%d",k+1,MAXRANK);
    
    label[0] = unit[0] = fmt[0] = coordsys[0] = 0;
    ret = DFSDgetdatastrs(label, unit, fmt, coordsys);
    ret = DFSDgetNT(&type);
    
    if (sds_time < 0)       /* search for: "AT TIME=" : not all have them!! */
      sds_time = string_real(label,"AT TIME=");
    
    if (k==0) {				/* first time: allocate */
      if (rank == 2) {
	i0=0;      /* start at the first dimension */
      } else if (rank == 3) {
	dprintf(1,"Shape: %c=%d %c=%d %c=%d\n",axis[0],shape[0],axis[1],shape[1],axis[2],shape[2]);
	i0=1;      /* start at the first dimension */
      } else
	error("Cannot handle rank=%d, must be 2 or 3",rank);
      i1=i0+1;     /* radius dim is always one higher than angle (in HDF) */
      coord = (float **) allocate(rank*sizeof(float *));
      for (i=0, size=1; i<rank; i++) {     /*  process all axes */
	size *= shape[i];
	coord[i] = (float *) allocate(shape[i] * sizeof(float));
	ret = DFSDgetdimscale(i+1, shape[i],coord[i]);
	if (ret<0) error("getting shape[%d]",i+1);
	cmin = cmax = coord[i][0];
	for (j=0; j<shape[i]; j++) {
	  cmin = MIN(cmin, coord[i][j]);
	  cmax = MAX(cmax, coord[i][j]);
	}
	dprintf(0,"HDF Dimension %d Order %c Size %d CoordinateValue Min %g Max %g %s\n",
		i+1,axis[i],shape[i],cmin,cmax,
		shape[i] == 1 ? "(** dimension will be skipped **)" : "");
      }
      dprintf(1,"i0=%d i1=%d\n",i0,i1);
      if (shape[i0] < 2)
	error("Sorry, don't know how to deal with 1D maps (phi axis has %d pixel)",shape[i0]);
      nzero = shape[i1];    /* an extra array for just 0's in case periodic=f ?? [i0] ?? */
      dprintf(1,"DFSDgetdata() will use a buffersize %d+%d\n",size,nzero);
      buffer1 = (float *) allocate((size+nzero) * sizeof(float));
      if (both) buffer2 = (float *) allocate((size+nzero) * sizeof(float));
    }
    if (k+1 == isel) {          /* if we got the right one, get data now */
      if (both && isel==2)
	buffer = buffer2;
      else
	buffer = buffer1;
      ret = DFSDgetdata(infile,rank,shape,buffer);
      dmin = dmax = buffer[0];
      for (i=1; i<size; i++) {
	dmin = MIN(dmin,buffer[i]);
	dmax = MAX(dmax,buffer[i]);
      }
      dprintf(0,"%d Datamin/max read: %g %g\n",isel,dmin,dmax);
      if (both) {
	isel++;
	if (isel>2) break;
      } else
	break;
    }
  } /* k */

  if (rank==3 && axis[0]=='z')
    zoff = z * shape[1]*shape[2];
  else 
    zoff = 0;
      

  image1 = (float **) allocate((shape[i0]+1) *sizeof(float *));
  for (i=0; i<shape[i0]+1; i++)        /* notice extra row of zero also pointed to */
    image1[i] = &buffer1[zoff+i*shape[i1]];

  if (both) {
    error("Cannot combine SDS yet");
    image2 = (float **) allocate(shape[i0] *sizeof(float *));
    for (i=0; i<shape[i0]; i++)
      image2[i] = &buffer2[i*shape[i1]];
    if (streq(zvar,"vx")) {
      for (i=0; i<shape[i0]; i++) {    /* phi */
	cosp = cos((double)coord[i0][i]);
	sinp = sin((double)coord[i0][i]);
	for (j=0; j<shape[i1]; j++) {     /* rad */
	  vr = image1[i][j];
	  vt = image2[i][j];
	  image1[i][j] = vr*cosp-vt*sinp;
	}
      }
    } else if (streq(zvar,"vy")) {
      for (i=0; i<shape[i0]; i++) {    /* phi */
	cosp = cos((double)coord[i0][i]);
	sinp = sin((double)coord[i0][i]);
	for (j=0; j<shape[i1]; j++) {     /* rad */
	  vr = image1[i][j];
	  vt = image2[i][j];
	  image1[i][j] = vr*sinp+vt*cosp;
	}
      }
    } else if (streq(zvar,"vm")) {
      for (i=0; i<shape[i0]; i++) {    /* phi */
	cosp = cos((double)coord[i0][i]);
	sinp = sin((double)coord[i0][i]);
	for (j=0; j<shape[i1]; j++) {     /* rad */
	  vr = image1[i][j];
	  vt = image2[i][j];
	  image1[i][j] = ((vr-vt)*cosp - (vr+vt)*sinp)*OOST;
	}
      }
    } else if (streq(zvar,"vp")) {
      for (i=0; i<shape[i0]; i++) {    /* phi */
	cosp = cos((double)coord[i0][i]);
	sinp = sin((double)coord[i0][i]);
	for (j=0; j<shape[i1]; j++) {     /* rad */
	  vr = image1[i][j];
	  vt = image2[i][j];
	  image1[i][j] = ((vr+vt)*cosp + (vr-vt)*sinp)*OOST;
	}
      }
    }
  } /* both */
  
  image = image1;
  /* the image can now be referred to as in:  image[i_phi][i_rad]  */

  nx = getiparam("nx");
  ny = getiparam("ny");
  setrange(xrange,getparam("xrange"),nx);
  setrange(yrange,getparam("yrange"),ny);
  nr = shape[i1];
  np = shape[i0];
  rads = coord[i1];
  phis = coord[i0];

  zero = image[np];
  for (i=0; i<nr; i++) zero[i] = 0.0;

  /* 2st method, see how many sections of the phi range fit in 2.pi */
  pratio = TWO_PI/(phis[np-1]-phis[0]);
  ncopy = (int) lrint(pratio);    /* nearest integer ? */
  dprintf(0,"Periodicity (%g) appears to be %d\n",pratio,ncopy);

  /* 2nd method, count steps to patch lower and upper coordinates into -pi .. pi */
  /* this is subject to rounding error, since we're using float's and add them   */
  deltap = phis[1]-phis[0];
  dprintf(1,"dPHI=%g Range=%g\n",deltap,phis[np-1]-phis[0]);
  if (deltap < 0) error("Cannot handle inverted PHI axis: dPHI=%g",deltap);
  for (nu=0, edge=phis[np-1]; ;nu++) {
    edge += deltap;
    dprintf(2,"Upper test: %d %g\n",nu,edge);
    if (edge > PI) break;
  }
  for (nl=0, edge=phis[0]; ;nl++) {
    edge -= deltap;
    dprintf(2,"Lower test: %d %g\n",nl,edge);
    if (edge < -PI) break;
  }
  npf = nl + np + nu;
  if (npf != ncopy*np)
    warning("Lower (%d) and Upper (%d) patch to fill %d into -pi..pi not exactly even (%d)",
	    nl,nu,np,ncopy*np);
  dprintf(1,"Nupper=%d Nlower=%d totals %d, compare to %d\n",
	  nu,nl,npf,ncopy*np);

  /* now create an index array such that image_filled[idx[i]] points to image[i] */
  /* and patch the upper and lower pieces to point to each other                 */
  /* also make a new 'phis' coordinate array for the full  plane, called 'phisf' */
  /* the extra +2 is a shadow location:                                          */
  /* 0 is to point just below -pi, npf+1 just above pi                           */
  idx = (int *) allocate((npf+2)*sizeof(int));
  phisf = (float *) allocate((npf+2)*sizeof(float));
  if (nl+nu > 0) {
    for (i=0; i<np; i++) {                              /* place the original array with the known offset */
      idx[i+nl+1] = i;
      phisf[i+nl+1] = phis[i];
    }
    for (i=0, j=0; i<nu; i++) {                         /* iteratively add in the upper parts */
      idx[i+nl+np+1] = Qperiodic ? j : np;
      phisf[i+nl+np+1] = phisf[i+nl+np] + deltap;
      j++;
      if (j==np) j=0;    
    }
    for (i=0, j=np-1; i<nl; i++) {                      /* iteratively add in the lower parts */
      idx[nl-i] = Qperiodic ? j : np;
      phisf[nl-i] = phisf[nl-i+1] - deltap;
      if (j==0) j=np;
      j--;
    }
    idx[0]     = Qperiodic ? idx[npf] : np;
    idx[npf+1] = Qperiodic ? idx[1]   : np;
    phisf[0]     = phisf[npf] - TWO_PI;
    phisf[npf+1] = phisf[1]   + TWO_PI;
    
    n = idx[1]-idx[0];
    if (Qperiodic && n != 1) {
      warning("First shadow cell covers %d, trying to fix it:",n);
      dprintf(0,"idx[1] = %d idx[0] was : %d  is: %d\n", idx[1],idx[0],idx[1]-1);
      idx[0] = idx[1]-1;
      phisf[0] += deltap;
    }
    n = idx[npf+1]-idx[npf];
    if (Qperiodic && n != 1) {
      warning("Last shadow cell covers %d, trying to fix it:",n);
      dprintf(0,"idx[N] = %d idx[N+1] was : %d  is: %d\n", idx[npf],idx[npf+1],idx[npf]+1);
      idx[npf+1] = idx[npf]+1;
      phisf[npf+1] -= deltap;
    }
  } else {
    warning("should not have to do this section ?? ");
    idx[0] = np-1;
    for (i=0; i<np; i++) idx[i+1] = i;
    idx[np+1] = 0;
  }
  
  for (i=0; i<npf+2; i++)
    dprintf(2,"idx[%d] = %d/%d    %g   %d %d\n",i,idx[i],np,phisf[i],nl,nu);

  outstr = stropen(getparam("out"),"w");
  create_image(&iptr, nx, ny);
  
  n1 = n2 = 0;
  dmin = dmax = 0.0;    /* just in case all points outside grid */
  for (j=0; j<ny; j++) {                  /* loop over all rows */
    y = yrange[0] + (j+0.5)*yrange[2];
    for (i=0; i<nx; i++) {              /* loop over all columns */
      x = xrange[0] + (i+0.5)*xrange[2];
      rad = sqrt(x*x + y*y);
      phi = atan2(y,x);               /* phi will be in the -pi : pi range */
      
      if (rad > rads[nr-1] || rad < rads[0]) {     /* outside radial grid... */	
	dprintf(3,"rad %g not in  %g : %g range\n",rad,rads[0],rads[nr-1]);
	MapValue(iptr,i,j) = 0.0;
	n1++;
	continue;
      }
      if (phi > phisf[npf+1] || phi < phisf[0]) {  /* should never happen */
	dprintf(0,"%g not in %g : %g range....???\n",phi,phisf[0],phisf[npf-1]);
	MapValue(iptr,i,j) = 0.0;
	n2++;
	continue;
      }
      for (ir=0; ir<nr-1; ir++)           /* simple search in Rad */
	if (rads[ir+1] > rad) break;
      jr = ir+1;
      rad1 = rads[ir];
      rad2 = rads[jr];

      for (ip=0; ip<npf+1; ip++)
	if (phisf[ip+1] > phi) break;
      /* get full angles */
      phi1 = phisf[ip];
      phi2 = phisf[ip+1];
      /* but index in original range to get at the data */
      ip = idx[ip];
#if 1
      if (Qperiodic) {
	jp = ip+1;
	if (jp>=np) 
	  jp=0;    /* this is a nasty hack but doesn't work for data where delta_shadow > 1 */
      } else
	jp = idx[ip+1];      /* something wrong with this :-) ??? */
#endif
      
      dprintf(2,"(x,y)%d %d (r,p) %g %g [@ %d:%d %d:%d] LL: %g %g\n",
	      i,j,rad,phi,ir,jr,ip,jp,rads[ir],phi1);
      
      c1 = (rad-rad2)/(rad1-rad2);
      c2 = (rad-rad1)/(rad2-rad1);
      c3 = (phi-phi2)/(phi1-phi2); 
      c4 = (phi-phi1)/(phi2-phi1); 
      if (c1<0 || c2<0 || c3<0 || c4<0 || c1>1 || c2>1 || c3>1 || c4>1) {  /* should not happen */
	dprintf(0,"c1234= %g %g %g %g\n",c1,c2,c3,c4);
	dprintf(0,"  phi: %g < %g < %g ? %d \n",phi1,phi,phi2,ip);
      } 
      a1 = image[ip][ir]*c1 + image[ip][jr]*c2;
      a2 = image[jp][ir]*c1 + image[jp][jr]*c2;
      if (both) { 
	error("cannot handle two images yet");
	if (flip) 
	  MapValue(iptr,i,j) = -a1*c3 + a2*c4;
	else if (mirror)
	  MapValue(iptr,i,j) = -a1*c3 - a2*c4;
      } else
	MapValue(iptr,i,j) = a1*c3 + a2*c4; 

      if (first) {
	dmin = dmax = MapValue(iptr,i,j);
	first = FALSE;
      } else {
	dmin = MIN(dmin, MapValue(iptr,i,j));
	dmax = MAX(dmax, MapValue(iptr,i,j));
      }
    } /* i */
  } /* j */
  dprintf(0,"%d/%d out of %d cartesian points outside grid in rad/phi\n",n1,n2,nx*ny);
  
  /* finish off the header */
  Xmin(iptr) = xrange[0];
  Ymin(iptr) = yrange[0];
  Dx(iptr) = xrange[2];
  Dy(iptr) = yrange[2];
  MapMin(iptr) = dmin;
  MapMax(iptr) = dmax;
  Unit(iptr) = unit;
  Time(iptr) = sds_time;
  
  set_headline(label);
  write_image(outstr,iptr);
  strclose(outstr);
  dprintf(0,"Datamin/max written: %g %g\n",dmin,dmax);
}

local void setrange(real *rval, string rexp, int nr)
{
  char *cptr, first[80];
  int n;
  
  cptr = strchr(rexp, ':');
  if (cptr) {
    n = cptr-rexp;
    if(n>79)error("Not enough string space in setrange[80]...");
    strncpy(first,rexp,n);
    first[n] = 0;
    if (nemoinpr(first,&rval[0],1) < 1) error("Parsing first part in %s",rexp);
    if (nemoinpr(cptr+1,&rval[1],1) < 1) error("Parsing second part in %s",rexp);
  } else {
    rval[0] = 0.0;
    if (nemoinpr(rexp,&rval[1],1) < 1) error("Parsing single %s",rexp);
  }
  rval[2] = (rval[1] - rval[0]) / nr;
}

local real string_real(string label, string key)
{
  char *cp = label;
  int klen = strlen(key);
  
  while (*cp) {
    if (strncmp(cp,key,klen)==0) {
      cp += klen;
      if (*cp == 0) 
	error("string_real conversion: %s %s",label,key);
      return atof(cp);
    }
    cp++;
  }
  warning("string_real conversion: %s %s",label,key);
  return -1.0;
}


/*
 * shift an array of size n by k elements (the slow way)
 */

local void fshift(int n, float *a, int k)
{
  int i, j;
  float tmp;
  
  for(j=0; j<k; j++) {
    tmp = a[n-1];
    for (i=n-1; i>0; i--) {
      a[i] = a[i-1];
    }
    a[0] = tmp;
  }
}

/* 
 * expands the array 2nd dimension from the back
 * such that image[nz][ny][nx] -> image[nz][ncopy*ny][nx]
 * 
 */

float *expand_array(float *a, int size, int n, int *idx)
{
  int x0, y0, z0, nx0, ny0, nz0;
  float *b = (float *) allocate(size*sizeof(float));

  if (n==2) {
    nx0 = idx[0];
    ny0 = idx[1];
    dprintf(0,"reorder: %d x %d\n",nx0,ny0);
    for (x0=0; x0<nx0; x0++)
      for (y0=0; y0<ny0; y0++)
	b[x0+nx0*y0] = a[y0+ny0*x0];
  } else {
    nx0 = idx[0];
    ny0 = idx[1];
    nz0 = idx[2];
    dprintf(0,"reorder: %d x %d x %d\n",nx0,ny0,nz0);
    for (x0=0; x0<nx0; x0++)
      for (y0=0; y0<ny0; y0++)
	for (z0=0; z0<nz0; z0++)
	b[x0+nx0*y0+nx0*ny0*z0] = a[z0+nz0*y0+nz0*ny0*x0];
  }		
  free(a);
  return b;
}

