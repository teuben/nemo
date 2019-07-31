/* 
 * CCDMOM: take a moment along an axis in an image/cube
 *
 *	quick and dirty:  8-jun-95		pjt
 *      19-oct-95   very dirty axis=3           pjt
 *	12-dec-98   0.3  more implementations   pjt
 *      31-dec-98   0.4  keep reduced axis, and replace it with new value PJT
 *      20-feb-01   0.4a implemented special mom=3 for a 3point poly fit to peak  PJT
 *                       also fixed ccdmom 
 *      25-mar-01   0.5  added mom=-1 as mean value along an axis, fixed bugs
 *      18-oct-05   0.6  added mom=-2 as dispersion around mean along an axis
 * ??   16-sep-11   ???  added clip= and  rngmsk= ???   code lost ???
 *      19-jul-12   0.7  allow peak (mom=3) to find 2nd peak
 *      27-nov-12   1.0  add oper=  to insert an operator (ie. out = in <oper> out )
 *       8-dec-12   1.1  allow mom=-3 for differentials (axis=3 only for now) in 2..Nz()
 *      13-feb-13   2.0  default integration, instead of just summing
 *      29-apr-13   2.1  add clumping definition  http://arxiv.org/abs/1304.1586  (mom=-4)
 *      12-jan-16   2.2a minmax computation was forgotten
 *      24-jun-16   2.3  mom=30,31,32,33
 *      30-apr-17   2.3c various for EDGE clarifications
 *      10-may-17   2.4  add cont subtraction for modes 30,31,32,33,...
 *      15-jun-17   2.5  pos=x,y triggering debug output
 *      21-jun-17   2.6  use abs values for
 *      25-sep-18   2.7  tinkering because of "bettermoments"
        29-jul-19   2.7c fix bug when no clip was given
 *                      
 * TODO : cumulative along an axis, sort of like numarray.accumulate()
 *        man page talks about clip= and  rngmsk=, where is this code?
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h>

/* this may slow the code down if you don't ever need the pos= keyword */
#define USE_POS

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "axis=3\n       Axis to take moment along (1=x 2=y 3=z)",
  "mom=0\n	  Moment to take [0=sum,1=mean loc,2=disp loc,3=peak loc,4=peak mom1,-1=mean val,-2=disp val,-3=clump]",
  "keep=f\n	  Keep moment axis in full length, and replace all values",
  "cumulative=f\n Cumulative axis (only valid for mom=0)",
  "oper=\n        Operator on output (enforces keep=t)",
  "peak=0\n       Find N-th peak in case of peak finding (mom=3) [0 means first as well]",
  "clip=\n        If used, clip values between -clip,clip or clip1,clip2",
  "rngmsk=f\n     Invalidate pixel when first moment falls outside range of valid axis [not impl]",
  "abs=f\n        Use absolute values in moment calculations (abs=t is better for low S/N)",
  "zero=f\n       Mask out values below zero for mom=30,31,32 (**experimental**)",
  "contsub=t\n    Fake a continuum subtraction for mom=30,31,32  (**experimental**)",
  "integrate=t\n  Use integration instead of just summing, only used for mom=0",
#ifdef USE_POS  
  "pos=\n         Given this (i,j) [0 based] location, debug output spectral information",
#else  
  "pos=\n         ** keyword disabled via the #ifdef USE_POS **",
#endif  
  "VERSION=2.7c\n 30-jul-2019 PJT",
  NULL,
};

string usage = "moment along an axis of an image";
string cvsid="$Id$";

local real peak_spectrum(int n, real *spec, int p);
local real peak_mom(int n, real *spec, int *smask, int peak, int mom, bool Qcontsub, bool Qabs, bool Qzero);
local real peak_axis(imageptr iptr, int i, int j, int k, int axis);
local int  peak_find(int n, real *data, int *mask, int npeak);
local int  peak_assign(int n, real *data, int *mask);
local bool out_of_range(real *clip, real x);
local void image_oper(imageptr ip1, string oper, imageptr ip2);



void nemo_main()
{
    stream  instr, outstr;
    string  oper;
    int     i,j,k,nx, ny, nz, nx1, ny1, nz1;
    int     pos[2];
    int     ii, axis, mom;
    int     nclip, apeak, apeak1, cnt;
    imageptr iptr=NULL, iptr1=NULL, iptr2=NULL;      /* pointer to images */
    real    tmp0, tmp1, tmp2, tmp00, newvalue, peakvalue, scale, offset;
    real    *spec, ifactor, cv, clip[2], m_min, m_max;
    int     *smask;
    bool    Qkeep = getbparam("keep");
    bool    Qoper = hasvalue("oper");
    int     npeak = getiparam("peak");
    bool    Qint  = getbparam("integrate"); 
    bool    Qclip = hasvalue("clip");
    bool    Qpos  = hasvalue("pos");
    bool    Qabs  = getbparam("abs");
    bool    Qzero = getbparam("zero");
    bool    Qcontsub = getbparam("contsub");

    if (Qoper) {
      Qkeep = TRUE;
      oper = getparam("oper");
    }

    instr = stropen(getparam("in"), "r");
    mom = getiparam("mom");
    if (mom < -4)  error("Illegal value mom=%d",mom);
    axis = getiparam("axis");
    if (axis < 0 || axis > 3) error("Illegal value axis=%d",axis);

    if ((mom%10==3 ) && axis!=3 && npeak>1) error("Nth-peak>1 finding only axis=3");

    if (Qpos)
      if (nemoinpi(getparam("pos"),pos,2) != 2)
	error("Need 2 integer (0 based) values for x,y position in pos=%s",getparam("pos"));

    if (Qclip) {
      nclip = nemoinpr(getparam("clip"),clip,2);
      if (nclip<1) error("error parsing clip=%s",getparam("clip"));
      if (nclip==1) {
	clip[1] =  clip[0];
	clip[0] = -clip[1];
      }
    }

    if (getbparam("cumulative"))
      axis = -axis;

    read_image( instr, &iptr);
    nx1 = nx = Nx(iptr);	
    ny1 = ny = Ny(iptr);
    nz1 = nz = Nz(iptr);

    if (mom==-4) {   /* clumping hack : only works for axis=3 */
      if (axis==3) {
	for (k=0; k<nz; k++) {
	  tmp0 = tmp1 = tmp2 = 0.0;
	  for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {	
	      cv = CubeValue(iptr,i,j,k);
	      if (Qclip && (clip[0]<=cv && cv<=clip[1])) {
		dprintf(1,"%d %d %d  %f %f %f\n",i,j,k,cv,clip[0],clip[1]);
		continue;
	      }
	      tmp0 += 1.0;
	      tmp1 += cv;
	      tmp2 += cv*cv;
	    }
	  }
	  if (tmp1 != 0.0)
	    printf("%d %f %.0f\n",k,(tmp2*tmp0)/(tmp1*tmp1),tmp0);
	  else
	    printf("%d 0.0 0\n",k);
	}
      } else 
	error("axis=%d not yet supported for mom=%d",axis,mom);
      free_image(iptr);
      return;
    }

    if (Qkeep) {
        dprintf(0,"Keeping %d*%d*%d cube\n",nx1,ny1,nz1);
        if (axis==1) {
	    spec = (real *) allocate(nx*sizeof(real));
	    smask = (int *) allocate(nx*sizeof(int));
        } else if (axis==2) {
	    spec = (real *) allocate(ny*sizeof(real));
	    smask = (int *) allocate(ny*sizeof(int));
        } else if (axis==3) {
	    spec = (real *) allocate(nz*sizeof(real));
	    smask = (int *) allocate(nz*sizeof(int));
        } else 
            error("Invalid axis: %d (Valid: 1,2,3)",axis);
    } else {
        if (axis==1) {
            nx1 = 1;    ny1 = ny;   nz1 = nz;
	    spec = (real *) allocate(nx*sizeof(real));
	    smask = (int *) allocate(nx*sizeof(int));
        } else if (axis==2) {
            nx1 = nx;   ny1 = 1;    nz1 = nz;
	    spec = (real *) allocate(ny*sizeof(real));
	    smask = (int *) allocate(ny*sizeof(int));
        } else if (axis==3) {
            nx1 = nx;   ny1 = ny;   nz1 = 1;
	    spec = (real *) allocate(nz*sizeof(real));
	    smask = (int *) allocate(nz*sizeof(int));
        } else if (axis < 0) {
	    nx1 = nx;   ny1 = ny;   nz1 = nz;
	    spec = NULL;
	    smask = NULL;
	} else
            error("Invalid axis: %d (Valid: 1,2,3)",axis);
        dprintf(1,"Reducing %d*%d*%d to a %d*%d*%d cube\n",
                   nx,ny,nz, nx1,ny1,nz1);
    }

    outstr = stropen(getparam("out"), "w");

    if (axis > 0) {
      create_cube(&iptr1,nx1,ny1,nz1);
      create_cube(&iptr2,nx1,ny1,nz1);
    } else {
      copy_image(iptr,&iptr1);
    }

    ifactor = 1.0;
    if (axis==1) {
      scale = Dx(iptr);
      offset = Xmin(iptr);
      if (Qint) ifactor *= ABS(Dx(iptr));
      for (k=0; k<nz; k++)
        for (j=0; j<ny; j++) {
	    tmp0 = tmp00 = tmp1 = tmp2 = 0.0;
	    cnt = 0;
	    peakvalue = CubeValue(iptr,0,j,k);
            for (i=0; i<nx; i++) {
	        spec[i] = CubeValue(iptr,i,j,k);
 	        if (Qclip && out_of_range(clip,CubeValue(iptr,i,j,k))) continue;
		cnt++;
    	        tmp0 += CubeValue(iptr,i,j,k);
		tmp00 += sqr(CubeValue(iptr,i,j,k));
    	        tmp1 += i*CubeValue(iptr,i,j,k);
    	        tmp2 += i*i*CubeValue(iptr,i,j,k);
		if (CubeValue(iptr,i,j,k) > peakvalue) {
		  apeak = i;
		  peakvalue = CubeValue(iptr,i,j,k);
		}
            }
	    if (cnt==0 || (tmp0==0.0 && tmp00==0.0)) {
	      newvalue = 0.0;
	    } else {
	      if (mom==-1)
		newvalue = tmp0/cnt;
	      else if (mom==-2)
		newvalue = sqrt(tmp00/cnt - sqr(tmp0/cnt));
	      else if (mom==0) 
		newvalue = tmp0 * ifactor;
	      else if (mom==1)
		newvalue = offset + scale*(tmp1/tmp0);
	      else if (mom==2)
		newvalue = scale*sqrt(tmp2/tmp0 - sqr(tmp1/tmp0));
	      else if (mom==3)
		newvalue = scale*(apeak + peak_axis(iptr,apeak,j,k,axis)) + offset;
	    }
	    for (i=0; i<nx1; i++)
	      CubeValue(iptr1,i,j,k) = newvalue;
        }

        /* TODO: fix up Qkeep headers */

        Xmin(iptr1) = Xmin(iptr) + 0.5*(nx-1)*Dx(iptr);
        Ymin(iptr1) = Ymin(iptr);
        Zmin(iptr1) = Zmin(iptr);
        Dx(iptr1) = nx * Dx(iptr);
        Dy(iptr1) = Dy(iptr);
        Dz(iptr1) = Dz(iptr);
        
        Namex(iptr1) = Namex(iptr); /* care: we're passing a pointer */
        Namey(iptr1) = Namey(iptr);
        Namez(iptr1) = Namez(iptr);

	if (Qoper) image_oper(iptr,oper,iptr1);

    } else if (axis==2) {                      
      scale = Dy(iptr);
      offset = Ymin(iptr);
      if (Qint) ifactor *= ABS(Dy(iptr));
      for (k=0; k<nz; k++)
        for (i=0; i<nx; i++) {
            tmp0 = tmp00 = tmp1 = tmp2 = 0.0;
	    cnt = 0;
	    peakvalue = CubeValue(iptr,i,0,k);
            for (j=0; j<ny; j++) {
	        spec[j] = CubeValue(iptr,i,j,k); 
 	        if (Qclip && out_of_range(clip,CubeValue(iptr,i,j,k))) continue;
		cnt++;
    	        tmp0 += CubeValue(iptr,i,j,k);
		tmp00 += sqr(CubeValue(iptr,i,j,k));
    	        tmp1 += j*CubeValue(iptr,i,j,k);
    	        tmp2 += j*j*CubeValue(iptr,i,j,k);
		if (CubeValue(iptr,i,j,k) > peakvalue) {
		  apeak = j;
		  peakvalue = CubeValue(iptr,i,j,k);
		}
            }
	    if (cnt==0 || (tmp0==0.0 && tmp00==0.0)) {
	      newvalue = 0.0;
	    } else {
	      if (mom==-1)
		newvalue = tmp0/cnt;
	      else if (mom==-2)
		newvalue = sqrt(tmp00/cnt - sqr(tmp0/cnt));
	      else if (mom==0)
		newvalue = tmp0 * ifactor;
	      else if (mom==1)
		newvalue = scale*(tmp1/tmp0) + offset;
	      else if (mom==2)
		newvalue = scale*sqrt(tmp2/tmp0 - sqr(tmp1/tmp0));
	      else if (mom==3)
		newvalue = scale*(apeak + peak_axis(iptr,i,apeak,k,axis)) + offset;
	    }
            for (j=0; j<ny1; j++)
                CubeValue(iptr1,i,j,k) = newvalue;
        }

        /* TODO: */

        Xmin(iptr1) = Xmin(iptr);
        Ymin(iptr1) = Ymin(iptr) + 0.5*(ny-1)*Dy(iptr);
        Zmin(iptr1) = Zmin(iptr);
        Dx(iptr1) = Dx(iptr);
        Dy(iptr1) = ny * Dy(iptr);
        Dz(iptr1) = Dz(iptr);
        
        Namex(iptr1) = Namex(iptr); /* care: we're passing a pointer */
        Namey(iptr1) = Namey(iptr);
        Namez(iptr1) = Namez(iptr);

	if (Qoper) image_oper(iptr,oper,iptr1);

    } else if (axis==3) {                       /* this one is well tested and has more options */
        scale = Dz(iptr);
	offset = Zmin(iptr);
	if (Qint) ifactor *= ABS(Dz(iptr));
    	for(j=0; j<ny; j++) {
      	  for(i=0; i<nx; i++) {                         /* loop over all X and Y positions */
    	    tmp0 = tmp00 = tmp1 = tmp2 = 0.0;
	    cnt = 0;
    	    for(k=0; k<nz; k++) {
	        spec[k] = CubeValue(iptr,i,j,k);
 	        if (Qclip && out_of_range(clip,spec[k])) continue;
		if (cnt==0) {
		  apeak = k;
		  peakvalue = spec[k];
		}
		cnt++;
    	        tmp0  += spec[k];
    	        tmp1  += k*spec[k];
    	        tmp2  += k*k*spec[k];
		tmp00 += sqr(spec[k]);
		if (spec[k] > peakvalue) {
		  apeak = k;
		  peakvalue = spec[k];
		}
		if (mom==-3) {
		  if (k==0)
		    CubeValue(iptr1,i,j,k) = 0;
		  else
		    CubeValue(iptr1,i,j,k) = spec[k]-spec[k-1];
		}
    	    } /* for(k) */
	    if (cnt==0 || (tmp0==0.0 && tmp00==0.0)) {
	      newvalue = 0.0;
	    } else {
	      if (mom==-1)
		newvalue = tmp0/cnt;
	      else if (mom==-2)
		newvalue = sqrt(tmp00/cnt - sqr(tmp0/cnt));
	      else if (mom==0)
		newvalue = tmp0 * ifactor;
	      else if (mom==1)
		newvalue = scale*(tmp1/tmp0) + offset;
	      else if (mom==2) {
		newvalue = tmp2/tmp0 - sqr(tmp1/tmp0);
	        if (newvalue <= 0.0)
		  newvalue = 0.0;
		else
		  newvalue = scale*sqrt(newvalue);
	      } else if (mom==3 || mom/10==3) {  /* mom=3, 30,31,32,33,34 */
		if (npeak == 0) {
		  if (mom==3) {
		      //newvalue = scale*(apeak + peak_axis(iptr,i,j,apeak,axis)) + offset;
		      newvalue = scale*(apeak + peak_spectrum(nz,spec,apeak)) + offset;
		  } else if (mom>=30) {
		      (void) peak_find(nz, spec, smask, 0);                  /* initialize smask */
		      newvalue = peak_mom(nz, spec, smask, 0, mom-30, Qcontsub, Qabs, Qzero);
		      if (mom==31) newvalue = scale*newvalue + offset;
		      if (mom==32) newvalue = scale*newvalue;
		  } else {
		      newvalue = 0.0;
		  }
		} else {
		  (void) peak_find(nz, spec, smask, 0);                  /* initialize smask */
		  apeak1 = peak_find(nz, spec, smask, 1);                /* first peak again */
		  if (apeak1 > 0) {
		    if  (apeak1 != apeak && (apeak!=0 && apeak!=nz-1)) { /* odd if it's not finding the same */
		      for (k=0; k<nz; k++)  
			printf("%d %g  %d\n",k,spec[k],smask[k]);
		      warning("peak_find no good (%d,%d) %d %d",i,j,apeak,apeak1);
		    }
		    newvalue = scale*(apeak1 + peak_spectrum(nz,spec,apeak1)) + offset;
		  } else
		    newvalue = 0.0;

		  if (npeak > 1) { 
		    for (k=2; k<=npeak; k++)
		      apeak1 = peak_find(nz,spec,smask,k);
		    if (mom==3) {
		      if (apeak1 > 0) {
			newvalue = scale*(apeak1 + peak_spectrum(nz,spec,apeak1)) + offset;
		      } else
			newvalue = 0;
		    }
		  }
		  peak_assign(nz, spec, smask);
#if 0
		  for (k=0; k<nz; k++)  
		      printf("%d %g  %g\n",k,spec[k],smask[k]==0 ? 0 : peakvalue/smask[k]);
#endif
		  dprintf(1,"MOM: mom/30\n",mom/30);
		  if (mom >= 30) {
		      newvalue = peak_mom(nz, spec, smask, npeak, mom-30,Qcontsub,Qabs,Qzero);
		      if (mom==31) newvalue = scale*newvalue + offset;
		      if (mom==32) newvalue = scale*newvalue;
#ifdef USE_POS
		      /* DEBUG */
		      if (Qpos && i==pos[0] && j==pos[1]) {
			printf("# spectrum at %d,%d (0 based pixels)\n",i,j);
			for (ii=0; ii<nz; ii++) {
			  printf("%d %g %g %d\n",ii,spec[ii],smask[ii]*spec[ii],smask[ii]);
			}
			printf("# mom = %g\n",newvalue);
		      }
#endif		      
		  }
		} /* npeak */
	      } /* mom */
	    } /* cnt */
	    if (mom>-3)
	      for (k=0; k<nz1; k++)
		CubeValue(iptr1,i,j,k) = newvalue;
	  } /* i */
    	} /* j */

        Xmin(iptr1) = Xmin(iptr);
        Ymin(iptr1) = Ymin(iptr);
        Zmin(iptr1) = Zmin(iptr) + 0.5*(nz-1)*Dz(iptr);
        Dx(iptr1) = Dx(iptr);
        Dy(iptr1) = Dy(iptr);
        Dz(iptr1) = nz * Dz(iptr);
        
        Namex(iptr1) = Namex(iptr); /* care: we're passing a pointer */
        Namey(iptr1) = Namey(iptr);
        Namez(iptr1) = Namez(iptr);

	Beamx(iptr1) = Beamx(iptr);
	Beamy(iptr1) = Beamy(iptr);

	if (Qoper) image_oper(iptr,oper,iptr1);
        
    } else if (axis == -1) {
      for (k=0; k<nz; k++)
        for (j=0; j<ny; j++) {
	  tmp0 = 0.0;
	  for (i=0; i<nx; i++) {
	    tmp0 += CubeValue(iptr,i,j,k);
	    CubeValue(iptr1,i,j,k) = tmp0;
	  }
        }
    } else if (axis == -2) {                      
      for (k=0; k<nz; k++)
        for (i=0; i<nx; i++) {
	  tmp0 = 0.0;
	  for (j=0; j<ny; j++) {
	    tmp0 += CubeValue(iptr,i,j,k);
	    CubeValue(iptr1,i,j,k) = tmp0;
	  }
        }
    } else if (axis == -3) {               
    	for(j=0; j<ny; j++)
	  for(i=0; i<nx; i++) {
    	    tmp0 = 0.0;
    	    for(k=0; k<nz; k++) {
	      tmp0 += CubeValue(iptr,i,j,k);
	      CubeValue(iptr1,i,j,k) = tmp0;
	    }
	  }
    } else
        error("Cannot do axis %d",axis);

    m_min = HUGE;
    m_max = -HUGE;
    for (k=0; k<Nz(iptr1); k++)
    for (j=0; j<Ny(iptr1); j++)
    for (i=0; i<Nx(iptr1); i++) {
      cv = CubeValue(iptr1,i,j,k);
      m_max = MAX(m_max, cv);
      m_min = MIN(m_min, cv);

    }
    MapMin(iptr1) = m_min;
    MapMax(iptr1) = m_max;
    write_image(outstr, iptr1);
}


/*
 * peak_spectrum:
 * return location of peak for (-1,y1) (0,y2) (1,y3)
 * as determined from the 2nd order polynomial going through
 * these 3 points.
 * (also finds valleys)
 *
 * Returns a number between -0.5 and 0.5 from the center
 * because it is guarenteed y2 is a local extremum
 * (though this is NOT checked here)
 *
 * It returns an exact fit, because of the 2nd order poly,
 * could try 5 points and do a lsqfit or gaussfit....
 */

local real peak_spectrum(int n, real *spec, int p)
{
  real y1, y2, y3;

  y1 = spec[p-1];
  y2 = spec[p];
  y3 = spec[p+1];
  if (y1+y3 == 2*y2) return 0.0;
  return 0.5*(y1-y3)/(y1+y3-2*y2);
}

/*
 * @todo:    have an option to use abs() for mom2 calculations,
 *           but this assumes Qcontsub has been applied
 */

local real peak_mom(int n, real *spec, int *smask, int peak, int mom, bool Qcontsub, bool Qabs, bool Qzero)
{
  int i;
  Moment m;
  bool Qfirst = TRUE;
  real cont = 0.0;

  if (Qzero) {
    for (i=0; i<n; i++) {
      if (smask[i] && spec[i] < 0.0) smask[i] = 0;
    }
  }

  ini_moment(&m, mom, n);
  if (Qcontsub) {
    for (i=0; i<n; i++) {
      if (smask[i]==peak) {
	if (Qfirst) {
	  cont = spec[i];
	  Qfirst = FALSE;
	}
	if (spec[i] < cont) cont = spec[i];
      }
    }
  }
  dprintf(1,"peak_mom %d %d %g\n",peak,mom,cont);
  for (i=0; i<n; i++)
    if (smask[i]==peak)
      accum_moment(&m, i, spec[i]-cont);
  if (mom==0)
    return sum_moment(&m);
  else if (mom==1)
    return mean_moment(&m);
  else if (mom==2)
    return sigma_moment(&m);
  else if (mom==3)
    //return skewness_moment(&m); 
    return h3_moment(&m); 
  else if (mom==4)
    //return kurtosis_moment(&m); 
    return h4_moment(&m); 
  else
    return 0.0;
}

/* deprecated ; use peak_spectrum */
 
local real peak_axis(imageptr iptr, int i, int j, int k, int axis)
{
  real y1, y2, y3;
  y2 = CubeValue(iptr,i,j,k);
  if (axis==1)
    i -= 1;
  else if (axis==2)
    j -= 1;
  else if (axis==3)
    k -= 1;
  y1 = CubeValue(iptr,i,j,k);
  if (axis==1)
    i += 2;
  else if (axis==2)
    j += 2;
  else if (axis==3)
    k += 2;
  y3 = CubeValue(iptr,i,j,k);

  if (y1+y3 == 2*y2) return 0.0;
  return 0.5*(y1-y3)/(y1+y3-2*y2);
}

/* 
 * this routine can be called multiple times
 * each time it will find a peak, and then walk down the peak
 * mask[] contains 0 if not part of a peak, 1 for the first peak
 * 2 for the 2nd etc.
 *
 * During the first call, the mask[] array is set to all 0's.
 * i.e. unassigned to a peak, you will need npeak=0
 *
 * Returns: peak location (an integer)  and mask[] was modified
 *
 * TODO:    after all peaks have been find, need to resolve the issue
 *          of possibly re-assigning membership to the neighbor peak
 */


local int peak_find(int n, real *data, int *mask, int npeak)
{
  int i, ipeak=0, apeak=-1;
  real peakvalue, oldvalue;

  dprintf(1,"peak_find %d\n",n);
  if (npeak==0) {               /* initialize by resetting the mask */
    for(i=0; i<n; i++)
      mask[i] = 0;              /* 0 means no peak assigned */
    return -1;
  }

  while (1) {                   /* iterate until a good peak found ? */
    ipeak++;
    dprintf(1,"iter ipeak=%d npeak=%d\n",ipeak,npeak);
    if (ipeak > npeak+3) break;

    for(i=0; i<n; i++) {        /* go over all good points and find the peak */
      if (mask[i]==0) {         /* for all good data not yet assigned to a peak */
	if (apeak<0) {          /* make sure we initialize first good value as peak */
	  peakvalue = data[i];   
	  apeak = i; 
	  continue;
	}
	if (data[i] > peakvalue) {   /* find largest values and remember where it was */
	  peakvalue = data[i];
	  apeak = i;
	}
      }
    }
    dprintf(1,"apeak=%d\n",apeak);
    if (apeak==0) {             /* never allow first point to be the peak ... */
      mask[apeak] = 0;
      apeak = -1;
      for (i=1, oldvalue=peakvalue; mask[i]==0 && i<n; i++) {   /* walk down, and mask out */
	if (data[i] > oldvalue) break;                         /* util data increase again */
	oldvalue = data[i];
	mask[i] = npeak;
      }
      continue;
    }
    if (apeak==(n-1)) {         /* ...or last point, since they have no neighbors */
      mask[apeak]= 0;
      apeak = -1;
      for (i=n-2; oldvalue=peakvalue, mask[i]==0 && i>0; i--) {  /* walk down, and mask out */
	if (data[i] > oldvalue) break;                        /* until data increase again */
	oldvalue = data[i];
	mask[i] = npeak;
      }
      // break;
      continue;
    }
    if (apeak > 0) {            /* done, found a peak, but mask down the peak */
      mask[apeak] = npeak;
      for (i=apeak+1, oldvalue=peakvalue; mask[i]==0 && i<n;  i++) {
	if (data[i] > oldvalue) break;
	oldvalue = data[i];
	mask[i] = npeak;
      }
      for (i=apeak-1, oldvalue=peakvalue; mask[i]==0 && i>=0; i-- ) {
	if (data[i] > oldvalue) break;
	oldvalue = data[i];
	mask[i] = npeak;
      }
      break;
    }

  } /* while() */
  return apeak;
}

static inline int outside(real x,int i0,int i1) {
  if(x<i0) return 0;
  if(x>i1) return 0;
  return 1;
}

local int peak_assign(int n, real *d, int *s)
{
  int i;
  real p;

  for (i=1; i<n-1; i++) {
    if (s[i] && s[i-1] && s[i]!=s[i-1]) {      /* look for neighbor peaks */
      if (d[i] < d[i-1]) {
	p = peak_spectrum(n,d,i);
	//dprintf(0,"assign-1: %d %d %d   %g   %g %g\n",i,s[i],s[i-1],p,d[i-1],d[i]);
	if (p>0) {
	  dprintf(1,"re-assign-1 @%d from %d to %d\n",i,s[i],s[i-1]);
	  s[i] = s[i-1];                       /* re-assign */
	}
      } else {
	p = peak_spectrum(n,d,i-1);
	//dprintf(0,"assign-2: %d %d %d   %g   %g %g\n",i,s[i],s[i-1],p,d[i-1],d[i]);	
	if (p<0) {
	  dprintf(1,"re-assign-2 @%d from %d to %d\n",i,s[i-1],s[i]);
	  s[i-1] = s[i];                       /* re-assign */
	}
      } 
    }
  }
}


local bool out_of_range(real *clip, real x)
{
  /* already certified that clip has a valid number, and wasn't blank */
  if (clip[0] == clip[1]) return FALSE;
  if (x < clip[0] || x > clip[1]) return FALSE;
  return TRUE;
}

#define CV(p,i,j,k)   CubeValue(p,i,j,k)
#define LOOP(i,n)     for(i=0;i<n;i++)

local void image_oper(imageptr ip1, string oper, imageptr ip2)
{
  int i,j,k, nx,ny,nz;
  nx = Nx(ip1);
  ny = Ny(ip1);
  nz = Nz(ip1);
  if(nx!=Nx(ip2)) error("image_oper: NX size %d != %d",nx,Nx(ip2));
  if(ny!=Ny(ip2)) error("image_oper: NY size %d != %d",ny,Ny(ip2));
  if(nz!=Nz(ip2)) error("image_oper: NZ size %d != %d",nz,Nz(ip2));

  if (*oper== '+')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) + CV(ip2,i,j,k);
  else if (*oper== '-')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) - CV(ip2,i,j,k);
  else if (*oper== '*')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) * CV(ip2,i,j,k);
  else if (*oper== '/')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) / CV(ip2,i,j,k);
  else 
    error("invalid operator: %s",oper);
 
}

