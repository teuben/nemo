/*
 * WCS: wrapper to worldpos/xypix WCS
 *
 *	2-mar-03	1.1	fixed pos->pix- conversion
 *      7-may-04        1.3     added hms=
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>

string defv[] = {
    "ref=\n     Reference coordinates (deg) [i.e. crval]",
    "refpix=\n  Reference pixels [i.e. crpix]",
    "inc=\n     Reference increment [i.e. cdelt]",
    "rot=\n     Rotation (deg) (from N through E)",
    "type=\n    Projection type (SIN, TAN, ARC, NCP, GLS, MER, AIT)",
    "format=%g\n Output format",
    "pix=\n     (optional) Input Pixel number coordinates",
    "pos=\n     (optional) Input Coordinates",
    "hms=t\n    Should Longitude in HMS (else DMS) ?",
    "VERSION=1.3\n  7-may-04",
    NULL,
};

string usage = "coordinate transforms";

bool is_sexa(string);

nemo_main() 
{
    double pix[2], pos[2], ref[2], refpix[2], inc[2], rot, ss, *outval;
    char proj[10];
    bool Qworldpos, Qhms;
    string type, fmt = getparam("format"), sexa;
    int  n, hh,dd,mm;

    Qhms = getbparam("hms");
    

    if (!hasvalue("pix") && !hasvalue("pos")) error("Need pix= or pos=");
    if (hasvalue("pix") && hasvalue("pos")) error("Need one of pix= or pos=");

    Qworldpos = hasvalue("pix");
    if (Qworldpos) {
        n = nemoinpd(getparam("pix"),pix,2);
	dprintf(1,"pix: %g %g\n",pix[0],pix[1]);
    } else {
	if (is_sexa(getparam("pos"))) {
	  n = nemoinpx(getparam("pos"),pos,2);
	  pos[0] *= 15;   /* assume that was hms :-) */
	} else {
	  n = nemoinpd(getparam("pos"),pos,2);
	}
	dprintf(1,"pos: %g %g\n",pos[0],pos[1]);
    }
    if (n!=2) error("%d: error parsing pix/pos=",n);

    
    if (is_sexa(getparam("ref"))) {
      n = nemoinpx(getparam("ref"),ref,2);
      ref[0] *= 15;   /* assume that was hms :-) */
    } else {
      n = nemoinpd(getparam("ref"),ref,2);
    }
    if (n!=2) error("%d: parsing ref",n);
    n = nemoinpd(getparam("refpix"),refpix,2);
    if (n!=2) error("%d: parsing refpix",n);
    
    if (is_sexa(getparam("inc"))) {
      n = nemoinpx(getparam("inc"),inc,2);
      ref[0] *= 15;   /* assume that was hms :-) */
    } else {
      n = nemoinpd(getparam("inc"),inc,2);
    }
    if (n!=2) error("%d: parsing inc",n);
    dprintf(1,"ref: %g %g  refpix: %g %g  inc: %g %g\n",
	    ref[0],ref[1],refpix[0],refpix[1],inc[0],inc[1]);

    rot = getdparam("rot");
    type = getparam("type");
    if (*type == '-') 
        strcpy(proj,type);
    else
        sprintf(proj,"-%s",type);
    strtoupper(proj);
    if (Qworldpos) {
        n = my_worldpos(pix,ref,refpix,inc,rot,proj,pos);
        outval = pos;
    } else {
        n = my_xypix(pos,ref,refpix,inc,rot,proj,pix);
        outval = pix;
    }
    if (n!=0) error("%d: returned from worldpos/xypix",n);
    printf(fmt,outval[0]);
    printf(" ");
    printf(fmt,outval[1]);
    printf(" ");
    if (Qworldpos) {
      if (Qhms) {
	to_hms(outval[0],&hh,&mm,&ss);
	printf("%3d:%02d:%06.3f ",hh,mm,ss);
      } else {
	to_dms(outval[0],&hh,&mm,&ss);
	printf("%03d:%02d:%06.3f ",hh,mm,ss);
      }
      to_dms(outval[1],&dd,&mm,&ss);
      printf("%03d:%02d:%06.3f",dd,mm,ss);
    }
    printf("\n");
}

to_hms(double dval, int *dd, int *mm, double *ss)
{
  int sign = SGN(dval);
  dval = ABS(dval)/15.0;
  *dd = (int) floor(dval);
  dval = (dval-(*dd))*60.0;
  *mm = (int) floor(dval);
  *ss = (dval-(*mm))*60.0;
  *dd *= sign;
}

to_dms(double dval, int *dd, int *mm, double *ss)
{
  int sign = SGN(dval);
  dval = ABS(dval);
  *dd = (int) floor(dval);
  dval = (dval-(*dd))*60.0;
  *mm = (int) floor(dval);
  *ss = (dval-(*mm))*60.0;
  *dd *= sign;
}




/*
 * 	interface to the Public Domain worldpos/xypix routines
 */

int worldpos(double xpix, double ypix, double xref, double yref,
      double xrefpix, double yrefpix, double xinc, double yinc, double rot,
      char *type, double *xpos, double *ypos);

int xypix(double xpos, double ypos, double xref, double yref, 
      double xrefpix, double yrefpix, double xinc, double yinc, double rot,
      char *type, double *xpix, double *ypix);


int my_worldpos(double *pix,double *ref,double *refpix,double *inc,
            double rot,string proj,double *pos)
{
    return worldpos(pix[0],pix[1], ref[0],ref[1], refpix[0],refpix[1],
        inc[0],inc[1], rot, proj, &pos[0], &pos[1]);
}

int my_xypix(double *pos,double *ref,double *refpix,double *inc,
            double rot,string proj,double *pix)
{
    return xypix(pos[0],pos[1], ref[0],ref[1], refpix[0],refpix[1],
        inc[0],inc[1], rot, proj, &pix[0], &pix[1]);
}

strtoupper(char *text)
{
    char *cp = text;
    while (*cp) {
        if (islower(*cp)) *cp = toupper(*cp);
        cp++;
    }
}

bool is_sexa(string s) 
{
  char *cp = strchr(s, ':');
  return (cp != 0);
}
    

#ifdef NEMO
#include "worldpos.c"
#endif
