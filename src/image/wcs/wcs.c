/*
 * WCS: wrapper to worldpos/xypix WCS
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>

string defv[] = {
    "ref=\n     Reference coordinates (deg)",
    "refpix=\n  Reference pixels",
    "inc=\n     Reference increment",
    "rot=\n     Rotation (deg) (from N through E)",
    "type=\n    Projection type (SIN, TAN, ARC, NCP, GLS, MER, AIT)",
    "format=%g\n Output format",
    "pix=\n     (optional) Input Pixel number coordinates",
    "pos=\n     (optional) Input Coordinates",
    "VERSION=1.0\n  14-oct-94",
    NULL,
};

string usage = "coordinate transforms";

nemo_main() 
{
    double pix[2], pos[2], ref[2], refpix[2], inc[2], rot, *outval;
    char proj[10];
    bool Qworldpos;
    string type, fmt = getparam("format");
    int  n;
    

    if (!hasvalue("pix") && !hasvalue("pos")) error("Need pix= or pos=");
    if (hasvalue("pix") && hasvalue("pos")) error("Need one of pix= or pos=");

    Qworldpos = hasvalue("pix");
    if (Qworldpos)
        n = nemoinpd(getparam("pix"),pix,2);
    else
        n = nemoinpd(getparam("pos"),pos,2);
    if (n!=2) error("%d: error parsing pix/pos=",n);

    n = nemoinpd(getparam("ref"),ref,2);
    if (n!=2) error("%d: parsing ref",n);
    n = nemoinpd(getparam("refpix"),refpix,2);
    if (n!=2) error("%d: parsing refpix",n);
    n = nemoinpd(getparam("inc"),inc,2);
    if (n!=2) error("%d: parsing inc",n);
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
    printf("\n");
   
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
    return worldpos(pos[0],pos[1], ref[0],ref[1], refpix[0],refpix[1],
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
    

#include "worldpos.c"

