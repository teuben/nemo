/*
    testpos.c: Program to test the two functions in file worldpos.c 
    Copyright (C) 1994
    Associated Universities, Inc. Washington DC, USA.
   
    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Library General Public License as published by
    the Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.
   
    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
    License for more details.
   
    You should have received a copy of the GNU Library General Public License
    along with this library; if not, write to the Free Software Foundation,
    Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
   
    Correspondence concerning AIPS should be addressed as follows:
           Internet email: aipsmail@nrao.edu
           Postal address: AIPS Group
                           National Radio Astronomy Observatory
                           520 Edgemont Road
                           Charlottesville, VA 22903-2475 USA
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

int worldpos(double xpix, double ypix, double xref, double yref,
      double xrefpix, double yrefpix, double xinc, double yinc, double rot,
      char *type, double *xpos, double *ypos);

int xypix(double xpos, double ypos, double xref, double yref, 
      double xrefpix, double yrefpix, double xinc, double yinc, double rot,
      char *type, double *xpix, double *ypix);

main() 
{
    double xpix, ypix;
    double xpos, ypos;
    double xref, yref, xrefpix, yrefpix, xinc, yinc, rot;
    char ctypes[9][5] ={"???","-SIN","-TAN","-ARC","-NCP", "-GLS", "-MER", 
			    "-AIT", "-STG"};
    char proj[5], msg[100], temp[60];
    int i, j, k, allsky, s_xypix, s_worldpos;
    double xpos2, ypos2;
    double dxpos, dypos, tolerance;
    
    printf ("\n    Circularity tests of various projection types\n");
    printf ("        (output from program \'testpos\')\n");
    printf ("\n    xypix() does (RA,Dec)-->(x,y) transformation,\n");
    printf ("    worldpos() does (x,y)-->(RA2,Dec2),\n");
    printf ("    then we compare (RA2,Dec2) with (RA,Dec).\n");

    for (i = 0; i < 2; i++) {
	for (k = 0; k < 3; k++) {
	    switch (i) {
	    case 0: /* small-field */
		switch (k) {
		case 0:
		    xref = 30.0; yref = 40.0; 
		    xpos = xref + 0.4; ypos = yref + 0.3;
		    break;
		case 1:
		    xref = yref = 0.0;
		    printf ("\nNOTE: NCP error below is correct behavior!\n");
		    printf ("        (NCP undefined for yref==0)\n");
		    xpos = 359.999 - 0.4; ypos = yref - 0.3;
		    break;
		case 2:
		    xref = 359.8; yref = 40.0;
		    xpos = 360.2; ypos = yref + 0.3;
		    break;
		default:
		    continue;
		}
		xrefpix = yrefpix = 500.0;
		xinc = yinc = 0.001;
		rot = -13.0;
		printf ("\n    Small-field cases:\n\n");
		break;
	    case 1: /* all-sky */
		xref = yref = 0.0;
		switch (k) {
		case 0:
		    xpos = 317.0; ypos = -71.0;
		    break;
		case 1:
		    xpos = 181.0; ypos = 87.0;
		    break;
		case 2:
		    xpos = 181.0; ypos = 67.0;
		    break;
		default:
		    continue;
		}
		xrefpix = yrefpix = 500.0;
		xinc = yinc = 0.5;
		rot = 0.0; 
		printf ("\n    All-sky cases:\n\n");
		break;
	    default: 
		exit(13);
	    }
	    
	    printf ("   pos=%6.1lf %6.1lf deg\n", xpos, ypos);
	    printf ("   ref=%6.1lf %6.1lf deg\n", xref, yref);
	    printf ("refpix=%6.1lf %6.1lf pix\n", xrefpix, yrefpix);
	    printf ("   inc=%6.3lf %6.3lf deg\n", xinc, yinc);
	    printf ("   rot=       %6.1lf deg\n", rot);
	    printf ("\n");
	    
	    tolerance = 1e-9 / (3.14159265 / 180.0); /* 1 nanoradian */
	    printf ("i  code    xpix     ypix     Remarks (tolerance=%.1lg_deg)\n",
		    tolerance);
	    printf ("-  ----   ------   ------   ---------\n");
	    
	    for (j = 0; j < 9; j++) {
		allsky = ((j==6) || (j==7));
		if (((i==0) && allsky) 
		    || ((i==1) && !allsky)) continue;
		strcpy (proj, ctypes[j]);
		s_xypix = xypix (xpos, ypos, 
				 xref, yref, xrefpix, yrefpix, 
				 xinc, yinc, rot, proj,
				 &xpix, &ypix);
		s_worldpos = worldpos (xpix, ypix,
				       xref, yref, xrefpix, yrefpix, 
				       xinc, yinc, rot, proj,
				       &xpos2, &ypos2);
		dxpos = xpos2 - xpos;
		dypos = ypos2 - ypos;
		
		strcpy (msg, "");
		if (s_xypix != 0) {
		    sprintf (temp, " s_xypix=%d;", s_xypix);
		    strcat (msg, temp);
		}
		if (s_worldpos != 0) {
		    sprintf (temp, " s_worldpos=%d", s_worldpos);
		    strcat (msg, temp);
		}
		if (fabs(dxpos) > tolerance) {
		    sprintf (temp, " DX=%lg_deg (xpos2=%lg);", dxpos, xpos2);
		    strcat (msg, temp);
		}
		if (fabs(dypos) > tolerance) {
		    sprintf (temp, " DY=%lg_deg;", dypos);
		    strcat (msg, temp);
		}
		if (strcmp (msg, "") == 0) strcat (msg, " OK");
		if (strcmp(proj,"???")==0) strcat (msg, " (default linear)");
		printf ("%1d %5s %8.3lf %8.3lf  %s\n",
			j, proj, xpix, ypix, msg);
	    }
	}
    }
    exit(EXIT_SUCCESS);
}


#ifdef NEMO
#include "worldpos.c"
#endif

