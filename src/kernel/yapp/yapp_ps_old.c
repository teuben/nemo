/*
 * YAPPPS: Yet Another Plotting Package, using postscript
 *    Peter Teuben   March 87   IAS Princeton, NJ
 *	10-Apr-87   some improvements, needs still more hacking   PJT
 *	25-apr-87   bug: maximum number of lineto's before a stroke~1500
 *			THIS IS NOT A BUG, IT IS A FEATURE OF LP !!
 *	 6-May-87   a /s5 macro did not generate a stroke
 *		    header per page + page number added
 *	15-jun-86   a more formal version now, calibrated sizes etc.-system crash
 *	17-jun-87   support for psfig/TeX !!!
 *	25-jun-87   experiments with grayscale fillers
 *	 8-jul-87   proper def. of cell position into pl_matrix()
 *
 * Remarks: to send the file to laserprinter:   'lpr -h yapp.ps'
 *	    to check status of ps-file in case
 *		no output has appeared:		'tail /usr/spool/lw/lw-log'
 *		this has to be done on the laserprinter host machine!!
 */
 
#include <stdio.h>
#include <stdinc.h>		/* nemo */
#include <string.h>
#include <time.h>		/* just to get a time stamp */
#include <sys/types.h>
#include <sys/timeb.h>

#define PS_XMAX   540		/* AppleLaser size in points */
#define PS_YMAX   720

#define PICSCALE (20.0 * 72.0 / 2.54 / 540.0)	/* assuming 1inch=72pt */

#define BUFSIZE  1000		/* seems to have a finite buffer of lineto's */

#define PTPCM	(72.0/2.54)	/* points per cm; for conversions */
#define FUDGE   1.40		/* fudge factor to make pltext() larger by */

/*	Global Static constants, used by the plotting package */


local double	ax,bx,ay,by;			/* transformation parameters */
local double	xp,yp;				/* current pen position	 */
local double	gxmin, gxmax, gymin, gymax;
local double	dxymax,width,height;
local stream	o;				/* ps-file */
local int	vecstate =   0;			/* status of stroke's */
local int	textjust =  -1;
local int	oldps;				/* font size */
local char	yappfile[20]="yapp.ps";		/* ps file name */
local char	font[80]="Times-Roman";		/* default font name */
local char	dateid[80]="";
local char 	versionid[80]="Version 1.3c 17-Jun-87 P.J. Teuben";
local char	*tptr;				/* pointer to time string in ASCII */
local long	lt;				/* local time */
local int	page = 1;
local int	linetos = 0;			/* count for buffer overflow */
local bool	verbose = TRUE;			/* headers? etc */

/*==========================================================================*/
/*
 * PLINIT : begin of all graphics
 */
 
plinit(pltdev, xmin, xmax, ymin, ymax)
char	*pltdev;		/* output device name (ignored) */
double	xmin, xmax;	/* user plotting area */
double	ymin, ymax;
{
	if (strcmp(pltdev,"ps-")==0) 
		verbose=FALSE;			/* silent mode */
	else
		verbose=TRUE;			/* more garbage on plot (header) */
	printf ("YAPP with PostScript; output will be on 'yapp.ps' and can\n");
	printf ("be send to laserprinter using 'lpr -h yapp.ps'\n");
	printf ("`tail /usr/spool/lw/lw-log` on lp-hostmachine for status\n\n");
	if (verbose)
		printf ("Verbose is on (default)\n");
	else
		printf ("Verbose is off\n");
	printf ("%s\n",versionid);


	lt = time(0);		 /* get the seconds since 1970 */
	tptr=ctime(&lt);	 /* and convert to ASCII */
	strcat (dateid,tptr);	 /* put it in a string for plot-header */

	printf ("Dateid=%s\n",dateid);		

		
	if ((o=fopen(yappfile,"w"))==NULL) {
		printf ("Error: %s cannot be opened\n",yappfile);
		exit(1);
	}

	vec_open();		/* initialize yapp.ps file */
	
				/******   initialize screen   *****/
	gxmin = 0.0;
	gxmax = PS_XMAX;
	gymin = 0.0;
	gymax = PS_YMAX;

	ax = (gxmax-gxmin)/(xmax-xmin);
	bx = gxmin - ax*xmin;
	ay = (gymax-gymin)/(ymax-ymin);
	by = gymin - ay*ymin;

	if (ymax-ymin < xmax-xmin) {
		dxymax = xmax - xmin;
		width = 1.0;
		height= (ymax-ymin) / dxymax;
	}
	else {
		dxymax = ymax - ymin;
		width = (xmax-xmin) / dxymax;
		height= 1.0;
	}
}


/*
 * PLFRAME : 	defines new page/screen
 */
 
plframe()
{
	fprintf(o,"showpage\n");
	page++;
	if (verbose) {
		fprintf(o," 76 685 mt (Page - %d) show\n",page);
		fprintf(o,"548 685 mt (%s) rshow\n",dateid);	
	}
}


/*
 * PLFLUSH: output pending graphics
 */
 
plflush()
{
	/*	does nothing yet */
} 


/*
 * PLSTOP : end of all graphics
 */
 
plstop()
{
	vec_close();
}


/*
 * PLLINE, PLMOVE, PLPOINT
 */
 
plline(x,y)
double x,y;
{
	double convx(), convy();
	
	fprintf(o,"%.1f %.1f lt\n",convx(x),convy(y));
	vecstate = 2;
	if (linetos++ > BUFSIZE) {		/* prevent buffer overflow */
		plmove (x,y);
		linetos = 0;			/* could be done in vec_stroke? */
	}
}

plmove(x,y)
double x,y;
{

	double convx(),convy();
	
	if (vecstate==2) vec_stroke();
/*	if (vecstate==2 || vecstate==0) fprintf(o,"newpath\n");   */
	fprintf(o,"%.1f %.1f mt\n",convx(x),convy(y));
	vecstate = 1;
}
	
plpoint(x,y)
double x,y;
{
	double convx(),convy();
	
	if (vecstate==2) vec_stroke();
	fprintf (o,"%.1f %.1f s5\n",convx(x), convy(y));
	vecstate = 1;
}


/*
 *  PLCIRCLE, PLCROSS, PLBOX:
 */
 
plcircle (x,y,r)
double x,y,r;
{
	int	npnts, i;
	double	theta, sin(), cos();

/*  Old PS independent solution --------------------------- */
	npnts = MAX(2400 * r / dxymax, 6.0);
	plmove (x+r,y);
	for (i=1; i<=npnts; i++) {
		theta = TWO_PI * ((double) i)/((double) npnts);
		plline (x + r*cos(theta), y+r*sin(theta));
	}
 
/* PS only */
/*	fprintf(o,"%.1f %.1f {2 copy newpath %.1f 0 360 arc closepath stroke}\n",
		x, y, 10.0*r);
   DOES not work yet in this way
   */
}


plcross(x,y,s)
double x,y,s;
{
	if (s>0.0) {
		plmove (x-s,y);
		plline (x+s,y);
		plmove (x,y-s);
		plline (x,y+s);
	} else {
		s = s / 1.4142;
		plmove (x-s,y-s);
		plline (x+s,y+s);
		plmove (x-s,y+s);
		plline (x+s,y-s);
	}
}


plbox (x,y,s)
double x,y,s;
{
	if (s>0.0) {
		plmove (x-s,y-s);
		plline (x+s,y-s);
		plline (x+s,y+s);
		plline (x-s,y+s);
		plline (x-s,y-s);
	} else {
		s = s / 1.4142;
		plmove (x-s,y);
		plline (x,y-s);
		plline (x+s,y);
		plline (x,y+s);
		plline (x-s,y);
	}
}

/*
 * PLTEXT : 
 */
 
pltext (msg, x, y, hgt, ang)
char *msg;
double x,y;
double hgt;
double ang;
{
	int fontsize;
	
	fontsize = hgt * PTPCM * FUDGE;
	vec_size ( fontsize );		/* change font size */

	fprintf(o,"gsave\n");		/* save graphics state, because we will
					   use a funny translation now */
	fprintf(o,"%.1f %.1f translate\n",convx(x),convy(y));
	fprintf(o,"%.1f rotate newpath\n",ang);
	fprintf(o,"0 0 moveto\n");

	if (textjust==-1)
		fprintf(o,"(%s) lshow\n",msg);
	else if (textjust==1)
		fprintf(o,"(%s) rshow\n",msg);
	else
		fprintf(o,"(%s) cshow\n",msg);

	fprintf(o,"grestore\n");		/* restore graphics state */
}

/*
 * PLJUST:
 */
 
pljust(jus)
int jus;
{	
	textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}


/*
 * PLLTYPE: linetype: thickness and line pattern dash/dot/solid
 */
 
plltype (lwid, lpat)
int lwid;			/* line width */
int lpat;			/* line pattern */
{

   if (vecstate==2) vec_stroke();

   if (lwid>0) {
   	fprintf(o,"%.1f setlinewidth\n",0.48*(double)lwid);
   }

   if (lpat>0) {
	fprintf(o,"[");
	switch(lpat) {
		case 0:
		case 1: break;			/* solid */
		case 4: fprintf(o,"9 3");	/* long dashed */
			break;
		case 2: fprintf(o,"3");		/* dotted */
			break;
		case 3: fprintf(o,"6 3");	/* short dashed */
			break;
		case 6:	fprintf(o,"3 3 6 3");	/* dotdashed */
			break;
	}
	fprintf(o,"] 0 setdash\n");
   }
   vecstate=1;
}


/*
 * PLSWAP: does nothing
 */
 
plswap()
{
}

/*
 * PLXSCALE, PLYSCALE:  
 *			at the moment they are identity transformations
 */
 
double plxscale(x,y)
double x,y;
{
	return(x);	/* nothing */
}

double plyscale(x,y)
double x,y;
{
	return(y);	/* nothing */
}

/*
 * PL_MATRIX:    paint a matrix:   0=white 1=black (if fmin<fmax)
 *				  0=black 1=white (if fmin>fmax) = normal PS def'n
 *		I.e. normally a negative image --- hmm, maybe we should convert again
 *	frame = matrix (1D array)    (  f[ix][iy] == *(f+ix*ny+iy)  )
 *	nx,ny = size of matrix	     ( must be contguous in memory due to above definition )
 *	xmin,ymin = position (in cm) of middle of lower left matrix cell
 *	cell = cellsize (cm)
 *	fmin = what should be white (if fmin<fmax)
 *	fmax = what should be black 
 *	findex = power of transfer function (0,infinity) 1=linear
 */

pl_matrix (frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex,blank)
double *frame, xmin, ymin, cell, fmin, fmax, findex, blank;
int nx, ny;
{
	double x,y,f,grayscale,ds,pow();
	int    ix,iy;

	ds=cell*PTPCM;			/* convert cm -> pixels */

	printf ("GRAY_MATRIX cm: %d * %d matrix  ll corner: %f %f cell %f\n",nx,ny,xmin,ymin,cell);
	printf ("GRAY_MATRIX DC: %d * %d matrix  ll corner: %f %f cell %f\n",nx,ny,
		convx(xmin),convy(ymin),ds);	
	if (fmin==fmax)
		error("PL_MATRIX: equal min and max");
			
	grayscale = 1.0/(fmax-fmin);	/* normally positive for a negative image */
	

	for (ix=0, x=xmin-0.5*cell; ix<nx; ix++, x+=cell) {
	   for (iy=0, y=ymin-0.5*cell; iy<ny; iy++, y+=cell) {
		f = *(frame + ix*ny + iy);
		
					/* apply linear grayscale + cutoff */
		if (grayscale>0.0)
			f = (f-fmin)*grayscale;
		else
			f = (fmin-f)*grayscale;
		if (f>1.0)			/* upper cutoff */
			f=1.0;
		if (f<0.0)			/* lower cutoff */
			f=0.0;
		f = pow(f,findex);		/* transfer function */
		
		vec_paint (convx(x),convy(y),ds,f);		/* paint this square */
	   }  
	}
}		


vec_paint (x,y,dx,gray)
double x,y,dx,gray;		/* coordinates must be in DC (x:0-540  y:0-720) */
{
	char line[80];

	if (gray<0.01)		/* what will be white, can be ignored */
		return(0);

	dx += 0.5;		/* kludge to hopefully solve the interference */
		
	sprintf (line,"%.1f %.0f %.0f %.2f f",dx,x,y,1.0-gray);
	vec_write (line);
}



/*
 *  some lower level PS related things 
 *
 */
 
vec_open()
{
	fprintf(o,"%%!PSAdobe-1.0\n");
	fprintf(o,"%%%%BoundingBox: 0 0 %d %d\n",PS_XMAX,PS_YMAX);	/* for TEX/psfig */
	fprintf(o,"%%%%Pages:1\n");
	fprintf(o,"%%%%Creator: Peter Teuben\n");
	fprintf(o,"%%%%Title: YAPP.PS\n");
	fprintf(o,"%%%%DocumentsFonts Times-Bold\n");
	fprintf(o,"%%%%EndComments\n");
	fprintf(o,"%%%%EndProlog\n");
	fprintf(o,"%%%%Page: 0 1\n");
	fprintf(o,"%%\n");
	fprintf(o,"%% some aliases \n");
	fprintf(o,"/off {36 add} def\n");	/* actually not used anymore */
	fprintf(o,"/mt {moveto} def\n");
	fprintf(o,"/lt {lineto} def\n");
	fprintf(o,"/s {gsave stroke grestore newpath} def\n");
	fprintf(o,"/lshow {5 0 8#040 4 3 roll widthshow} def\n");
	fprintf(o,"/cshow {dup stringwidth pop 2 div neg 0 rmoveto show} def\n");
	fprintf(o,"/rshow {dup stringwidth pop neg 0 rmoveto show} def\n");
	fprintf(o,"/Times-Roman findfont 12 scalefont setfont\n");
	fprintf(o,"/starbody {newpath 0 360 arc closepath fill} def \n");
	fprintf(o,"/s5 {0.33 starbody} def \n");
	fprintf(o,"/box {newpath 0 0 mt 0 1 lt 1 1 lt 1 0 lt closepath} def\n");    /* filling */
	fprintf(o,"/f {gsave setgray translate dup scale box fill grestore} def\n"); 
#ifdef TESTBED
	vec_size(18);					
	fprintf(o,"314 685 mt (YAPP.PS) cshow\n");	
#endif
	vec_size(10);
#ifdef TESTBED
	fprintf(o," 76 685 mt (%s) show\n",versionid);	
#else
	if (verbose)
		fprintf(o," 76 685 mt (Page - %d) show\n",page);
#endif
	if (verbose)
		fprintf(o,"548 685 mt (%s) rshow\n",dateid);	
	fprintf(o,"0.5 setlinewidth\n");
	fprintf(o,"newpath\n");	

}

vec_close()
{
	fprintf(o,"stroke\n");
	fprintf(o,"showpage\n");
	fprintf(o,"%% end of yapp.ps\n");

	fclose(o);

}

vec_size(newps)
int newps;
{
	if (newps!=oldps)
		fprintf (o,"/%s findfont %d scalefont setfont\n", font, newps);
	oldps = newps;
#ifdef DEBUG
	printf ("vec_font size=%d\n",newps);
#endif	
}

vec_write(command)
char *command;
{
	fprintf(o,"%s\n",command);
}

vec_stroke()
{
	fprintf(o,"s\n");
}

vec_font(fontname)
char *fontname;
{
	strcpy (font,fontname);
	fprintf(o,"/%s findfont %d scalefont setfont\n",fontname,oldps);
}


/* conversion from xm -> pixels */

double convx(x)
double x;
{
	x = ax*x + bx;		/* temp sol: x and y same */
	return( x*PICSCALE+25);
}

double convy(y)
double y;
{
	y = ax*y + bx;		/* temp sol: x and y same to be square */
	return( y*PICSCALE+100);	/* about to center the 20*20 on page */
}


#ifdef TESTBED

main(argc, argv)
int argc;
string argv[];
{
    int i, j;
    double x,y,dx,dy,gray;

    
    plinit("***", 0.0, 20.0, 0.0, 20.0);
    plmove(0.0, 0.0);
    plline(20.0, 0.0);
    plline(20.0, 20.0);
    plline(0.0, 20.0);
    plline(0.0, 0.0);
    plline(20.0, 20.0);
    plmove(20.0, 0.0);
    plline(0.0, 20.0);
    plltype(12, 0);
    plmove(4.0, 18.0);
    plline(16.0, 18.0);
    for (i = 1; i <= 4; i++) {
	plltype(2*i, 1);
        plmove(1.0, 13.0 - i);
        plline(3.0, 13.0 - i);
        plpoint(3.5, 13.0 - i);
	plltype(1, i);
	for (j = 1; j <= 4; j++) {
	    plmove(1.5, 13.0 - i - 0.2*j);
	    plline(1.5 + j, 13.0 - i - 0.2*j);
	}
    }
    plltype(1, 1);
    plcircle(15.0, 9.0, 0.5);
    plcircle(16.0, 9.0, 0.25);
    plcircle(17.0, 9.0, 0.125);
    plcircle(18.0, 9.0, 0.0625);
    plbox(16.0, 8.0, 0.4);
    plbox(17.0, 8.0, 0.2);
    plbox(18.0, 8.0, -0.2);
    plcross(16.0, 7.0, 0.4);
    plcross(17.0, 7.0, 0.2);
    plcross(18.0, 7.0, -0.2);
    pltext("Foo Bar!", 8.0, 5.0, 0.5, 0.0);
    pltext("Fum Bar!", 8.0, 3.0, 0.25, 0.0);
    for (i = 0; i <= 4; i++)
	pltext(" testing angles", 16.0, 10.0, 0.2, 45.0*i);
    plmove(10.0, 8.5);
    plline(10.0, 11.5);
    pljust(-1);
    pltext("left justified",  10.0,  9.0, 0.25, 0.0);
    pljust(0);
    pltext("centered",        10.0, 10.0, 0.25, 0.0);
    pljust(1);
    pltext("right justified", 10.0, 11.0, 0.25, 0.0);

    plframe();
    plmove(0.0, 0.0);
    plline(20.0, 0.0);
    plline(20.0, 20.0);
    plline(0.0, 20.0);
    plline(0.0, 0.0);

    pljust(0);
    pltext("This is page 2", 10.0,10.0,0.25,0.0);

#define IMAX 100
#define ISTEP 20.0/IMAX

    plmove (0.0,0.0);
    for (i=0; i<IMAX; i++)
       plline(i*ISTEP, i*ISTEP);

#define NCOLORS 32
    dx = (PS_XMAX-50)/ NCOLORS;	/* use offset */
    printf ("Stepsize = %f for %d colors\n",dx,NCOLORS);
    for (i=0; i<NCOLORS; i++) {
    	x=i*dx;
    	gray = i/(NCOLORS-1.0);	
    	vec_paint (x+50.0,300.0,dx,gray);
    }

    plstop();
}

#endif
