/*
 * YAPP_PSG: Yet Another Plotting Package, using postscript w/ greyscales
 *
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
 *      30-mar-97   2.0 resurrection from yapp_ps_old.c since nothing else 
 *			has greyscales
 *		    greyscale is ok for SINGLEPREC, but contours not !!!
 *      16-mar-05   added blankval for pl_matrix
 *		    
 *
 */
 
#include <stdinc.h>
#include <yapp.h>

#include <time.h>		/* to get a time stamp */
#include <sys/types.h>
#if 0
	/* was ok for old BSD systems, but SYSV doesn't have/need this */
#include <sys/timeb.h>		/* e.g. does not exist on SGI */
#endif

#define PS_XMAX   540		/* AppleLaser size in points */
#define PS_YMAX   720

#define PICSCALE (20.0 * 72.0 / 2.54 / 540.0)	/* assuming 1inch=72pt */

#define BUFSIZE  1000		/* seems to have a finite buffer of lineto's */

#define PTPCM	(72.0/2.54)	/* points per cm; for conversions */
#define FUDGE   1.40		/* fudge factor to make pltext() larger by */

/*	Global Static constants, used by the plotting package */


local real	ax,bx,ay,by;			/* transformation parameters */
local real	xp,yp;				/* current pen position	 */
local real	gxmin, gxmax, gymin, gymax;
local real	dxymax,width,height;
local stream	o;				/* ps-file */
local int	vecstate =   0;			/* status of stroke's */
local int	textjust =  -1;
local int	oldps;				/* font size */
local char	yappfile[20]="yapp.ps";		/* ps file name */
local char	font[80]="Times-Roman";		/* default font name */
local char	dateid[80]="";
local char 	versionid[80]="Version 2.0 30-mar-97 P.J. Teuben";
local char	*tptr;				/*  -> time string in ASCII */
local long	lt;				/* local time */
local int	page = 1;
local int	linetos = 0;			/* count for buffer overflow */
local bool	verbose = TRUE;			/* headers? etc */

local real convx(real), convy(real);
local void vec_paint(real,real,real,real);
local void vec_open(void);
local void vec_close(void);
local void vec_size(int);
local void vec_write(string);
local void vec_stroke(void);
local void vec_font(string);


/*==========================================================================*/
/*
 * PLINIT : begin of all graphics
 */
 
plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
	if (strcmp(pltdev,"ps-")==0) 
		verbose=FALSE;			/* silent mode */
	else
		verbose=TRUE;			/* more garbage on plot (header) */
	dprintf(0,"YAPP with PostScript; output will be on 'yapp.ps'\n");
	dprintf(0,"%s\n",versionid);


	lt = time(0);		 /* get the seconds since 1970 */
	tptr=ctime(&lt);	 /* and convert to ASCII */
	strcat (dateid,tptr);	 /* put it in a string for plot-header */

	dprintf(1,"Dateid=%s\n",dateid);		

		
	if ((o=fopen(yappfile,"w"))==NULL)
		error("%s cannot be opened",yappfile);

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
 
plmove(real x, real y)
{
	if (vecstate==2) vec_stroke();
/*	if (vecstate==2 || vecstate==0) fprintf(o,"newpath\n");   */
	fprintf(o,"%.1f %.1f mt\n",convx(x),convy(y));
	vecstate = 1;
}
	
plline(real x, real y)
{
	fprintf(o,"%.1f %.1f lt\n",convx(x),convy(y));
	vecstate = 2;
	if (linetos++ > BUFSIZE) {		/* prevent buffer overflow */
		plmove (x,y);
		linetos = 0;			/* could be done in vec_stroke? */
	}
}

plpoint(real x, real y)
{
	if (vecstate==2) vec_stroke();
	fprintf (o,"%.1f %.1f s5\n",convx(x), convy(y));
	vecstate = 1;
}


/*
 *  PLCIRCLE, PLCROSS, PLBOX:
 */
 
plcircle (real x, real y, real r)
{
	int	npnts, i;
	real	theta;

/*  Old PS independent solution --------------------------- */
	npnts = MAX(2400 * r / dxymax, 6.0);
	plmove (x+r,y);
	for (i=1; i<=npnts; i++) {
		theta = TWO_PI * ((real) i)/((real) npnts);
		plline (x + r*cos(theta), y+r*sin(theta));
	}
 
/* PS only */
/*	fprintf(o,"%.1f %.1f {2 copy newpath %.1f 0 360 arc closepath stroke}\n",
		x, y, 10.0*r);
   DOES not work yet in this way
   */
}


plcross(real x, real y, real s)
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


plbox (real x, real y, real s)
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
 
pltext (string msg, real x, real y, real hgt, real ang)
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
 
pljust(int jus)
{	
	textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}


/*
 * PLLTYPE: linetype: thickness and line pattern dash/dot/solid
 */
 
plltype (int lwid, int lpat)
{
   if (vecstate==2) vec_stroke();

   if (lwid>0) {
   	fprintf(o,"%.1f setlinewidth\n",0.48*(real)lwid);
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
 
real plxscale(real x, real y)
{
	return x;
}

real plyscale(real x, real y)
{
	return y;
}

pl_screendump(string fname)
{
	error("pl_screendump: not implemented");
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
 *      blankval = value where no paint is applied
 */

pl_matrix (real *frame, int nx, int ny,
           real xmin, real ymin, real cell, real fmin, real fmax, real findex, real blankval)
{
	real x,y,f,grayscale,ds;
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
		if (f==blankval) continue;
		
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


local void vec_paint (real x, real y, real dx, real gray)
		/* coordinates must be in DC (x:0-540  y:0-720) */
{
	char line[80];

	if (gray<0.01)		/* what will be white, can be ignored */
		return;

	dx += 0.5;		/* kludge to hopefully solve the interference */
		
	sprintf (line,"%.1f %.0f %.0f %.2f f",dx,x,y,1.0-gray);
	vec_write (line);
}



/*
 *  some lower level PS related things 
 *
 */

local void vec_open()
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

local void vec_close()
{
	fprintf(o,"stroke\n");
	fprintf(o,"showpage\n");
	fprintf(o,"%% end of yapp.ps\n");

	fclose(o);

}

local void vec_size(int newps)
{
	if (newps!=oldps)
		fprintf (o,"/%s findfont %d scalefont setfont\n", font, newps);
	oldps = newps;
#ifdef DEBUG
	printf ("vec_font size=%d\n",newps);
#endif	
}

local void vec_write(string command)
{
	fprintf(o,"%s\n",command);
}

local void vec_stroke()
{
	fprintf(o,"s\n");
}

local void vec_font(string fontname)
{
	strcpy (font,fontname);
	fprintf(o,"/%s findfont %d scalefont setfont\n",fontname,oldps);
}


/* conversion from xm -> pixels */

local real convx(real x)
{
	x = ax*x + bx;		/* x and y same */
	return x*PICSCALE+25;
}

local real convy(real y)
{
	y = ax*y + bx;		/* x and y same to be square */
	return y*PICSCALE+100;	/* to center the 20*20 cm on a page */
}

