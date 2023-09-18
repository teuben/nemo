/*
 * YAPP: Yet Another Plotting Package.
 *
 *      15-jun-2018:    cloned off yapp_pgplot                PJT
 *      26-oct-2022:    pure giza now, finally                PJT
 *
 *  Linking will need:    -lgiza
 *  See also "mknemo giza" for a self-built version of giza
 *
 *  @todo    placement of numbers at axis tickmarks seems a bit off
 *           pl_matrix doesn't work
 *           pl_contour doesn't work (but contour does, so ccdplot is ok)
 *           builtin circle doesn't work, but giza_circle is ok
 */

#include <stdinc.h>
#include <yapp.h>
#include <giza.h>

extern string yapp_string;	/* a kludge, see: getparam.c */
extern int debug_level;         /* see dprintf.c   from DEBUG env.var. */

local real   dxymax;    /* size of user window */
local int    iterm;     /* terminal number */

#define MAXCOLOR 256
local int ncolors=0;
local float red[MAXCOLOR];		/* RGB color tables		    */
local float blue[MAXCOLOR];
local float green[MAXCOLOR];
local void cms_rgbsetup();

/*
 * PLINIT: initalize the plotting package.
 */

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
    double width, height, x1,x2,y1,y2, zero, one;
    int   dummy, nx, ny, ask, units;

    warning("testing YAPP_GIZA");

    iterm = 1;
    if (yapp_string == NULL || *yapp_string == 0) {
	if (pltdev == NULL || *pltdev == 0)
        	iterm = 0;
	else
	    yapp_string = pltdev;
    }

    //  define a 20x20 cm grid
    int gid = giza_open_device_size (yapp_string, "?", 200.0, 200.0, GIZA_UNITS_MM);

    int interactive = giza_device_has_cursor();

    // normally 0..20 and 0..20
    giza_set_window_equal_scale(xmin, xmax, ymin, ymax);     
    giza_set_character_height(0.75);
    // use the whole sheet
    giza_set_viewport(0.0, 1.0, 0.0, 1.0);

    // dummy frame setting
    giza_print_id();
    //giza_label("Xlabel", "Ylabel", "PlotTitle");


    char font[128];
    giza_get_font(font, 128);
    dprintf(0,"font: %s\n", font);
    strcpy(font,"Sans Serif");
    strcpy(font,"Arial");
    giza_set_font(font);
    giza_set_fill(2);   // 1=solid 2=hollow 3=hatch 4=crosshatch
   
    cms_rgbsetup();
    plcolor(1);     /* set initial color to the forground color */
    
    return 0;
}

/*
 * PLSWAP: does nothing.
 */

int plswap() { return 0; }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

real plxscale(real x, real y)
{
    return x;
}

real plyscale(real x, real y)
{
    return y;
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

int plltype(int lwid, int lpat)
{
    int lw, ls;

    if (iterm==0) return 0;       /* no graphics output requested */
    dprintf(1,"ltype: %d %d\n", lwid,lpat);
    if (lwid > 0)
      giza_set_line_width(1.0*lwid);      
    if (lpat > 0)
      giza_set_line_style(1.0*lpat);
      
    return 0;
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

int plline(real x, real y)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    dprintf(1,"line: %g %g\n", x, y);
    giza_draw(x, y);
    return 0;
}

int plmove(real x, real y)
{
   if (iterm==0) return 0;       /* no graphics output requested */

   dprintf(1,"move: %g %g\n", x, y);   
   giza_move(x,y);
   return 0;
}

int plpoint(real x, real y)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    dprintf(1,"point: %g %g\n", x, y);    
    giza_single_point(x, y, 1);
    return 0;
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

int plcircle(real x, real y, real r)
{
    if (iterm==0) return 0;       /* no graphics output requested */

#if 1
    dprintf(1,"circle: %g %g %g\n", x,y,r);
    giza_circle(x,y,r);
#else
    int npnts, i;
    real theta;

    // this also seems broken in giza, but should work
    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
        theta = TWO_PI * ((real) i) / ((real) npnts);
        plline(x + r * cos(theta), y + r * sin(theta));
    }
#endif
    return 0;
}

int plcross(real x, real y, real s)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    if (s > 0.0) {
        plmove(x - s, y);
        plline(x + s, y);
        plmove(x, y - s);
        plline(x, y + s);
    } else {
        s = s / 1.4142;
        plmove(x - s, y - s);
        plline(x + s, y + s);
        plmove(x - s, y + s);
        plline(x + s, y - s);
    }
    return 0;
}

int plbox(real x, real y, real s)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    if (s > 0.0) {
        plmove(x - s, y - s);
        plline(x + s, y - s);
        plline(x + s, y + s);
        plline(x - s, y + s);
        plline(x - s, y - s);
    } else {
        s = s * 1.4142;
        plmove(x - s, y);
        plline(x, y - s);
        plline(x + s, y);
        plline(x, y + s);
        plline(x - s, y);
    }
    return 0;
}

/*
 * PLJUST: specify justification of strings and numbers.
 * Imports: jus: 
 */

static float fjust = 0.0;   /* pgplot default: left justified */

int pljust(int jus)       /* -1, 0, 1 for left, mid, right just */
{
    if (iterm==0) return 0;       /* no graphics output requested */

    fjust = (jus < -1 ? 0.0 : (jus > 1 ? 1.0 : 0.5));
    return fjust;
}

/*
 * PLTEXT: plot a text string.
 */

int pltext(string msg, real x, real y, real hgt, real ang)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    real just = (fjust+1)/2;  //  -1 .. 1  =>  0 .. 1
    dprintf(1,"text: %g %g %g %g %s\n", x,y,ang,just,msg);
    giza_ptext(x, y, ang, just, msg);
    return 0;
}

/*
 * PLFLUSH: output any pending graphics.
 */

int plflush() 
{ 
    if (iterm==0) return 0;

    giza_flush_device();
    return 0;
}

/*
 * PLFRAME: advance to next frame.
 */

int plframe()
{
    if (iterm==0) return 0;       /* no graphics output requested */

    giza_change_page();
    return 0;
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

int plstop()
{
    if (iterm==0) return 0;       /* no graphics output requested */

    giza_close_device();
    return 0;
}

int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin,
	  real cell,real fmin,real fmax,real findex,real blankval)

{
    int ix,iy,ix0,ix1,iy0,iy1;
    float gray0, gray1, tr[6], *data, *dp;
    real dval,dfac;
    
    if (iterm==0) return 0;       /* no graphics output requested */
    dprintf(0,"matrix: %d %d (gray not working)\n", nx,ny);

    // double affine[6] = {0, 0, 1, 0, 0, 1};
    // double affine[6] = {1, 0, 0, 1, 0, 0};
    double affine[4] = {1, 0, 0, 1};
    giza_render_gray(nx, ny, frame, 0, nx-1, 0, ny-1, fmin, fmax, 0, affine);
     
#if 0
    dprintf(1,"PL_MATRIX development version for YAPP_PGPLOT\n");
    if (findex < 0.0) {
        dprintf(0,"Swapping min and max to : %g %g\n",fmax,fmin);
        findex = -findex;
        dval = fmin;
        fmin = fmax;
        fmax = dval;
    }
    if (fmin==fmax) {
        warning("Cannot greyscale a uniform image yet");
        return 1;
    } else {
        dfac = 1.0/(fmax-fmin);
    }
    ix0 = 1;
    iy0 = 1;
    ix1 = nx;
    iy1 = ny;
    gray0 = 0.0;
    gray1 = 1.0;
    tr[2] = tr[4] = 0.0;    /* off axis matrix scaling */
    tr[1] = 16.0/nx;        /* linear scaling to make it fit in 16x16cm */
    tr[5] = 16.0/ny;
    tr[0] = 2.0 - 8.0/nx;   /* offset to make pgplot fit the ll corner */
    tr[3] = 2.0 - 8.0/ny;
    data = (float *) allocate(sizeof(float)*nx*ny);
    for(iy=0, dp=data; iy<ny; iy++)
        for(ix=0; ix<nx; ix++) {
            dval =  *(frame + ix*ny + iy);
	    if (dval == blankval) {
	      *dp++ = pow(fmin,findex) - 1;
	    } else {
	      if (fmin < fmax) {
                if (dval < fmin) dval=fmin;
                if (dval > fmax) dval=fmax;
                dval = (dval-fmin)*dfac;
	      } else {
                if (dval < fmax) dval=fmax;
                if (dval > fmin) dval=fmin;
                dval = (dval-fmax)*dfac + 1.0;
	      }
            *dp++ = pow(dval,findex);
	    }
        }
    //pggray_(data,&nx,&ny,&ix0,&ix1,&iy0,&iy1,&gray1,&gray0,tr);
    free(data);
#endif    
    return 0;
}

/* not functional yet */

int pl_contour(real *frame,int nx,int ny, int nc, real *c)
{
    int ix,iy,ix0,ix1,iy0,iy1;
    float tr[6], *data, *dp, *cnt;
    real dval,dfac;
    
    if (iterm==0) return 0;       /* no graphics output requested */

    dprintf(0,"PL_CONTOUR development version for YAPP_PGPLOT\n");
    ix0 = 1;
    iy0 = 1;
    ix1 = nx;
    iy1 = ny;
    tr[2] = tr[4] = 0.0;    
    tr[1] = 16.0/nx;
    tr[5] = 16.0/ny;
    tr[0] = 2.0 - 8.0/nx;
    tr[3] = 2.0 - 8.0/ny;
    data = (float *) allocate(sizeof(float)*nx*ny);
    for(iy=0, dp=data; iy<ny; iy++)
        for(ix=0; ix<nx; ix++)
            *dp++ = *(frame + ix*ny + iy);
    cnt = (float *) allocate(sizeof(float)*nc);
    for(ix=0; ix<nc; ix++)
        cnt[ix] = c[ix];
    //pgcont_(data,&nx,&ny,&ix0,&ix1,&iy0,&iy1,cnt,&nc,tr);
    free(data);
    free(cnt);
    return 0;
}


int pl_screendump(string fname)
{
  printf("pl_screendump(%s): Not implemented for yapp_giza\n",fname);
  return 0;
}

int pl_getpoly(float *x, float *y, int n)
{
    int nn, delay, loc, k, symbol;
    float xold,yold, xnew,ynew;
    char ch[10];

    printf("Define a polygon:\n");
    printf("   LEFT   = define first/next vertex of polygon\n");
    printf("   MIDDLE = delete previous vertex point from polygon\n");
    printf("	    (this option cannot remove already plotted vertex lines\n");
    printf("   RIGHT  = close polygon\n");
#if 0
    nn = 0;                           /* count points in polygon */
    for(;;) {
        k = pgcurse_(&xnew, &ynew, ch, 1);
        if(k==0) {
            warning("Device has no cursor...");
            return 0;
        }
        if (ch[0]=='A') {                     /* left button = DEFINE POLYGON */
            xold = xnew;
            yold = ynew;
            dprintf (2,"Button A at x=%f y=%f\n",xold,yold);
            if (nn==0) {
                pgmove_(&xold,&yold);
                k=1;
                symbol=2;
                pgpoint_(&k,&xold,&yold,&symbol);
#if 0
                set_marker_symbol(43);      /* a 'plus' symbol */
                marker_abs_2(xold,yold);
        
#endif
            } else
                pgdraw_(&xold,&yold);
            if (nn<n) {
               x[nn] = xold;
               y[nn] = yold;
               nn++;
            } else
                warning("No more space to store this data");
        } else if (ch[0]=='D') {                  /* middle  = DELETE PREVIOUS POINT */
            if (nn>0)
                nn--;
            if (nn!=0)
               pgmove_(&x[nn-1],&y[nn-1]);    /* reset */
        } else if (ch[0]=='X') {                          /* right = QUIT */
            break;
        } else
            warning("Unknown return char %c from PGCURSE (expected A,D,X)",ch[0]);
    }   /* for(;;) */
    dprintf (2,"Button C pressed to finish up polygon (%d)\n",nn);
    if (nn>0) {
        pgdraw_(&x[0],&y[0]);       /* close polygon */
    }
    return nn<3 ? 0 : nn;
#else
    return 0;
#endif
}

int pl_cursor(real *x, real *y, char *c)
{
    char inf[8], ans[8];
    int len, inf_len, ans_len;
    permanent float xsave, ysave;
#if 0                
    strcpy(inf,"CURSOR");
    inf_len = strlen(inf);
    ans_len = 1;
    pgqinf_(inf, ans, &len, inf_len, ans_len);
    if (ans[0]=='n' || ans[0]=='N') return 0;
    pgcurs_(&xsave, &ysave, c, 1);
    *x = xsave;
    *y = ysave;
#endif    
    return 1;
}



/*
 * CMS_RGBSETUP: default color table initialization.
 */

#define	BLACK		0
#define RED		1
#define GREEN		2
#define BLUE		3
#define CYAN 		4	/* (Green + Blue) */	
#define MAGENTA 	5	/* (Red + Blue)   */	
#define YELLOW		6	/* (Red + Green)  */	
#define ORANGE		7
#define GREEN_YELLOW	8
#define GREEN_CYAN	9
#define BLUE_CYAN	10
#define BLUE_MAGENTA	11
#define RED_MAGENTA	12
#define DARK_GRAY	13
#define LIGHT_GRAY	14
#define	WHITE   	15

local void cms_rgbsetup()
{
  dprintf(1,"cms_rgbsetup\n");
  ncolors = WHITE+1;
  /* default PGPLOT  colors: although, defined here, plpalette is not called */
  red[BLACK]=0.00;         green[BLACK]=0.00;          blue[BLACK]=0.00;
  red[RED]=1.00;           green[RED]=0.00;            blue[RED]=0.00;
  red[GREEN]=0.00;         green[GREEN]=1.00;          blue[GREEN]=0.00;
  red[BLUE]=0.00;          green[BLUE]=0.00;           blue[BLUE]=1.00;
  red[CYAN]=0.00;          green[CYAN]=1.00;           blue[CYAN]=1.00;
  red[MAGENTA]=1.00;       green[MAGENTA]=0.00;        blue[MAGENTA]=1.00;
  red[YELLOW]=1.00;        green[YELLOW]=1.00;         blue[YELLOW]=0.00;
  red[ORANGE]=1.00;        green[ORANGE]=0.50;         blue[ORANGE]=0.00;
  red[GREEN_YELLOW]=0.50;  green[GREEN_YELLOW]=1.00;   blue[GREEN_YELLOW]=0.00;
  red[GREEN_CYAN]=0.00;    green[GREEN_CYAN]=1.00;     blue[GREEN_CYAN]=0.50;
  red[BLUE_CYAN]=0.00;     green[BLUE_CYAN]=0.50;      blue[BLUE_CYAN]=1.00;
  red[BLUE_MAGENTA]=0.50;  green[BLUE_MAGENTA]=0.00;   blue[BLUE_MAGENTA]=1.00;
  red[RED_MAGENTA]=1.00;   green[RED_MAGENTA]=0.00;    blue[RED_MAGENTA]=0.50;
  red[DARK_GRAY] =0.33;    green[DARK_GRAY]=0.33;      blue[DARK_GRAY]=0.33;
  red[LIGHT_GRAY]=0.66;    green[LIGHT_GRAY]=0.66;     blue[LIGHT_GRAY]=0.66;
  red[WHITE]=1.00;         green[WHITE]=1.00;          blue[WHITE]=1.00;
  // WHITE is the last color
}

/*
 * PLCOLOR: specify new plotting color as an integer between 0 and ncolors-1;
 * normally values outside this range are mapped to the nearest endpoint.
 * For PGPLOT however, color 0 is the background, 1 the foreground
 */

void plcolor(int color)
{
    if (color < 0)
	color = 0;
    else if (color > ncolors - 1)
        color = 1;
    //pgsci_(&color);
    dprintf(1,"color: %d/%d\n",color,ncolors);
    giza_set_colour_index(color);
}

/*
 * PLNCOLORS: return current value of local variable ncolors.
 */

int plncolors()
{
    return ncolors;
}

/*
 * PLPALETTE: re-initialize color table from user-supplied values.
 */

void plpalette(real *r, real *g, real *b, int nc)
{
    int i;
#if 0
    if (nc > MAXCOLOR)
	error("plpalette: cannot define more than %d colors", MAXCOLOR);
    ncolors = nc;
    for (i = 0; i < ncolors; i++) {
	red[i] = r[i];
	green[i] = g[i];
	blue[i] = b[i];
	dprintf(1,"->PGSCR_(%d,%g,%g,%g) \n",i,red[i],green[i],blue[i]);
        pgscr_(&i,&red[i],&green[i],&blue[i]);
    }
    red[ncolors] = green[ncolors] = blue[ncolors] = 0.0;    /* terminate */
    plcolor(1);			/* reset default color to foreground */
    dprintf(0,"Setting default color to: %g %g %g\n",
            red[1], green[1], blue[1]);
#endif    
}

/*
 * PLLUT:  read a new RGB palette from a LUT color table
 *	   The lut must be a simple ascii file, with 1st, 2nd and 3rd column being
 *	   the R, G and B response, a real number between 0.0 and 1.0
 */
#define NOCOLOR(x)  ((x)<0||(x)>1)
void pllut(string fname, bool compress)
{
    stream cstr;
    int ncolors, nskip=0;
    real red[MAXCOLOR], green[MAXCOLOR], blue[MAXCOLOR];
    float r, g, b, r_old, g_old, b_old;
    char line[256];

    if (fname==NULL || *fname == 0) return;
    cstr = stropen(fname, "r");
    ncolors = 0;
    while (fgets(line,256,cstr)) {
       if(ncolors>=MAXCOLOR) error("(%d/%d): Too many colors in %s",
              ncolors+1, MAXCOLOR, fname);
       sscanf(line,"%f %f %f",&r, &g, &b);
       if (NOCOLOR(r) || NOCOLOR(g) || NOCOLOR(b)) {
          warning("Skipping RGB=%f %f %f",r,g,b);
	  continue;
       }
       if (ncolors>0 && r==r_old && g==g_old && b==b_old && compress) {
          nskip++;
       } else {
          dprintf(1,"LUT(%d): %f %f %f\n",ncolors,r,g,b);
          red[ncolors]   = r;   r_old = r;
          green[ncolors] = g;   g_old = g;
          blue[ncolors]  = b;   b_old = b;
          ncolors++;
       }
    }
    strclose(cstr);
    dprintf(0,"Colortable %s: %d entries; %d redundant entries skipped\n",
            fname, ncolors,nskip);
    plpalette(red, green, blue, ncolors);
}

