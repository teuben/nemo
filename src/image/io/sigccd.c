/*
 *   SIGCCD:   turn a table into a ccd file - for WASP only
 *
 *	12-nov-98	coded, quick and dirty as usual (from fitsccd)
 *         dec-00       added nx= (aphid has nx=16, wasp has nx=128)
 *       1-jan-04       get_line changed interface, fixed to nemo_file_lines
 */

#include <stdinc.h>
#include <getparam.h>
#include <image.h>


string defv[] = {
    "in=???\n		Input sig (table) file",
    "out=???\n		Output image file",
    "nx=16\n            first dimension???",
    "dcol=2\n           Column for data",
    "nmax=100000\n      Allocation space for piped I/O",
    "VERSION=1.1b\n	12-apr-04 PJT",
    NULL,
};

string usage = "convert (WASP sig)tables into images";

#if !defined(MAX_LINELEN)
#define MAX_LINELEN 10000
#endif

void nemo_main()
{
    stream outstr;
    char line[MAX_LINELEN];
    stream instr;
    double freq, freqold, data[3];
    int nmax, nx, ny, i, j, ierr;
    imageptr iptr;
    real rmin, rmax;
    string name = getparam("in");
    int dcol = getiparam("dcol") - 1;
    char namex[32], namey[32];

    
    nmax = nemo_file_lines(name, getiparam("nmax"));
    nx = getiparam("nx");
    ny = nmax / (nx+1);
    if (nmax % (nx+1))
        warning("odd sized file: nlines=%d nx=%d ny=%d?\n",nmax,nx,ny);
    dprintf(0,"%s: Found %d x %d image\n",name,nx,ny);
    create_image(&iptr,nx,ny);
    outstr = stropen(getparam("out"),"w");
    strcpy(namex,"Lag");
    strcpy(namey,"Frequency");

    instr = stropen(name,"r");
    for (j=0; j<ny; j++) {
        if (get_line(instr, line) <= 0)
	  error("blank line,or no more data found for row j=%d\n",j+1);
        freq = atof(line);
        dprintf(1,"%d : freq=%g\n",j,freq);
        if (j==1) {    /* at the 2nd point, figure out the linear delta */
            Xmin(iptr) = 0.0;
            Ymin(iptr) = freqold;
            Dx(iptr) = 1.0;
            Dy(iptr) = freq - freqold;
        }
        freqold = freq;
        for (i=0; i<nx; i++) {
	    if (get_line(instr, line) <= 0)
	      error("blank line, or no more data found for col i=%d\n",i+1);
            patch_line(line);
            if ((ierr=nemoinpd(line,data,3)) != 3)
                error("(%d,%d) badly formatted line: %s; err=%d\n",i,j,line,ierr);
            MapValue(iptr,i,j) = data[dcol];
            if (i==0 && j==0) {
                rmin = rmax = data[dcol];
            } else {
                rmin = MIN(rmin, data[dcol]);
                rmax = MAX(rmax, data[dcol]);
            }
        }
    }
    MapMin(iptr) = rmin;
    MapMax(iptr) = rmax;
    Namex(iptr) = namex;
    Namey(iptr) = namey;
    
    write_image(outstr,iptr);
    strclose(outstr);
}

/* darn gipsy doesn't like tab's */

patch_line(char *line)
{
    while (*line) {
        if (*line == '\t') *line = ' ';
        line++;
    }
}
