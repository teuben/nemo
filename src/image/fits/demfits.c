/*
 * DEM: read a 7.5' DEM (Digital Elevation Maps) from USGS
 *      and frame it into something we want
 *
 *      11-dec-96   V0.1    trial run for HatCreek new sites        pjt
 *                          since no real fits file written yet, you
 *                          need a shell script with this to go.
 *      12-dec-96   V0.2    attempt to check more, and merge back 1deg maps
 */

#include <stdinc.h>
#include <getparam.h>
#include <strlib.h>

string defv[] = {
    "in=???\n           Input DEM file",
    "out=???\n          Output FITS file",
    "naxis=0,0\n        Pre-requested size in X,Y",
    "crval=0,0\n        Pre-requested coordinates of lower left corner",
    "newfmt=f\n		Use new format (has data elements 17..28 in A)",
    "c=f\n		Read C?",
    "tab=\n             Output table if IX,IY,IZ\n",
    "VERSION=0.2\n      12-dec-96 PJT",
    NULL,
};

string usage = "Read DEM (Digital Elevation Maps) from USGS";

local void   read_dem(string, string, string, bool, int *, double *);

local string get_string(stream fp,  int len);
local byte   get_byte(stream fp, int len);
local short  get_short(stream fp, int len);
local float  get_float(stream fp, int len);
local double get_double(stream fp, int len);


nemo_main()
{
    string infile, outfile, tabfile;
    bool Qnewfmt = getbparam("newfmt");
    int naxis[2];
    double crval[2];

    infile = getparam("in");
    outfile = getparam("out");
    if (hasvalue("tab"))
        tabfile = getparam("tab");
    else
        tabfile = 0;
    (void) nemoinpi(getparam("naxis"),naxis,2);
    (void) nemoinpd(getparam("crval"),crval,2);

    read_dem(infile, outfile, tabfile, Qnewfmt, naxis, crval);
}



#define UNIT_RAD    0
#define UNIT_FEET   1
#define UNIT_METER  2
#define UNIT_ARCSEC 3

struct record_a {
    char   file_name[40+1], 
           free_text[40+1],
           filler1[55+1],
           process_code[1+1],
           filler2[1+1],
           sec_ind[3+1];
    char   MC_origin_code[4+1];
    short  DEM_level_code;
    short  elev_pattern;
    short  planim_ref;
    short  planim_zone;
    double projpar[15];
    short  cunits;
    short  eunits;
    short  npoly;
    double corners[4][2];
    double minmax[2];
    double angle;
    short  accuracy;
    float  step[3];
    short  size[2];     /* m rows, n columns ;  here m==1, n > 1 */

    /* old format stopped here */

    short  lpcint;
    byte   scintu1;
    short  spcint;
    byte   scintu2;
    short  sdate;
    short  rdate;
    char   rflag[1];
    byte   dvflag;

    byte   vaflag;
    byte   vdatum;
    byte   hdatum;
    short  dedition;
} A;

struct record_b {
    short row;
    short col;
    short m;            /* 'row' check */
    short n;            /* 'col' check */
    double xgo;
    double ygo;
    double elev;
    double emin;
    double emax;
    short *data;        /* NOT USED len: m */
} B;

struct record_c {
    short has_stat2;
    short rmse2[3];
    short sam_size2;     /* INTEGER*3 : check if this isn't a short ??? */
    short has_stat5;
    short rmse5[3];
    short sam_size5;
} C;

#define ROUNDUP(a,b) ((b)*(((a)+(b)-1)/(b)))

void read_dem(string fname, string oname, string tname,
              bool Qnewfmt, 
              int *naxis, double *crval)
{
    stream instr, outstr, tabstr;
    int i, j, n, count, offset, col, naxis1;
    double xmin, xmax, ymin, ymax;
    short sdata[2048];

    instr = stropen(fname,"r");
    outstr = stropen(oname,"w");
    if (tname)
        tabstr = stropen(tname,"w");
    else
        tabstr = NULL;

    strcpy(A.file_name,      get_string(instr,40));
    strcpy(A.free_text,      get_string(instr,40));
    strcpy(A.filler1,        get_string(instr,55));
    strcpy(A.process_code,   get_string(instr,1));
    strcpy(A.filler2,        get_string(instr,1));
    strcpy(A.sec_ind,        get_string(instr,3));

    printf("File name:      %s\n",A.file_name);
    printf("other info:     %s\n",A.free_text);
    printf("Process code:   %s\n",A.process_code);


    strcpy(A.MC_origin_code, get_string(instr,4));
    
    A.DEM_level_code = get_short(instr,6);
    A.elev_pattern = get_short(instr,6);
    A.planim_ref = get_short(instr,6);
    A.planim_zone = get_short(instr,6);
    for (i=0; i<15; i++)
        A.projpar[i] = get_double(instr,24);
    A.cunits = get_short(instr,6);
    A.eunits = get_short(instr,6);
    A.npoly = get_short(instr,6);
    for (i=0; i<4; i++)
    for (j=0; j<2; j++)
        A.corners[i][j] = get_double(instr,24);
    A.minmax[0] = get_double(instr,24);
    A.minmax[1] = get_double(instr,24);
    A.angle = get_double(instr,24);
    A.accuracy = get_short(instr,6);
    for (i=0; i<3; i++)
        A.step[i] = get_float(instr,12);
    A.size[0] = get_short(instr,6);
    A.size[1] = get_short(instr,6);

    if (A.size[0] != 1) error("Cannot handle A.size[0] = %d",A.size[0]);
    if (A.angle != 0.0) warning("Cannot handle rotated (%g) maps",A.angle);

    printf("DATAMIN: %g\n",A.minmax[0]);
    printf("DATAMAX: %g\n",A.minmax[1]);
    printf("CDELT1: %g\n",A.step[0]/3600.0);
    printf("CDELT2: %g\n",A.step[1]/3600.0);
    printf("CRPIX1: 1.0\n");
    printf("CRPIX2: 1.0\n");
    printf("CRVAL1: %g\n",A.corners[0][0]/3600.0);
    printf("CRVAL2: %g\n",A.corners[0][1]/3600.0);

    /* old format (data-items 1-16) stopped here */

    if (Qnewfmt) {                          /* new items (17-28) */
        A.lpcint  = get_short(instr,5);
        A.scintu1 = get_byte(instr,1);
        A.spcint  = get_short(instr,5);
        A.scintu2 = get_byte(instr,1);

        A.sdate   = get_short(instr,4);
        A.rdate   = get_short(instr,4);
        strcpy(A.rflag, get_string(instr,1));
        A.dvflag  = get_byte(instr,1);
        A.vaflag   = get_byte(instr,1);
        A.vdatum   = get_byte(instr,1);
        A.hdatum   = get_byte(instr,1);
        A.dedition = get_short(instr,4);

        dprintf(1,"sdate/rdate=%d %d dvflag=%d\n",
        	A.sdate, A.rdate, A.dvflag);
    }
    fseek(instr, ROUNDUP(ftell(instr),1024), SEEK_SET);
    dprintf(2,"0x%x: finished Reading A\n",ftell(instr));

    for (col=0; col<A.size[1]; col++) {


        B.row = get_short(instr,6);
        B.col = get_short(instr,6);
        B.n   = get_short(instr,6);
        B.m   = get_short(instr,6);
        B.xgo = get_double(instr,24);
        B.ygo = get_double(instr,24);
        B.elev = get_double(instr,24);
        B.emin  = get_double(instr,24);
        B.emax  = get_double(instr,24);

        if (naxis[0] == 0 || naxis[1] == 0) {
            if (B.n != A.size[1]) 
                error("Cannot handle non-square %d * %d maps", B.n, A.size[1]);
            warning("Setting mapsize %d * %d",B.n,A.size[1]);
            naxis[0] = A.size[1];
            naxis[1] = B.n;
        }

        for (i=0; i<naxis[1]; i++) sdata[i] = 0;
    	if (B.m != 1) error("Cannot handle B.m = %d",B.m);

        if (col==0) {
            xmin = xmax = B.xgo;
            ymin = B.ygo;
            ymax = ymin + (B.n - 1)*A.step[1];
            if (crval[0] != 0.0) {
	    	count = (xmin-crval[0])/A.step[0];
            	dprintf(0,"Writing %d blank columns first\n",count);
            	for (i=0; i<count; i++)
                    (void) fwrite(sdata,sizeof(short),naxis[1],outstr);
                naxis1 = count;
            } else {
                count = 0;
                naxis1 = 0;
            }
        } else {
            xmin = MIN(xmin, B.xgo);
            ymin = MIN(ymin, B.ygo);
            xmax = MAX(xmax, B.xgo);
            ymax = MAX(ymax, B.ygo + (B.n - 1)*A.step[1]);
        }

        count = 0;
        if (crval[1] != 0)
            offset = (B.ygo - crval[1])/A.step[1];
        else
            offset = 0;
        n = B.n;


        dprintf(1,"   %d %d %d %d Xgo=%g Ygo=%g offset=%d\n",
			B.row, B.col, B.n, B.m, B.xgo,B.ygo,offset);


        for (i=0; i<146 && count<n; i++, count++) {
            sdata[offset+count] = get_short(instr,6);
            if(tabstr) fprintf(tabstr,"%d %d %d\n",
                    col,offset+count,sdata[offset+count]);
        }
        fseek(instr, ROUNDUP(ftell(instr),1024), SEEK_SET);
        dprintf(2,"0x%x: finished Reading B1 (%d,%d,%d,%d) elev=%g [%g - %g]\n",
                ftell(instr),B.row,B.col,B.n,B.m,B.elev,B.emin,B.emax);

        while (count < n) {
            for (i=0; i<170 && count<n; i++, count++) {
                sdata[offset+count] = get_short(instr,6);
                if(tabstr) fprintf(tabstr,"%d %d %d\n",
                        col,offset+count,sdata[offset+count]);
            }
            fseek(instr, ROUNDUP(ftell(instr),1024), SEEK_SET);
            dprintf(3,"0x%x: finished Reading another B...\n",ftell(instr));
        }

        dprintf(2,"0x%x: finished Reading B; count=%d\n",ftell(instr),count);
        fwrite(sdata,sizeof(short),naxis[1],outstr);
        naxis1++;
    }
    
    count = naxis[0] - naxis1;
    dprintf(0,"Writing %d blank columns last\n",count);
    for (i=0; i<naxis[1]; i++) sdata[i] = 0;
    for (i=0; i<count; i++)
        (void) fwrite(sdata,sizeof(short),naxis[1],outstr);

    printf("XMIN: %g\n",xmin);
    printf("XMAX: %g\n",xmax);
    printf("YMIN: %g\n",ymin);
    printf("YMAX: %g\n",ymax);
    printf("NAXIS1: %g\n",1+(xmax-xmin)/A.step[0]);
    printf("NAXIS2: %g\n",1+(ymax-ymin)/A.step[1]);

}

local char buffer[256];

local string get_string(stream fp,  int len)
{
    int pos = ftell(fp);
    int n = fread(buffer,1,len,fp);
    if (n != len) error("string: read %d, expected %d",n,len);
    buffer[len] = 0;
    dprintf(4,"0x%x: string=%s\n",pos,buffer);
    return buffer;
}

local byte get_byte(stream fp, int len)
{
    int pos = ftell(fp);
    int n = fread(buffer,1,len,fp);
    if (n != len) error("string: read %d, expected %d",n,len);
    buffer[len] = 0;
    dprintf(4,"0x%x: byte=%s\n",pos,buffer);
    return (byte) atoi(buffer);
}

local short get_short(stream fp, int len)
{
    int pos = ftell(fp);
    int n = fread(buffer,1,len,fp);
    if (n != len) error("string: read %d, expected %d",n,len);
    buffer[len] = 0;
    dprintf(4,"0x%x: short=%s\n",pos,buffer);
    return (short) atoi(buffer);
}

local float get_float(stream fp, int len)
{
    float fval;
    int pos = ftell(fp);
    int n = fread(buffer,1,len,fp);
    if (n != len) error("string: read %d, expected %d",n,len);
    buffer[len] = 0;
    dprintf(4,"0x%x: float=%s\n",pos,buffer);
#if 0
    return (float) atof(buffer);
#else
    nemoinpf(buffer,&fval,1);
    return fval;
#endif
}

local double get_double(stream fp, int len)
{
    double dval;
    int pos = ftell(fp);
    int n = fread(buffer,1,len,fp);
    if (n != len) error("string: read %d, expected %d",n,len);
    buffer[len] = 0;
    dprintf(4,"0x%x: double=%s\n",pos,buffer);
#if 0
    return (double) atof(buffer);
#else
    nemoinpd(buffer,&dval,1);
    return dval;
#endif
}


