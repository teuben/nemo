/*
 *  MILLER - read Miller ``amazing'' plot files
 *
 *	11-jul-94   V1.0 	allright, I gave it a try		pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n      Input file name - Dick Millers plot file",
    "out=\n        Output file name - snapshot",
    "headline=\n   Random mumblage for humans",
    "VERSION=1.0\n 12-jul-94 PJT",
    NULL,
};

string usage = "convert Miller (ascii) plot file to snapshot format";


#define MAXLINE  256
#define MAXTIMES 256
#define MAXVALS   32

stream instr, outstr=NULL;
int snapcnt=0;

int nbody=0, nbody_max=0, bits=0;
real tsnap=0.0;
Body *btab = NULL;

real times[MAXTIMES];
int ntimes=0;

/* (global) Miller variables : */
int  step, part, pactiv, vspill, cspill, kpart, kact, lim, nmax;
real crtim, rconf, ke, pe, vmax, asq, won;

/* extern's */
extern string *burststring ARGS((string, string));

/* local's */
local bool in_header();
local int get_double(int, double*);
local int get_int(int, int*);
local int get_id(char*);
     

void nemo_main(void)
{
    int i;
    char name[20];

    instr = stropen(getparam("in"), "r");
    if (hasvalue("out"))
        outstr = stropen(getparam("out"), "w");
    if (hasvalue("headline")) set_headline(getparam("headline"));
    if (outstr) put_history(outstr);
    bits = TimeBit | PhaseSpaceBit | KeyBit;
    
    while (in_header()) {	/* loop reading data frames */
        snapcnt++;
        if (outstr) put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
    
    if (outstr) strclose(outstr);
}

local bool in_header()
{
    double dval, dvals[32];
    int i, ival, ivals[32];
    char run_id[80];
    Body *bp;

    i = get_id(run_id);		/* one id string */
    if (i<1) return FALSE;
    set_headline(run_id);
    get_int(30,ivals);		/* 30 integers, over 2 lines normally */
    get_double(32,dvals);	/* 32 floats, over 4 lines normally */

    step = ivals[0];
    part = ivals[1];
    pactiv = ivals[2];
    vspill = ivals[3];
    cspill = ivals[4];
    kpart = ivals[5];
    kact = ivals[6];
    lim = ivals[7];
    nmax = ivals[13];

    crtim = dvals[55-32];
    rconf = dvals[56-32];
    ke = dvals[57-32];
    pe = dvals[58-32];
    vmax = dvals[61-32];
    asq = dvals[62-32];
    won = dvals[63-32];
    
    dprintf(1,"step=%d part=%d pactiv=%d vspill=%d cspill=%d\n",
                step, part, pactiv, vspill, cspill);
    dprintf(1,"kpart=%d kact=%d lim=%d nmax=%d\n",
                kpart, kact, lim, nmax);
    dprintf(1,"ctim=%g rconf=%g ke=%g pe=%g vmax=%g asq=%g won=%g\n",
                crtim, rconf, ke, pe, vmax, asq, won);

    dprintf(0,"[Read %d/%d bodies (sampled from %d) at step %d]\n",
                kact,kpart, part, step);

    /* set our local variables for the snapshot to be output */
    
    nbody =  kact;
    tsnap = step;

    if (btab==NULL) {
            btab = (Body *) allocate(nbody*sizeof(Body));
            nbody_max = nbody;
    }
    if (nbody > nbody_max) {
            warning("Reallocating snapshot from %d to %d",nbody,nbody_max);
            btab = (Body *) reallocate(btab,nbody*sizeof(Body));
            nbody_max = nbody;
    }

    for (i=0, bp=btab; i<kpart; i++) {
        if (get_int(4,ivals) != 4) return FALSE;
        if (ivals[0] == 65535 &&
            ivals[1] == 65535 &&
            ivals[2] == 65535 &&
            ivals[3] == 65535) continue;
        Pos(bp)[0] = ivals[0]/1024.0;
        Pos(bp)[1] = ivals[1]/1024.0;
        Pos(bp)[2] = ivals[2]/1024.0;
        Vel(bp)[0] = 0.0;
        Vel(bp)[1] = 0.0;
        Vel(bp)[2] = 0.0;
        Key(bp) = ivals[3];
        bp++;
    }

    return TRUE;
}

/*
 * GET_DOUBLE:  get 'n' double values from a header, not needed 
 * 	        to be on one line. 
 *
 * NOTE: not sure if the line on which last number (n) was
 *       read is still accessible, since I use fgets() in
 *       reading non-header information 
 *	 The extra \n seems to have done the trick.... not sure
 *	 about portability.... ought to use cashed fgets()+nemoinp()
 */

local int get_double(int n, double *d)
{
    int i, k=0;
    char line[MAXLINE];
    while (k<n) {
        if (fgets(line, MAXLINE, instr) == NULL) {
            if (k==0) return 0;
            warning("Unexpected EOF in header for snapshot %d",snapcnt+1);
            return k;
        }
        if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
        i = nemoinpd(line,&d[k],n-k);
        if (i==0) return k;
        if (i<0) error("Error parsing %s",line);
        k += i;
    }
    return k;
}

local int get_int(int n, int *d)
{
    int i, k=0;
    char line[MAXLINE];
    while (k<n) {
        if (fgets(line, MAXLINE, instr) == NULL) {
            if (k==0) return 0;
            warning("Unexpected EOF in header for snapshot %d",snapcnt+1);
            return k;
        }
        if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
        i = nemoinpi(line,&d[k],n-k);
        if (i==0) return k;
        if (i<0) error("Error parsing %s",line);
        k += i;
    }
    return k;
}

local int get_id(char *s)
{
  int i, k=0;
  char *cp, line[MAXLINE];


  for(;;) {
    if (fgets(line, MAXLINE, instr) == NULL) {
        if (k==0) return 0;
        warning("Unexpected EOF in header for snapshot %d",snapcnt+1);
        return k;
    }

    if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0; /* trim end */
    cp = line;
    while (*cp == ' ') cp++;	/* skip leading blanks */
    if (*cp) break;
  }
  strcpy(s,cp);
  return strlen(s);
}


