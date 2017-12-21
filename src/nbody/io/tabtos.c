/*
 *  TABTOS - generalized tables to snapshot converter
 *
 *           atos(ph) was too cumbersome, as shown by the following 
 *           sequence to read in pure SPH data:
 *
 *           (in, out, var and color are useful shell variables)
 *
 *               atosph $in - | tsf - allline=t |\
 *               sed s/$var/Aux/g  | rsf - - |\
 *               snapcenter - - tab=f |\
 *               snapxyz - $out "color=$color" 
 *
 *           moreover, tabtos is 2.0 times faster than atosph,
 *		probably due to fscanf()  vs. fgets()+nemoinp()
 *
 *      26-aug-93   V1.0    Created - while in Tokyo           PJT
 *      30-aug-93   V1.1    Added options=scan                 PJT
 *	10-dec-93   V1.2    renamed time= to times=	       pjt
 *	22-feb-94   V1.2a   ansi headers
 *	 2-sep-94   V1.2b   fixed bug in get_double(), didn't return #
 *	26-oct-94   V1.2c   added options=comment 	pjt
 *       2-nov-94   V1.2d   added options=wrap,spill for Sellwood	pjt
 *			    (@Rutgers)
 *	13-dec-94   V1.2f   patch nbody if not enough nbody could be read pjt
 *			    (needs error=1)
 *	15-mar-95   V1.3    allow x,y,z   vx,vy,vz for pos,vel		pjt
 *	 3-aug-95   V1.3a   allow TABs in files				pjt
 *	25-mar-97	b   SINGLEPREC fix				pjt
 *       7-jul-98       c   fixed bug skipping blocks                   pjt
 *      19-aug-00       d   fixed TABs correctly                        pjt
 *	23-sep-01       e   ->nemo_file_lines
 *      24-jan-02       f   to block10, time to use indexed parameter   pjt
 *       3-feb-02       g   fixed bug if time in header, nbody not      pjt 
 *      11-may-02       h   longer lines to deal with DCR's 'ss2bt' files   PJT
 *      28-may-03       i   allow 'm' also for 'mass', also allow - for skip  pjt
 *      21-sep-03       j   fixed bug in calling nemo_file_lines              pjt
 *      29-jul-05   V1.4    added auto-incrementing time option (Amsterdam, cafe amaricain) pjt
 *      14-nov-06   V1.5    add first # comments to the NEMO output history.
 *      18-Jan-12   V1.5a   add 'dens' array                           jcl
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n      Input file name - ascii (structured) table",
    "out=\n        Output file name - snapshot",
    "header=\n     Header scalars (nbody, ndim, time, skip, -)",
    "block1=\n     What's in first block of nbody lines - mass",
    "block2=\n     - pos(ndim)",
    "block3=\n     - vel(ndim)",
    "block4=\n     - phase(2*ndim)",
    "block5=\n     - acc(ndim)",
    "block6=\n     - phi",
    "block7=\n     - aux",
    "block8=\n     - key",
    "block9=\n     - skip",
    "block10=\n    - skip",
    "nbody=\n      Number of bodies, if needed",
    "ndim=\n	   Dimension of pos,vel,acc arrays, if needed",
    "times=\n      Time(s) of snapshot, if needed, or override",
    "options=\n    Other processing options (scan|comment|wrap|spill|time)",
    "nskip=0\n     Number of lines skipped before each (header+block1+...)",
    "headline=\n   Random mumblage for humans",
    "VERSION=1.5c\n 16-nov-2017 PJT",
    NULL,
};

string usage = "convert ASCII tabular information to snapshot format";

string cvsid = "$Id$";


#define MAXLINE  1024
#define MAXTIMES 1024
#define MAXVALS    32
#define MAXHIST   128

local stream instr, outstr;

local int nbody=0, nbody_max=0, bits=0;
local real tsnap=0.0;
local Body *btab = NULL;

local string *header;
local int nheader=0;

#define MAXBLOCK 10
local string block[MAXBLOCK];
local int nblocks=0;

local real times[MAXTIMES];
local int ntimes=0;

local int nskip=0;
local int linecnt=0;
local int snapcnt=0;
local bool scan, wrap, spill, auto_time, Qhis, Qcom;


local int nhist = 0;
local string shist[MAXHIST];

/* extern's */
extern string *burststring(string, string);
extern bool scanopt(string, string);

/* local's */
local void check_options(void);
local bool get_header(void);
local bool get_block(int, string );
local int get_double(int, double *);
local int get_nbody(void);
local void do_scan(stream);     
local void tab2space(string);

void nemo_main(void)
{
    int i;
    char name[20];
    bool ok;

    check_options();

    header = burststring(getparam("header"),", ");
    nheader = xstrlen(header,sizeof(string)) - 1;

    for (i=0; i<MAXBLOCK; i++) {
        sprintf(name,"block%d",i+1);
        block[i] = getparam(name);
        if (hasvalue(name))
	    nblocks++;
        else
            break;
    }
    /* additional error check if not in scan mode */
    if (!scan) {
        if (!hasvalue("out")) error("Need output filename; out=");
        if (nblocks==0) error("Need at least one block1=; maximum is block9=");
    }

    ntimes = nemoinpr(getparam("times"),times,MAXTIMES);
    if (ntimes < 0) error("Parsing times=, perhaps too many? (<%d)",MAXTIMES);

    instr = stropen(getparam("in"), "r");
    if (scan) {
        do_scan(instr);
        strclose(instr);
        stop(0);
    }
    
    outstr = stropen(getparam("out"), "w");
    if (hasvalue("headline")) set_headline(getparam("headline"));

    while (get_header()) {	/* loop reading snapshots */
        snapcnt++;
        if (nbody==0) {         /* if unknown, try and estimate  */
            nbody = get_nbody();
            dprintf(0,"Assuming nbody=%d\n",nbody);
        }
        if (ntimes > 0) {		/* override or set time */
            if (snapcnt <= ntimes)
                tsnap = times[snapcnt-1];
            else
                warning("Too few snapshot times= specified, good luck");
        }
	if (auto_time)
	  tsnap += 1.0;
        dprintf(0, "[reading %d bodies at time %f]\n", nbody, tsnap);
        if (btab==NULL) {
            btab = (Body *) allocate(nbody*sizeof(Body));
            nbody_max = nbody;
	}
        if (nbody > nbody_max) {
            warning("Reallocating snapshot from %d to %d",nbody,nbody_max);
            btab = (Body *) reallocate(btab,nbody*sizeof(Body));
            nbody_max = nbody;
        }
        if (nbody <= 0) {
            warning("Bad nbody, possibly wrong number of block's specified");
            break;
        }

        bits = TimeBit;			/* reset the bits */
        for (i=0; i<nblocks; i++) {	/* process for each block */
            ok = get_block(i+1,block[i]);
            if (!ok) {
		break;
	    }
        }
        if (ok) {
	  if (nhist > 0) {
	    for (i=0; i<nhist; i++)
	      app_history(shist[i]);
	    nhist = 0;
	  }
	  put_history(outstr);
	  put_snap(outstr, &btab, &nbody, &tsnap, &bits);	      /* output */
        } else
	  break;
        if (ntimes > 0 && snapcnt == ntimes) break;
    }
    strclose(outstr);
} /* nemo_main */

void check_options(void)
{
    string options;
    bool scanopt();

    // if (!hasvalue("options")) return;

    options=getparam("options");
    scan = scanopt(options,"scan");
    Qcom = scanopt(options,"comment");
    wrap = scanopt(options,"wrap");
    spill = scanopt(options,"spill");
    auto_time = scanopt(options,"time");
    Qhis = scanopt(options,"history");
    if (Qhis) Qcom = TRUE;

    nskip = getiparam("nskip");
}

local bool get_header(void)
{
    int i, ndim=0;
    double dval, *dvals;
    char line[MAXLINE];

    if (feof(instr)) return FALSE;

    dprintf(0,"get_header(nskip=%d)\n",nskip);

    /* @TODO:  nskip appears not to work??? */
    for (i=0; i<nskip; i++) {
        if (fgets(line,MAXLINE,instr) == NULL) {
            if (i==0) return FALSE;
            error("in_header(%d): unexpected EOF at line %d",i+1,linecnt);
        }
    }

    dvals = (double *) allocate(nheader * sizeof(double));

    if (nheader > 0) {
        if (get_double(nheader,dvals) != nheader) return FALSE;
        for (i=0; i<nheader; i++) {
            if (streq(header[i],"nbody"))
                nbody = dvals[i];
            else if (streq(header[i],"time"))
                tsnap = dvals[i];
            else if (streq(header[i],"ndim"))
                ndim = dvals[i];
            else if (streq(header[i],"skip"))
                dval = dvals[i];
            else if (streq(header[i],"-"))
                dval = dvals[i];
            else
                warning("Skipping unknown header name %s",header[i]);
        }
    } else
	if (feof(instr)) return FALSE;

    if (nbody==0) {
        if (hasvalue("nbody")) nbody=getiparam("nbody");
        ndim = hasvalue("ndim") ? getiparam("ndim") : NDIM;
    }

    if (ndim>0 && ndim != NDIM)
	error("got ndim = %d, not %d", ndim, NDIM);
    return TRUE;
}

local bool get_block(int id,string options)
{
    int i, j, n, need, ngot, nvals;
    char line[MAXLINE];
    double dvals[MAXVALS];
    string *o;
    Body *bp;
    bool skip = streq(options,"skip");   /* if to skip the whole block */
    int colmass, colpos, colvel, colphase, colphi, colacc, colaux, coldens, colkey;
    int colx, coly, colz, colvx, colvy, colvz, colax, colay, colaz;

    dprintf(1,"Block%d: %s",id,options);
    
    if (options==NULL || *options==0) return FALSE;
    o = burststring(options,", \t");
    n = xstrlen(o,sizeof(string)) - 1;
    colmass=colpos=colvel=colphase=colphi=colacc=colaux=colkey=coldens=0;
    colx=coly=colz=colvx=colvy=colvz=colax=colay=colaz=0;

    need = 0;
    for (i=0; i<n; i++) {
        if (streq(o[i], "skip") || streq(o[i], "-")) {
            need += 1;
        } else if (streq(o[i], "mass")) { 
            colmass = need+1;   need += 1;      bits |= MassBit;
        } else if (streq(o[i], "m")) { 
            colmass = need+1;   need += 1;      bits |= MassBit;
        } else if (streq(o[i], "phase")){ 
            colphase= need+1;   need += 2*NDIM; bits |= PhaseSpaceBit;
	} else if (streq(o[i], "pos"))  { 
            colpos  = need+1;   need += NDIM;   bits |= PhaseSpaceBit;
	} else if (streq(o[i], "x"))  { 
            colx    = need+1;   need += 1;      bits |= PhaseSpaceBit;
	} else if (streq(o[i], "y"))  { 
            coly    = need+1;   need += 1;      bits |= PhaseSpaceBit;
	} else if (streq(o[i], "z"))  { 
            colz    = need+1;   need += 1;      bits |= PhaseSpaceBit;
	} else if (streq(o[i], "vel"))  { 
            colvel  = need+1;   need += NDIM;   bits |= PhaseSpaceBit;
	} else if (streq(o[i], "vx"))  { 
            colvx   = need+1;   need += 1;      bits |= PhaseSpaceBit;
	} else if (streq(o[i], "vy"))  { 
            colvy   = need+1;   need += 1;      bits |= PhaseSpaceBit;
	} else if (streq(o[i], "vz"))  { 
            colvz   = need+1;   need += 1;      bits |= PhaseSpaceBit;
        } else if (streq(o[i], "phi"))  { 
            colphi  = need+1;   need += 1;      bits |= PotentialBit;
        } else if (streq(o[i], "acc"))  { 
            colacc  = need+1;   need += NDIM;   bits |= AccelerationBit;
	} else if (streq(o[i], "ax"))  { 
            colax   = need+1;   need += 1;      bits |= AccelerationBit;
	} else if (streq(o[i], "ay"))  { 
            colay   = need+1;   need += 1;      bits |= AccelerationBit;
	} else if (streq(o[i], "dens"))  { 
            coldens = need+1;   need += 1;      bits |= DensBit;
        } else if (streq(o[i], "az"))  { 
	    colaz   = need+1;   need += 1;      bits |= AccelerationBit;
        } else if (streq(o[i], "aux"))  { 
            colaux  = need+1;   need += 1;      bits |= AuxBit;
        } else if (streq(o[i], "key"))  { 
            colkey  = need+1;   need += 1;      bits |= KeyBit;
        } else warning("column name %s not understood",o[i]);
    }
    if (skip)
        dprintf(1,"%d: skipping block ; need=%d\n",id,need);
    else
        dprintf(1,
     "%d: mass(%d),pos(%d),vel(%d),phase(%d),phi(%d),acc(%d),aux(%d),dens(%d),key(%d)\n",
		id,colmass,colpos,colvel,colphase,colphi,colacc,colaux,coldens,colkey);

    
    for (j=0, nvals=0, bp=btab; j<nbody; j++, bp++) {
    	while (nvals < need) {

	  do {
            if (fgets(line, MAXLINE, instr) == NULL) {
	      if (j==0 && nvals==0) {
		warning("Block %d: line %d nvals=0",id,linecnt);
		return FALSE;
	      }
	      error("process(%s): unexpected EOF for particle %d",options,j+1);
	      /* in case error is bypassed: reset nbody */
	      nbody = j;
	      warning("Resetting nbody=%d",nbody);
	      return TRUE;
            } else {
	      if (Qcom && line[0]=='#') {
		dprintf(1,"COMMENT1: %s",line);
		if (Qhis) {
		  if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
		  shist[nhist++] = strdup(line);
		}
		continue;
	      } else
		break;
	    }
	  } while(line[0]=='#' || line[0]=='\n');    /* read until EOF or non-comment lines */

	  linecnt++;
	  if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;
#if 0            
	  if (skip) continue;
#else
	  if (skip) break;            
#endif            
	  if (Qcom) if (line[0] == '#') continue;
	  tab2space(line);
	  ngot = nemoinpd(line,&dvals[nvals],MAXVALS-nvals);
	  if (ngot < 0)
	    error("(%d) Line %d Parsing %s",ngot,linecnt,line);
	  else
	    nvals += ngot;
	  if (!wrap) break;
        }
        if (skip) continue;
        if (nvals > need) {
            if (wrap) {                 /* in wrap mode */
                nvals -= need;          
                if (spill)              /* only spill if asked for */
                    nvals = 0;
            } else                      /* else always spill */
                nvals = 0;
        } else if (nvals < need) {
	  dprintf(0," bad line:%s\n",line);
	  do {
	    fgets(line, MAXLINE, instr);
	    if (Qcom && line[0]=='#') {
	      dprintf(0,"COMMENT2: %s",line);
	      if (Qhis) {
		if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
		shist[nhist++] = strdup(line);
	      }
	    }
	  } while (line[0] == '#');
	      
	  dprintf(0,"next line:%s\n",line);
            error("Bad line %d (%d) Not enough values (need %d for %s) - try wrap",
                          linecnt,nvals,need,options);
        } else  /* exact match : nvals == need */
            nvals = 0;
                    
        if (colmass>0)   Mass(bp) = dvals[colmass-1];
        if (colpos>0)    SETV(Pos(bp),&dvals[colpos-1]);
        if (colx>0)      Pos(bp)[0] = dvals[colx-1];
        if (coly>0)      Pos(bp)[1] = dvals[coly-1];
        if (colz>0)      Pos(bp)[2] = dvals[colz-1];
        if (colvel>0)    SETV(Vel(bp),&dvals[colvel-1]);
        if (colvx>0)     Vel(bp)[0] = dvals[colvx-1];
        if (colvy>0)     Vel(bp)[1] = dvals[colvy-1];
        if (colvz>0)     Vel(bp)[2] = dvals[colvz-1];
        if (colphase>0) {
                         SETV(Pos(bp),&dvals[colphase-1]);
                         SETV(Vel(bp),&dvals[colphase-1+NDIM]);
                        }
        if (colphi>0)    Phi(bp) = dvals[colphi-1];
        if (colacc>0)    SETV(Acc(bp),&dvals[colacc-1]);
        if (colaux>0)    Aux(bp) = dvals[colaux-1];
        if (coldens>0)   Dens(bp)= dvals[coldens-1];
        if (colkey>0)    Key(bp) = dvals[colkey-1];
        if (j==0) dprintf(1,"%d: First line: %s\n",id,line);
        if (j==nbody-1) dprintf(1,"%d: Last line: %s\n",id,line);
    }
    dprintf(1,"%d: Linecount=%d\n",id,linecnt);
    return TRUE;
}


/*
 * GET_DOUBLE:  get 'n' double values from a header,
 *              no need to be on one line. 
 *
 * NOTE: not sure if the line on which last number (n) was
 *       read is still accessible, since I use fgets() in
 *       reading non-header information 
 *	 The extra \n seems to have done the trick.... not sure
 *	 about portability.... ought to use cashed fgets()+nemoinp()
 */

local int get_double(int n, double *d)
{
#if 1
						/* using fgets + nemoinp */
    int i, k=0;
    char line[MAXLINE];
    while (k<n) {
      do {
        if (fgets(line, MAXLINE, instr) == NULL) {
	  if (k==0) return 0;
	  warning("Unexpected EOF in header for snapshot %d",snapcnt+1);
	  return k;
        } else {
	  if (Qcom && line[0]=='#') {
	    dprintf(0,"COMMENT3: %s\n",line);
	    if (Qhis) {
	      if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
	      shist[nhist++] = strdup(line);
	    }
	    continue;
	  } else
	    break;
	}
      } while(line[0]=='#');    /* read until EOF or non-comment line */
      if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
      tab2space(line);
      i = nemoinpd(line,&d[k],n-k);
      if (i==0) return k;
      if (i<0) error("(%d) Error parsing %s",i,line);
      k += i;
    }
    return k;
#else
						/* using (slow) fscanf */
    int i;

    if (n <= 0) return n;
    for (i=0; i<n; i++) 
        if (fscanf(instr, " %lf\n", &d[i]) != 1) return 0;
    return n;
#endif
}

/*
 * GET_NBODY:   estimate nbody if none supplied  from the
 *		number of lines. Only works if no header
 *		given of course, since we don't know the
 *		number of lines a header consumes (see prev.
 *		routine get_double()
 *              Also doesn't work well when options=wrap would be needed
 */

local int get_nbody(void)
{
    if (hasvalue("header"))
        error("Need value for nbody=, or specify it in header=");
    nbody = nemo_file_lines(getparam("in"),0);
    if (nbody <= 0) 
        error("Cannot determine nbody, try nbody= or header=");
    if (nblocks == 0) 
        error("No block's specified");
    if (ntimes > 0) {
        if (nbody % ntimes != 0)
            warning("File-length (%d) not multiple of number of times (%d)",
                    nbody, ntimes);
        nbody /= ntimes;
    }
    if (nbody % nblocks != 0)
        warning("File-length (%d) not multiple of number of blocks (%d)",
            nbody, nblocks);
    nbody /= nblocks;
    return nbody;
}

/*
 * DO_SCAN: scan file, and report on regularity in number of columns
 */

void do_scan(stream instr)
{
    char line[MAXLINE];
    string *sp;
    int ncol, ncol_old=-1;

    warning("Scanning....");
    dprintf(0,"Line\t#Columns\n");
    while (fgets(line, MAXLINE, instr) != NULL) {
        linecnt++;
        if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;        
        sp = burststring(line," ,\t");
        ncol = xstrlen(sp,sizeof(string)) - 1;
        if (ncol_old != ncol) {
            printf("%d\t%d\n",linecnt,ncol);
            ncol_old = ncol;
        }
        freestrings(sp);
    }
    printf("%d\t%d\n",linecnt,ncol);
}


/*
 * small helper function, replaces tabs by spaces before processing.
 * this prevents me from diving into gipsy parsing routines and  fix
 * the problem there
 * PJT - June 1998.	 - see also tabmath.c
 */
        
local void tab2space(char *cp)
{
    while (*cp) {
        if (*cp == '\t') *cp = ' ';
        cp++;
    }
}

