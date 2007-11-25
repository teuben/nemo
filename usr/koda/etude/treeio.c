/****************************************************************************/
/* TREEIO.C: I/O routines for hierarchical N-body code. Public routines:    */
/* inputdata(), startoutput(), output(), savestate(), restorestate().       */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000           */
/*    revised for NEMO by Peter Teuben - october 2007                       */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "treecode.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>

#if defined(USE_NEMO_IO)
#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#endif

/*
 * Prototypes for local routines.
 */

local void outputdata(void);                    /* write N-body data        */
local void diagnostics(void);                   /* eval N-body diagnostics  */
local void in_int(stream, int *);               /* input integer value      */
local void in_real(stream, real *);             /* input real value         */
local void in_vector(stream, vector);           /* input vector of reals    */
local void out_int(stream, int);                /* output integer value     */
local void out_real(stream, real);              /* output real value        */
local void out_vector(stream, vector);          /* output vector of reals   */

/*
 * Diagnostic output variables.
 */

local real mtot;                                /* total mass of system     */
local real etot[3];                             /* Etot, KE, PE of system   */
local matrix keten;                             /* kinetic energy tensor    */
local matrix peten;                             /* potential energy tensor  */
local vector cmpos;                             /* center of mass position  */
local vector cmvel;                             /* center of mass velocity  */
#if defined(THREEDIM)
local vector amvec;
#else
local real amvec;                               /* angular momentum vector  */
#endif

/*
 * INPUTDATA: read initial conditions from input file.
 */

void inputdata(void)
{
    stream instr;
    int i, ndim;
    bodyptr p;
    real *bsca;
    vector *bvec;

    instr = stropen(infile, "r");               /* open input stream        */
#if defined(USE_NEMO_IO)
    get_history(instr);
    get_set(instr,SnapShotTag);
      get_set(instr,ParametersTag);
        get_data(instr,NobjTag,IntType,&nbody,0);
        if (get_tag_ok(instr,TimeTag))
          get_data(instr,TimeTag,RealType,&tnow,0);
        else
          tnow = 0.0;
      get_tes(instr,ParametersTag);
      get_set(instr,ParticlesTag);
        bsca = (real *) allocate(nbody*sizeof(real));
        bodytab = (bodyptr) allocate(nbody * sizeof(body));

        get_data(instr,MassTag,RealType,bsca,nbody,0);
        for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
          Mass(p) = bsca[i];

        if (get_tag_ok(instr,PhaseSpaceTag)) {
          bvec = (vector *) allocate(nbody*sizeof(vector)*2);
          get_data(instr,PhaseSpaceTag,RealType,bvec,nbody,2,NDIM,0);
          for (p = bodytab, i=0; p < bodytab+nbody; p++, i++) {
            SETV(Pos(p),bvec[i]);
            i++;
            SETV(Vel(p),bvec[i]);
          }       
        } else {
          bvec = (vector *) allocate(nbody*sizeof(vector));
          get_data(instr,PosTag,RealType,bvec,nbody,NDIM,0);
          for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
            SETV(Pos(p),bvec[i]);
          
          get_data(instr,VelTag,RealType,bvec,nbody,NDIM,0);
          for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
            SETV(Vel(p),bvec[i]);
        }

        get_data(instr,SmoothTag,RealType,bsca,nbody,0);
        for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
          Hknl(p) = bsca[i];

        free(bvec);
        free(bsca);
      get_tes(instr,ParticlesTag);
    get_tes(instr,SnapShotTag);

#else
    in_int(instr, &nbody);                      /* read number of bodies    */
    if (nbody < 1)
        error("inputdata: nbody = %d is absurd\n", nbody);
    in_int(instr, &ndim);                       /* read number of dims      */
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real(instr, &mscale);                    /* mass scale in Msun       */
    in_real(instr, &rscale);                    /* radius scale in kpc      */
    in_real(instr, &tscale);                    /* time scale in year       */
    in_real(instr, &tnow);                      /* read starting time       */
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
                                                /* allocate body array      */
    for (p = bodytab; p < bodytab+nbody; p++)   /* loop over all bodies     */
        in_real(instr, &Mass(p));               /* read mass of each        */
    for (p = bodytab; p < bodytab+nbody; p++)
        in_vector(instr, Pos(p));               /* read position of each    */
    for (p = bodytab; p < bodytab+nbody; p++)
        in_vector(instr, Vel(p));               /* read velocity of each    */
    for (p = bodytab; p < bodytab+nbody; p++)
        in_real(instr, &Hknl(p));               /* read radius of each      */
    fclose(instr);                              /* close input stream       */
    if (scanopt(options, "reset-time"))         /* reset starting time?     */
        tnow = 0.0;                             /* then set it to zero      */
    for (p = bodytab; p < bodytab+nbody; p++)   /* loop over new bodies     */
        Type(p) = BODY;                         /* initialize type field    */
#endif
}

/*
 * STARTOUTPUT: begin output to log file.
 */

void startoutput(void)
{
    printf("\n%s\n", headline);                 /* print headline, params   */
    printf("\n%8s%8s%8s", "nbody", "freq", "eps");
#if !defined(QUICKSCAN)
    printf("%8s", "theta");
#endif
    printf("%8s%8s%8s\n", "usequad", "freqout", "tstop");
    printf("%8d%8.2f%8.4f", nbody, freq, eps);
#if !defined(QUICKSCAN)
    printf("%8.2f", theta);
#endif
    printf("%8s%8.2f%8.4f\n", usequad ? "true" : "false", freqout, tstop);
    if (! strnull(options))                     /* print options, if any    */
        printf("\n\toptions: %s\n", options);
    if (! strnull(savefile))                    /* was state file given?    */
        savestate(savefile);                    /* save initial data        */
}

/*
 * FORCEREPORT: print staristics on tree construction and force calculation.
 */

void forcereport(void)
{
    printf("\n\t%7s%7s%7s%8s%8s%8s%8s%8s\n",
           "rsize", "tdepth", "ftree", "actmax",
	   "nbbtot", "nbctot", "actcmx", "actbmx");
    printf("\t%7.1f%7d%7.3f%8d%8d%8d%8d%8d\n",
           rsize, tdepth, (nbody + ncell - 1) / ((real) ncell),
           actmax, nbbcalc, nbccalc, actcmax, actbmax);
    printf("\n\t%10s%10s%10s%10s%10s%10s\n",
	   "CPUtree", "CPUngbr", "CPUdens", "CPUsphf", "CPUgrav", "CPUfc");
    printf("\t%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
	   cputree, cpungbr, cpudens, cpusphf, cpugrav,
	   cpungbr + cpudens + cpusphf + cpugrav );
}

/*
 * OUTPUT: compute diagnostics and output body data.
 */

void output(void)
{
    real cmabs, amabs, teff;

    diagnostics();                              /* compute std diagnostics  */
    ABSV(cmabs, cmvel);                         /* find magnitude of cm vel */
#if defined(THREEDIM)
    ABSV(amabs, amvec);                         /* find magnitude of J vect */
#else
    amabs = ABS(amvec);                         /* find magnitude of J vect */
#endif
    printf("\n %8s%8s%8s%8s%8s%8s%8s%8s\n",
           "time", "|T+U|", "T", "-U", "-T/U", "|Vcom|", "|Jtot|", "CPUtot");
    if (selfgrav)
      printf(" %8.3f%8.5f%8.5f%8.5f%8.5f%8.5f%8.5f%8.3f\n",
	     tnow, ABS(etot[0]), etot[1], -etot[2], -etot[1]/etot[2],
	     cmabs, amabs, cputime());
    else
      printf(" %8.3f%8.5f%8.5f%8s%8s%8.5f%8.5f%8.3f\n",
	     tnow, ABS(etot[0]), etot[1], "NONE", "NONE",
	     cmabs, amabs, cputime());
    teff = tnow + (freq > 0 ? 0.125/freq : 0);  /* anticipate slightly...   */
    if (! strnull(outfile) && teff >= tout)     /* time for data output?    */
        outputdata();
    if (! strnull(savefile))                    /* was state file given?    */
        savestate(savefile);                    /* save data for restart    */
}

/*
 * OUTPUTDATA: output body data.
 */

void outputdata(void)
{
    static int iout=-1;
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;
    int i;
    static int first=1;
    real *bsca;
    vector *bvec, cm[2];
 

    sprintf(namebuf, "%s_%.3d", outfile, ++iout);
                                                /* construct output name    */
    if (iout == 0 && !strnull(infile))          /* if start from halfway    */
	while (stat(namebuf, &buf) == 0)        /* increment file number    */
	    sprintf(namebuf, "%s_%.3d", outfile, ++iout);
    outstr = stropen(namebuf, "w");             /* create & open for output */
#if defined(USE_NEMO_IO)
    if (first)
      put_history(outstr);
    put_set(outstr,SnapShotTag);
      put_set(outstr,ParametersTag);
        put_data(outstr,NobjTag,IntType,&nbody,0);
        put_data(outstr,TimeTag,RealType,&tnow,0);
      put_tes(outstr,ParametersTag);
      put_set(outstr,ParticlesTag);
        bsca = (real *) allocate(nbody*sizeof(real));
        bvec = (vector *) allocate(nbody*sizeof(vector));
        for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
          bsca[i] = Mass(p);
        if (first)
          put_data(outstr,MassTag,RealType,bsca,nbody,0);
        for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
          SETV(bvec[i],Pos(p));
        put_data(outstr,PosTag,RealType,bvec,nbody,NDIM,0);
        for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
          SETV(bvec[i],Vel(p));
        put_data(outstr,VelTag,RealType,bvec,nbody,NDIM,0);
        if (scanopt(options, "out-phi")) {
          for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
            bsca[i] = Phi(p);
          put_data(outstr,PotentialTag,RealType,bsca,nbody,0);
        }
        if (scanopt(options, "out-acc")) {
          for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
            SETV(bvec[i],Acc(p));
          put_data(outstr,AccelerationTag,RealType,bvec,nbody,NDIM,0);
        }
        for (p = bodytab, i=0; p < bodytab+nbody; p++, i++)
          bsca[i] = Hknl(p);
	put_data(outstr,SmoothTag,RealType,bsca,nbody,0);
      put_tes(outstr,ParticlesTag);
      put_set(outstr,DiagnosticsTag);
        put_data(outstr,EnergyTag,RealType,etot,3,0);
        SETV(cm[0],cmpos);
        SETV(cm[1],cmvel);
        put_data(outstr,CMPhaseSpaceTag,RealType,cm,2,NDIM,0);
      put_tes(outstr,DiagnosticsTag);
    put_tes(outstr,SnapShotTag);
    if (first) first = 0;
    free(bsca);
    free(bvec);
    strclose(outstr);
#else
    out_int(outstr, nbody);                     /* write number of bodies   */
    out_int(outstr, NDIM);                      /* number of dimensions     */
    out_real(outstr, mscale);                   /* mass scale in Msun       */
    out_real(outstr, rscale);                   /* radius scale in kpc      */
    out_real(outstr, tscale);                   /* time scale in year       */
    out_real(outstr, tnow);                     /* and current time value   */
    for (p = bodytab; p < bodytab+nbody; p++)   /* loop over all bodies     */
        out_real(outstr, Mass(p));              /* output mass of each      */
    for (p = bodytab; p < bodytab+nbody; p++)
        out_vector(outstr, Pos(p));             /* output positions         */
    for (p = bodytab; p < bodytab+nbody; p++)
        out_vector(outstr, Vel(p));             /* output velocities        */
    for (p = bodytab; p < bodytab+nbody; p++)
	out_real(outstr, Hknl(p));              /* output SPH radius        */
    if (scanopt(options, "out-phi"))            /* potentials requested?    */
        for (p = bodytab; p < bodytab+nbody; p++)
            out_real(outstr, Phi(p));           /* output potentials        */
    if (scanopt(options, "out-acc"))            /* accelerations requested? */
        for (p = bodytab; p < bodytab+nbody; p++)
            out_vector(outstr, Acc(p));         /* output accelerations     */
    fclose(outstr);                             /* close up output file     */
#endif
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
    tout += 1.0 / freqout;                      /* schedule next output     */
}

/*
 * DIAGNOSTICS: compute set of dynamical diagnostics.
 */

local void diagnostics(void)
{
    register bodyptr p;
    real velsq, tmpr;
    vector tmpv;
    matrix tmpt;

    mtot = 0.0;                                 /* zero total mass          */
    etot[1] = etot[2] = 0.0;                    /* zero total KE and PE     */
    CLRM(keten);                                /* zero ke tensor           */
    CLRM(peten);                                /* zero pe tensor           */
    CLRV(cmpos);                                /* zero c. of m. position   */
    CLRV(cmvel);                                /* zero c. of m. velocity   */
#if defined(THREEDIM)
    CLRV(amvec);                                /* zero am vector           */
#else
    amvec=0.0;                                  /* zero am vector           */
#endif
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all particles  */
        mtot += Mass(p);                        /* sum particle masses      */
        DOTVP(velsq, Vel(p), Vel(p));           /* square vel vector        */
        etot[1] += 0.5 * Mass(p) * velsq;       /* sum current KE           */
        etot[2] += 0.5 * Mass(p) * Phi(p);      /* and current PE           */
        MULVS(tmpv, Vel(p), 0.5 * Mass(p));     /* sum 0.5 m v_i v_j        */
        OUTVP(tmpt, tmpv, Vel(p));
        ADDM(keten, keten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));           /* sum m r_i a_j            */
        OUTVP(tmpt, tmpv, Acc(p));
        ADDM(peten, peten, tmpt);
#if defined(THREEDIM)
        CROSSVP(tmpv, Vel(p), Pos(p));          /* sum angular momentum     */
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec, amvec, tmpv);
#else
        CROSSVP(tmpr, Vel(p), Pos(p));          /* sum angular momentum     */
        tmpr = tmpr * Mass(p);
        amvec = amvec + tmpr;
#endif
        MULVS(tmpv, Pos(p), Mass(p));           /* sum cm position          */
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           /* sum cm momentum          */
        ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2];                /* sum KE and PE            */
    DIVVS(cmpos, cmpos, mtot);                  /* normalize cm coords      */
    DIVVS(cmvel, cmvel, mtot);
}

/*
 * IN_INT, IN_REAL, IN_VECTOR: low level input routines.
 */

local void in_int(stream str, int *iptr)
{
#if !defined(BINARYIO)
    if (fscanf(str, "%d", iptr) != 1)
        error("in_int: input conversion error\n");
#else
    if (fread((void *) iptr, sizeof(int), 1, str) != 1)
        error("in_int: fread failed\n");
#endif
}

local void in_real(stream str, real *rptr)
{
    double tmp;

#if !defined(BINARYIO)
    if (fscanf(str, "%lf", &tmp) != 1)
        error("in_real: input conversion error\n");
    *rptr = tmp;
#else
    if (fread((void *) rptr, sizeof(real), 1, str) != 1)
        error("in_real: fread failed\n");
#endif
}

local void in_vector(stream str, vector vec)
{
    double tmpx, tmpy, tmpz;

#if !defined(BINARYIO)
    if (fscanf(str, "%lf%lf%lf", &tmpx, &tmpy, &tmpz) != 3)
        error("in_vector: input conversion error\n");
    vec[0] = tmpx;
    vec[1] = tmpy;
    vec[2] = tmpz;
#else
    if (fread((void *) vec, sizeof(real), NDIM, str) != NDIM)
        error("in_vector: fread failed\n");
#endif
}

/*
 * OUT_INT, OUT_REAL, OUT_VECTOR: low level output routines.
 */

#define IFMT  " %d"                             /* output format for ints   */
#define RFMT  " %14.7E"                         /* output format for reals  */

local void out_int(stream str, int ival)
{
#if !defined(BINARYIO)
    if (fprintf(str, IFMT "\n", ival) < 0)
        error("out_int: fprintf failed\n");
#else
    if (fwrite((void *) &ival, sizeof(int), 1, str) != 1)
        error("out_int: fwrite failed\n");
#endif
}

local void out_real(stream str, real rval)
{
#if !defined(BINARYIO)
    if (fprintf(str, RFMT "\n", rval) < 0)
        error("out_real: fprintf failed\n");
#else
    if (fwrite((void *) &rval, sizeof(real), 1, str) != 1)
        error("out_real: fwrite failed\n");
#endif
}

local void out_vector(stream str, vector vec)
{
#if !defined(BINARYIO)
    if (fprintf(str, RFMT RFMT RFMT "\n", vec[0], vec[1], vec[2]) < 0)
        error("out_vector: fprintf failed\n");
#else
    if (fwrite((void *) vec, sizeof(real), NDIM, str) != NDIM)
        error("out_vector: fwrite failed\n");
#endif
}

/*
 * SAVESTATE: write current state to disk file.
 */

#define safewrite(ptr,len,str)                  \
    if (fwrite((void *) ptr, len, 1, str) != 1) \
        error("savestate: fwrite failed\n")

void savestate(string pattern)
{
    char namebuf[256];
    stream str;
    int nchars;

    sprintf(namebuf, pattern, nstep & 1);       /* construct alternate name */
    str = stropen(namebuf, "w!");
    nchars = strlen(getargv0()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getargv0(), nchars * sizeof(char), str);
    nchars = strlen(getversion()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getversion(), nchars * sizeof(char), str);
    safewrite(&freq, sizeof(real), str);
#if !defined(QUICKSCAN)
    safewrite(&theta, sizeof(real), str);
#endif
    safewrite(&usequad, sizeof(bool), str);
    safewrite(&eps, sizeof(real), str);
    nchars = strlen(options) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(options, nchars * sizeof(char), str);
    safewrite(&tstop, sizeof(real), str);
    safewrite(&freqout, sizeof(real), str);
    safewrite(&tnow, sizeof(real), str);
    safewrite(&tout, sizeof(real), str);
    safewrite(&nstep, sizeof(int), str);
    safewrite(&rsize, sizeof(real), str);
    safewrite(&nbody, sizeof(int), str);
    safewrite(bodytab, nbody * sizeof(body), str);
    fclose(str);
}

/*
 * RESTORESTATE: restore state from disk file.
 */

#define saferead(ptr,len,str)                  \
    if (fread((void *) ptr, len, 1, str) != 1) \
        error("restorestate: fread failed\n")

void restorestate(string file)
{
    stream str;
    int nchars;
    string program, version;

    str = stropen(file, "r");
    saferead(&nchars, sizeof(int), str);
    program = (string) allocate(nchars * sizeof(char));
    saferead(program, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    version = (string) allocate(nchars * sizeof(char));
    saferead(version, nchars * sizeof(char), str);
    if (! streq(program, getargv0()) ||         /* check program, version   */
          ! streq(version, getversion()))
        printf("warning: state file may be outdated\n\n");
    saferead(&freq, sizeof(real), str);
#if !defined(QUICKSCAN)
    saferead(&theta, sizeof(real), str);
#endif
    saferead(&usequad, sizeof(bool), str);
    saferead(&eps, sizeof(real), str);
    saferead(&nchars, sizeof(int), str);
    options = (string) allocate(nchars * sizeof(char));
    saferead(options, nchars * sizeof(char), str);
    saferead(&tstop, sizeof(real), str);
    saferead(&freqout, sizeof(real), str);
    saferead(&tnow, sizeof(real), str);
    saferead(&tout, sizeof(real), str);
    saferead(&nstep, sizeof(int), str);
    saferead(&rsize, sizeof(real), str);
    saferead(&nbody, sizeof(int), str);
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
    saferead(bodytab, nbody * sizeof(body), str);
    fclose(str);
}
