/*
 *  convert Martel 'pos-vel' dumps to snapshot - in single precision
 *    19-mar-99    written
 *    30-jun-01    optional swap, in case endianisms are different
 *                 also added coordsystem
 *     1-jul-01    1.1a ieck, data  are not stored (pos,vel) per particle
 *                 but first all X's, then all Y's, etc.etc.
 */
 
#define SinglePrec
#define SINGLEPREC
#undef  DOUBLEPREC

#include <nemo.h>
#include <unfio.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input pos-vel file",
    "out=???\n      Output snapshot file",
    "nbody=262144\n Maximum number of particles",
    "headline=\n    Random verbiage",
    "swap=f\n       Swap bytes in each float if endianism differs",
    "VERSION=1.1a\n 1-jul-01 pjt",
    NULL,
};

string usage = "convert Martel \"pos-vel\" dump files to snapshot";

extern void bswap(void *vdat, int len, int cnt);
 
void nemo_main()
{
    stream instr, outstr;
    real *fbuf, *fnewbuf, tsnap;
    bool Qswap = getbparam("swap");
    long nread;
    int nbody, nmax, i,j;
    int cs = CSCode(Cartesian, NDIM, 2);

    nmax = getiparam("nbody");
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    put_history(outstr);

    fbuf = (real *) allocate(nmax*6*sizeof(real));

    unfswap(Qswap);
    for (;;) {    
        nread = unfread(instr, (char *)fbuf, nmax*6*sizeof(real));
	if (Qswap) bswap((char *)fbuf, sizeof(real), nread/sizeof(real));
        if (nread<1) break;
	fnewbuf = (real *) allocate(nread*sizeof(real));
        nbody = nread/6/sizeof(real);
        if (nread % (6*sizeof(real)))
            error("nbody=%d does not divide into %s",nbody,6*sizeof(real));
	dprintf(1,"Found %d particles\n",nbody);
	for (i=0; i<6; i++) {
	  for (j=0; j<nbody; j++) {
	    fnewbuf[i+j*6] = fbuf[j+i*nbody];
	  }
	}
        tsnap = 0.0;
        put_set(outstr,SnapShotTag);
        put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr,ParametersTag);
        put_set(outstr,ParticlesTag);
	put_data(outstr,CoordSystemTag, IntType, &cs, 0);
        put_data(outstr,PhaseSpaceTag,RealType,fnewbuf,nbody,2,NDIM,0);
        put_tes(outstr,ParticlesTag);
        put_tes(outstr,SnapShotTag);
	free(fnewbuf);
    }
    free(fbuf);
    strclose(instr);
    strclose(outstr);
}




