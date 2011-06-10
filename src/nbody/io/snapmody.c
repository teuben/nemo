/*
 *  convert mody's 'pos-vel' dumps to snapshot - in single precision
 *
 *    23-feb-09  PJT  created ; 
 *    25-feb-09  PJT  added pin= (but not processing) , wrote out few more pars
 */
 
#define SinglePrec
#define SINGLEPREC
#undef  DOUBLEPREC

#include <nemo.h>
#include <unfio.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input snapshot file",
    "out=???\n      Output mody moutXX.bin pos-vel file",
    "out2=\n        Output of history/headlines",
    "swap=f\n       Swap bytes in each float if endianism differs",
    "header=t\n     Use standard MODY header?",
    "VERSION=1.2\n  10-jun-2011 pjt",
    NULL,
};

string usage = "convert snapshot to Mody's pos-vel file";

extern void bswap(void *vdat, int len, int cnt);
 
void nemo_main()
{
    stream instr, outstr, out2str;
    string hl;
    real *fbuf, *mbuf, tsnap, mass0, mass1, mass2;
    bool Qswap = getbparam("swap");
    bool Qhead = getbparam("header");
    long nread;
    int nbody, i, nbad;
    int cs = CSCode(Cartesian, NDIM, 2);
    int ipar[5];
    float rpar[5], posvel[6];
    

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");


    /* read snapshot */

    get_history(instr);
    hl = ask_headline();
    get_set(instr,SnapShotTag);
      get_set(instr,ParametersTag);
       get_data(instr, NobjTag, IntType, &nbody, 0);
       get_data_coerced(instr, TimeTag, RealType, &tsnap, 0);
      get_tes(instr,ParametersTag);

      fbuf = (real *) allocate(nbody*6*sizeof(real)); 
      mbuf = (real *) allocate(nbody*sizeof(real));   

      get_set(instr,ParticlesTag);
        get_data(instr,CoordSystemTag, IntType, &cs, 0);
        get_data_coerced(instr,MassTag,RealType,mbuf,nbody,0);
        get_data_coerced(instr,PhaseSpaceTag,RealType,fbuf,nbody,2,NDIM,0);
      get_tes(instr,ParticlesTag);
    get_tes(instr,SnapShotTag);


    ipar[0] = nbody;
    ipar[1] = 0;
    ipar[2] = 0;
    /* ipar[3,4] still unused */


    mass1 = mbuf[0];
    mass0 = mass1*nbody;
    rpar[0] = mass0;
    rpar[1] = tsnap;
    rpar[2] = 0;
    /* rpar[3,4] still unused */
    for (i=0, nbad=0, mass2=0.0; i<nbody; i++) {
      mass2 += mbuf[i];
      if (mbuf[i] != mass1) nbad++;
    }
    if (nbad) 
      warning("Found %d masses not equal to first one; total=%g, expected %g",
	      nbad, mass2, mass0);
			      

    if (Qhead) {
      nread = unfwrite(outstr, (char *)ipar, 5*sizeof(int));
      nread = unfwrite(outstr, (char *)rpar, 5*sizeof(float));
    }

    for (i=0; i<nbody; i++) {
      mbuf[i] = mass1;                              
      nread = unfwrite(outstr, (char *) &fbuf[i*6], 6*sizeof(real)); 
      if (nread<1) error("Early EOF for i=%d",i);
    }
    
    strclose(instr);
    strclose(outstr);

    if (hasvalue("out2")) {
      outstr = stropen(getparam("out2"),"w");
      put_history(outstr);
      strclose(outstr);
    }
}




