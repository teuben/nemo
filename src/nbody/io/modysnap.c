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
    "in=???\n       Input mody moutXX.bin pos-vel file",
    "out=???\n      Output snapshot file",
    "pin=\n         Input mbody poutXX.bin potential file (optional)",
    "swap=f\n       Swap bytes in each float if endianism differs",
    "headline=\n    Random verbiage",
    "VERSION=1.1\n  25-feb-09 pjt",
    NULL,
};

string usage = "convert Mody's pos-vel files to snapshot";

extern void bswap(void *vdat, int len, int cnt);
 
void nemo_main()
{
    stream instr, outstr;
    real *fbuf, *mbuf, tsnap, mass0, mass1, tdyn;
    bool Qswap = getbparam("swap");
    long nread;
    int nbody, i, mond_ind, model;
    int cs = CSCode(Cartesian, NDIM, 2);
    int ipar[5];
    float rpar[5], posvel[6];
    

    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    put_history(outstr);

    if (hasvalue("pin")) 
      warning("pin=%s not yet processed",getparam("pin"));


    nread = unfread(instr, (char *)ipar, 5*sizeof(int));
    dprintf(0,"ipar: nbody=%d mods=%d\n",
	    ipar[0],ipar[1]);

    nread = unfread(instr, (char *)rpar, 5*sizeof(float));
    dprintf(0,"rpar: totms=%f tnow=%f tdyn=%f vir=%f\n",
	    rpar[0],rpar[1],rpar[2],rpar[3]);

    nbody = ipar[0];
    model = ipar[1];
    mond_ind = ipar[2];
    /* ipar[3,4] still unused */

    mass0 = rpar[0];
    mass1 = mass0/nbody;
    tsnap = rpar[1];
    tdyn  = rpar[2];
    /* rpar[3,4] still unused */

    /* fbuf holds the array of posvel's
     * mbuf holds the (all the same ) masses 
     */
    fbuf = (real *) allocate(nbody*6*sizeof(real)); 
    mbuf = (real *) allocate(nbody*sizeof(real));   

    for (i=0; i<nbody; i++) {
      mbuf[i] = mass1;                              
      nread = unfread(instr, (char *) &fbuf[i*6], 6*sizeof(real)); 
      if (nread<1) error("Early EOF for i=%d",i);
    }
    
    /* write out snapshot for NEMO */

    put_set(outstr,SnapShotTag);
      put_set(outstr,ParametersTag);
       put_data(outstr, NobjTag, IntType, &nbody, 0);
       put_data(outstr, TimeTag, RealType, &tsnap, 0);
       put_data(outstr, "mass",      RealType, &mass0, 0);
       put_data(outstr, "tdyn",      RealType, &tdyn, 0);
       put_data(outstr, "model",     IntType, &model, 0);
       put_data(outstr, "mond_ind",  IntType, &mond_ind, 0);
      put_tes(outstr,ParametersTag);
      put_set(outstr,ParticlesTag);
        put_data(outstr,CoordSystemTag, IntType, &cs, 0);
        put_data(outstr,MassTag,RealType,mbuf,nbody,0);
        put_data(outstr,PhaseSpaceTag,RealType,fbuf,nbody,2,NDIM,0);
      put_tes(outstr,ParticlesTag);
    put_tes(outstr,SnapShotTag);

    strclose(instr);
    strclose(outstr);
}




