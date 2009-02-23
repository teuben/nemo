/*
 *  convert mody's 'pos-vel' dumps to snapshot - in single precision
 *    23-feb-09    created
 */
 
#define SinglePrec
#define SINGLEPREC
#undef  DOUBLEPREC

#include <nemo.h>
#include <unfio.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input mody pos-vel file",
    "out=???\n      Output snapshot file",
    "headline=\n    Random verbiage",
    "swap=f\n       Swap bytes in each float if endianism differs",
    "VERSION=1.0\n  23-feb-09 pjt",
    NULL,
};

string usage = "convert Mody's pos-vel files to snapshot";

extern void bswap(void *vdat, int len, int cnt);
 
void nemo_main()
{
    stream instr, outstr;
    real *fbuf, *mbuf, tsnap, mass1;
    bool Qswap = getbparam("swap");
    long nread;
    int nbody, i,j;
    int cs = CSCode(Cartesian, NDIM, 2);
    int ipar[5];
    float rpar[5], posvel[6];
    

    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    put_history(outstr);


    nread = unfread(instr, (char *)ipar, 5*sizeof(int));
    dprintf(0,"ipar: %d  nbody=%d mods=%d\n",nread,ipar[0],ipar[1]);

    nread = unfread(instr, (char *)rpar, 5*sizeof(float));
    dprintf(0,"rpar: %d  totms=%f tnow=%f tdyn=%f vir=%f\n",nread,rpar[0],rpar[1],rpar[2],rpar[3]);

    nbody = ipar[0];
    tsnap = rpar[1];
    mass1 = rpar[0]/nbody;

    fbuf = (real *) allocate(nbody*6*sizeof(real));
    mbuf = (real *) allocate(nbody*sizeof(real));

    for (i=0; i<nbody; i++) {
      mbuf[i] = mass1;
      nread = unfread(instr, (char *) &fbuf[i*6], 6*sizeof(real));
      if (nread<1) error("Early EOF for i=%d",i);
    }
    
    put_set(outstr,SnapShotTag);
    put_set(outstr,ParametersTag);
    put_data(outstr, NobjTag, IntType, &nbody, 0);
    put_data(outstr, TimeTag, RealType, &tsnap, 0);
    put_tes(outstr,ParametersTag);
    put_set(outstr,ParticlesTag);
    put_data(outstr,CoordSystemTag, IntType, &cs, 0);
    put_data(outstr,MassTag,RealType,mbuf,nbody,0);
    put_data(outstr,PhaseSpaceTag,RealType,fbuf,nbody,2,NDIM,0);
    put_tes(outstr,ParticlesTag);
    put_tes(outstr,SnapShotTag);

    free(fbuf);
    strclose(instr);
    strclose(outstr);
}




