/*
 *  convert Heller SPH data to standard (single precision) snapshot
 *
 *  25-jun-98   PJT     Created         
 *  22-may-01   PJT     Updated for new format that includes star formation
 *                      (done while in Mexico w/ Heller et al)
 *  26-may-01   PJT     added aux= (alas not implemented yet)
 */
 
#define SINGLEPREC

#include <nemo.h>
#include <memio.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input rvc file",
    "out=???\n      Output snapshot file",
    "lsfm=\n        Override contents of lsfm boolean in data",
    "aux=\n         Select something to be written to the 'aux' slot",
    "headline=\n    Random verbiage",
    "VERSION=2.2\n  26-may-01 pjt",
    NULL,
};

string usage = "convert Heller's SPH dumps to snapshot";


local void scatter(int n, real *from, real *to, int skip);
 
nemo_main()
{
    stream instr, outstr;
    char buf[400], date[25];    /* character*24 */
    char *cp = buf;
    int *ip = (int *) cp;
    float *fp = (float *) cp;
    double *dp = (double *) buf;
    long nread, indx, foffset;
    int i, nbody, ndim, dmpindx, n1, n2, *key, lsfm, block;
    real *phase, *pp, tsnap, *mass, *aux, *phi, *acc;
    real *sfmt, *sfmz;
    float *rv, *rvp, *rbuf;
    double version;
    int pjt_version = 2;
    bool Qlsfm;
    string aux_slot = getparam("aux");	/* not used yet */

    Qlsfm = hasvalue("lsfm");
    if (Qlsfm)
      lsfm = getbparam("lsfm") ? 1 : 0;   /* care, fortran could use something else */

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    if (hasvalue("headline"))
        set_headline(getparam("headline"));
    put_history(outstr);

    for (;;) {    			/* loop over Heller snapshots */
        block = 1;

        /* header 1 */ 
	dprintf(1,"Reading %d header-1\n",block++);
        nread = unfread(instr, buf, 400);
        if (nread<1) break;

        dmpindx = ip[0];
	strncpy(date,&buf[4],24); date[24] = 0;
        ndim = ip[7];    
        if (ndim != NDIM)
            error("No support for ndim=%d (NDIM=%d)\n",ndim,NDIM);

	version =  *((double *) &fp[9]);
	if (!Qlsfm)
	  lsfm = ip[11];     /* fortran boolean */

        /* header 2 */
	dprintf(1,"Reading %d header-2\n",block++);
        nread = unfread(instr, buf, 400);
        if (nread<1) error("cannot read 2nd header");
        
        dprintf(0,"Time %g N1 %d N2 %d  DMPINDX=%d  NDIM=%d version=%g lsfm=%d\n",
               dp[0],  ip[2], ip[3], dmpindx, ndim, version, lsfm);

        n1 = ip[2];
        n2 = ip[3];            
        nbody = n1 + n2;
        tsnap = dp[0];
        
        rbuf = (real *) allocate(nbody*sizeof(real));

        key = (int *) allocate(nbody*sizeof(int));        

	dprintf(1,"Reading %d CLASS\n",block++);
        nread = unfread(instr, key, nbody * sizeof(int));
        if (nread<1) error("cannot read CLASS");

	dprintf(1,"Reading %d STATE\n",block++);
        nread = unfread(instr, rbuf, nbody * sizeof(real));     /* ignore */
        if (nread<1) error("cannot read STATE");

        mass = (real *) allocate(nbody*sizeof(real));
	dprintf(1,"Reading %d MASS\n",block++);
        nread = unfread(instr, mass, nbody * sizeof(real));
        if (nread<1) error("cannot read MASS");

        phase = (real *) allocate(NDIM*2*nbody*sizeof(real));

        for (i=0; i<ndim; i++) {
	    dprintf(1,"Reading %d POS(%d)\n",block++,i+1);
            nread = unfread(instr, rbuf, nbody * sizeof(real));
            if (nread<1) error("cannot read POS(%d)",i+1);
            scatter(nbody,rbuf,&phase[i],2*ndim);
        }
        for (i=0; i<ndim; i++) {
	    dprintf(1,"Reading VEL(%d)\n",i+1);
            nread = unfread(instr, rbuf, nbody * sizeof(real));
            if (nread<1) error("cannot %d read VEL(%d)",block++,i+1);
            scatter(nbody,rbuf,&phase[i+ndim],2*ndim);            
        }

        if (dmpindx > 1) {

	    dprintf(1,"Reading %d EPS\n",block++);
            nread = unfread(instr, rbuf, nbody * sizeof(real));
            if (nread<1) error("cannot read EPS");

            phi = (real *) allocate(nbody*sizeof(real));        
	    dprintf(1,"Reading %d PHI\n",block++);
            nread = unfread(instr, phi, nbody * sizeof(real));
            if (nread<1) error("cannot read PHI");

            acc = (real *) allocate(NDIM*nbody*sizeof(real));
            for (i=0; i<ndim; i++) {
	        dprintf(1,"Reading %d ACC(%d)\n",block++,i+1);
                nread = unfread(instr, rbuf, nbody * sizeof(real));
                if (nread<1) error("cannot read ACC(%d)",i+1);
                scatter(nbody,rbuf,&acc[i],ndim);                
            }
	    if (lsfm) {
	      /* read them, but don't do anything with them */
	      sfmt = (real *) allocate(NDIM*nbody*sizeof(real));
	      sfmz = (real *) allocate(NDIM*nbody*sizeof(real));
	      dprintf(1,"Reading %d SFMT\n",block++);
	      nread = unfread(instr, sfmt, nbody * sizeof(real));
	      if (nread<1) error("cannot read SFMT");
	      dprintf(1,"Reading %d SFMZ\n",block++);
	      nread = unfread(instr, sfmz, nbody * sizeof(real));
	      if (nread<1) error("cannot read SFMZ");
	    }

        }

        if (dmpindx > 2 && n2 > 0) {

	  /* warning("Discarding U, D, DIVV, UDOT1, UDOT2, H"); */
	    dprintf(1,"Reading %d U\n",block++);
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read U");

	    dprintf(1,"Reading %d D\n",block++);
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read D");

	    dprintf(1,"Reading %d DIVV\n",block++);
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read DIVV");

	    dprintf(1,"Reading %d UDOT1\n",block++);
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read UDOT1");

	    dprintf(1,"Reading %d UDOT2\n",block++);            
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read UDOT2");

	    if (pjt_version > 1) {
   	      dprintf(1,"Reading %d UDOT3\n",block++);
	      nread = unfread(instr, rbuf, n2 * sizeof(real));
	      if (nread<1) error("cannot read UDOT3");

  	      dprintf(1,"Reading %d UDOT4\n",block++);
	      nread = unfread(instr, rbuf, n2 * sizeof(real));
	      if (nread<1) error("cannot read UDOT4");
	    }

	    dprintf(1,"Reading %d H\n",block++);
            nread = unfread(instr, rbuf, n2 * sizeof(real));
            if (nread<1) error("cannot read H");

	    if (pjt_version > 1) {
	      dprintf(1,"Reading %d WGTMOL\n",block++);
	      nread = unfread(instr, rbuf, n2 * sizeof(real));
	      if (nread<1) error("cannot read WGTMOL");

   	      dprintf(1,"Reading %d FRCNEU\n",block++);
	      nread = unfread(instr, rbuf, n2 * sizeof(real));
	      if (nread<1) error("cannot read FRCNEU");

	      dprintf(1,"Reading %d TAU\n",block++);
	      nread = unfread(instr, rbuf, n2 * sizeof(real));
	      if (nread<1) error("cannot read TAU");
	    }

        }

	/* ok, done, now write selected pieces to a NEMO snapshot */

        put_set(outstr,SnapShotTag);
        put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr,ParametersTag);
        put_set(outstr,ParticlesTag);
        put_data(outstr,MassTag,RealType,mass,nbody,0);
        put_data(outstr,PhaseSpaceTag,RealType,phase,nbody,2,NDIM,0);
        put_data(outstr,KeyTag,IntType,key,nbody,0);
        if (dmpindx > 1) {
            put_data(outstr,PotentialTag,RealType,phi,nbody,0);
            put_data(outstr,AccelerationTag,RealType,acc,nbody,NDIM,0);
        }
        put_tes(outstr,ParticlesTag);
        put_tes(outstr,SnapShotTag);
        
        
    }
}



local void scatter(int n, real *from, real *to, int skip)
{
    int i, j;

    for (i=0, j=0; i<n; i++, j+=skip)
        to[j] = from[i];
}
