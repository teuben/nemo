/*
 *  convert xyp (Quinn/Balcells) format to snapshot
 *      25-may-94 initial version
 *      30-may-94 added nobytes to handle single/double prec
 *	 7-jun-94 option to output more info , also added CoordSys.
 */
 
#include <nemo.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "in=???\n       Input xyp file",
    "out=\n         Output snapshot file",
    "model=\n       Select by model numbers (0=all, 1=first)",
    "times=all\n    Select by time ranges",
    "headline=\n    Random verbiage",
    "nbody=\n       Max bodies allowed to allocate",
    "nobytes=4\n    read single (4) or double (8) precision files",
    "header=f\n     Verbose output of the header?",
    "VERSION=1.2\n  7-jun-94 pjt",
    NULL,
};

string usage = "convert xvc files to snapshot";

#ifndef MAXMODELS
#define MAXMODELS 1024
#endif

/*	offsets into the header to get at header elements */

#define XVP_NBODY  0
#define XVP_ITER   1
#define XVP_T      2
#define XVP_E      3
#define XVP_J      4
#define XVP_M	   5
#define XVP_GRAVC  7
#define XVP_EPS    8

#define XVP_SATEL   17
#define XVP_NDIM    18
#define XVP_TWOGAL  19

#define XVP_MODE    99
#define XVP_NMAS   100
 
int my_read(), nbody_guess();

nemo_main()
{
    stream outstr;
    string infile, times;
    bool Qheader;
    char buf[128*4];
    char *cp = buf;
    int *ip = (int *) cp;
    float *fp = (float *) cp;
    long nread,indx,foffset;
    int i, nbody, maxbody, ierr, nobytes;
    int models[MAXMODELS], nmodels, model, idx, mode, ndim, coordsys;
    float *phase, *pp, tsnap, gravc;
    float *x, *y, *z, *vx, *vy, *vz, *mass, *pot;

    nobytes = getiparam("nobytes");
    Qheader = getbparam("header");
    infile = getparam("in");
    if (hasvalue("out")) {
        outstr = stropen(getparam("out"),"w");
        if (hasvalue("headline"))
            set_headline(getparam("headline"));
        put_history(outstr);
    } else
        outstr = NULL;
    if (hasvalue("nbody"))
        maxbody = getiparam("nbody");
    else {
        maxbody = nbody_guess(infile,nobytes);
        dprintf(0,"Using nbody=%d as maximum size of snapshots\n",maxbody);
    }
    if (hasvalue("model")) {
        nmodels = nemoinpi(getparam("model"),models,MAXMODELS);
        if (nmodels<1) error("parsing model=%s",getparam("model"));
    } else
        nmodels = 0;
    x = (float *) allocate(maxbody*sizeof(float));
    y = (float *) allocate(maxbody*sizeof(float));
    z = (float *) allocate(maxbody*sizeof(float));
    vx = (float *) allocate(maxbody*sizeof(float));
    vy = (float *) allocate(maxbody*sizeof(float));
    vz = (float *) allocate(maxbody*sizeof(float));
    mass = (float *) allocate(maxbody*sizeof(float));
    pot  = (float *) allocate(maxbody*sizeof(float));
    coordsys = CSCode(Cartesian, NDIM, 2);    	/* NDIM == 3 */

    for (idx=1;;idx++) {      /* loop, while reading snapshots */
        if (nmodels>0) {
            if (idx > nmodels) break;
            model = models[idx-1];
            if (model<0) {
                warning("Illegal model %d requested",model);
                break;
            }
        } else
            model = idx;
        ierr = my_read(nobytes,model,infile,buf,x,y,z,vx,vy,vz,mass,pot);
        if (ierr != 0) break;
        nbody = fp[XVP_NBODY];
        mode  = fp[XVP_MODE];
        ndim = fp[XVP_NDIM];
        if (ndim!=0 && ndim!=3) warning("Cannot handle ndim=%d\n",ndim);
        
        my_header(Qheader, model, buf);

        if (outstr==NULL) continue;
                  
        pp = phase = (float *) allocate(NDIM*2*nbody*sizeof(float));
        for (i=0; i<nbody; i++) {
            *pp++ = x[i];    
            *pp++ = y[i];    
            *pp++ = z[i];
            *pp++ = vx[i];   
            *pp++ = vy[i];   
            *pp++ = vz[i];   
        }
        put_set(outstr,SnapShotTag);
        put_set(outstr,ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, FloatType, &tsnap, 0);
        put_data(outstr, "GravConst", FloatType, &gravc, 0);
        put_tes(outstr,ParametersTag);
        put_set(outstr,ParticlesTag);
        put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);        
        put_data(outstr,MassTag,FloatType,mass,nbody,0);
        put_data(outstr,PhaseSpaceTag,FloatType,phase,nbody,2,NDIM,0);
        if (mode) put_data(outstr,PotentialTag,FloatType,pot,nbody,0);
        put_tes(outstr,ParticlesTag);
        put_tes(outstr,SnapShotTag);
        free(phase);
    }
}


my_header(bool verbose, int model, float *header)
{
    int nbody, iter, mode, satel, ndim, twogal, nmasses;
    int i, n, idx;
    float mass, gmass, cmass = 0.0;

    nbody = header[XVP_NBODY];
    iter  = header[XVP_ITER];
    satel = header[XVP_SATEL];
    ndim  = header[XVP_NDIM];
    twogal= header[XVP_TWOGAL];
    mode  = header[XVP_MODE];
    nmasses = header[XVP_MODE+1];

    
    if (verbose) {
        
        printf("Model %d: nbody=%d iter=%d time=%g \n",
              model,nbody,iter,header[XVP_T]);
        printf(" Total mass=%g energy=%g ang.mom=%g\n",
              header[XVP_M], header[XVP_E], header[XVP_J]);
        printf(" Gravitational softening: %g\n",header[XVP_EPS]);
        printf(" Gravitational constant:  %g\n",header[XVP_GRAVC]);

        if (satel)
            printf(" Satellite is present (%d)\n",satel);
        if (ndim)
            printf(" ndim = %d\n",ndim);
        if (twogal) {
            n = header[XVP_TWOGAL+1] - 1;
            printf(" Two-galaxies model: %d\n",twogal);
            printf("  Gal1: nbody=%d mass=%g r0=%g rmax=%g\n", n, 
                header[XVP_TWOGAL+2],
                header[XVP_TWOGAL+3],
                header[XVP_TWOGAL+4]);
            printf("  Gal2: nbody=%d mass=%g r0=%g rmax=%g\n", nbody-n, 
                header[XVP_TWOGAL+8],
                header[XVP_TWOGAL+9],
                header[XVP_TWOGAL+10]);

        }
        if (mode) {
            printf(" XVP mode: %d\n",mode);
            printf(" Number of mass groups: %d\n",nmasses);
            for (i=0, idx=XVP_MODE+2; i<nmasses; i++, idx += 2) {
                if (i==0)
                    n = header[idx];
                else
                    n = header[idx] - n;
                mass = header[idx+1];
                gmass = n*mass;
                cmass += gmass;
                printf("  Group %d: %d   Mass: %g  Group Mass: %g  Cum Mass: %g\n",
                    i+1, n, mass, gmass, cmass);                            
            }
            n = header[idx-2];
            if (n != nbody) warning("%d mass groups nbody=%d n(last)=%d",
                nmasses,nbody,n);
            
            
        }
        
       
    } else {
       dprintf(0,
	  "model %d nbody %d iter %d time %g etot %g am %g mass %g mode %d\n",
            model, nbody, iter,
            header[XVP_T], header[XVP_E], header[XVP_J], header[XVP_M],
            mode);
    }
}
           



/*
 *  SUBROUTINE xvpread(model,filename,header,
 *                     x, y, z, vx, vy, vz, pmass, ppot,
 *                     ierr)
 */
 
int my_read(
    int nobytes,        /* 4 or 8, for single & double prec data */
    int model,          /* 1=first one */
    string filename,    /* filename */
    float *header,      /* secret header */
    float *x,
    float *y, 
    float *z, 
    float *vx, 
    float *vy, 
    float *vz, 
    float *pmass, 
    float *ppot)
{
    int ierr;
    
    if(nobytes==4)
        xvpread_(&model, filename, header,x,y,z,vx,vy,vz,pmass,ppot,&ierr,
            strlen(filename));
    else if(nobytes==8)
        xvprdr8_(&model, filename, header,x,y,z,vx,vy,vz,pmass,ppot,&ierr,
            strlen(filename));
    return ierr;
}
          
/*
 * nbody_guess: guess optimal 'nbody' from the input file
 *              this routine reads the first 4 bytes as a float,
 *              converts it to an integer, and assumes this is
 *              the number of particles in this snapshot series
 */

int nbody_guess(string infile, int nobytes)
{
    stream instr = stropen(infile,"r");
    int nbody;
    float fbuf;
    double dbuf;
    
    if (nobytes==4) {
        if (fread(&fbuf,sizeof(float),1,instr) != 1)
            warning("Problems estimating nbody= from %s",infile);
        nbody = fbuf;
    } else if (nobytes==8) {
        if (fread(&dbuf,sizeof(double),1,instr) != 1)
            warning("Problems estimating nbody= from %s",infile);
        nbody = dbuf;
    } else {
        warning("unsupported nobytes=%d",nobytes);
        nbody=10000;
    }
    strclose(instr);
    return nbody;
}
