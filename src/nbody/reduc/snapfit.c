/*
 * SNAPFIT: various fitting methods to snapshots
 *          
 *  Originally designed to match an observed 3D datacube to 
 *  a 6D model. 
 *      28-may-92   V0.0 Initial design - fairly interactive
 *
 */
#define THREEDIM

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <image.h>

typedef struct {    /* special fake cube particle to aid coding */
    vector r;
} cube, *cubeptr;

#define RAD2DEG   (180.0/PI)
#define DEG2RAD   (PI/180.0)

#define MAXANG    181

/* -------------------------------------------------------------------------- */

string defv[] = {
    "in=???\n         Input snapshot with positions",
    "cube=\n            Input cube, with X,Y,Z positions",
    "cutoff=\n          Cutoff applied to datacube, if needed.",
    "xvar=x\n           X-Observable to match",
    "yvar=y\n           Y-Observable to match",
    "zvar=vz\n          Z-Observable to match",
    "weight=1\n         Weight",
    "out=\n             Optional output",
    "minchi=\n          Min chi-quared to trigger output",
    "frame=model\n      Output in model|data units",
    "zmode=simult\n     zvar needs simultaneous|separate fit",
    "theta1=\n          Inclination Angles to test",
    "theta2=\n          Position Angle Angles to test",
    "maxreport=\n       Max. best chi-squares reported if fixed thetas",
    "VERSION=0.1\n      29-may-92 PJT",
    NULL,
};

string usage = "fit snapshots to some model";

/* -------------------------------------------------------------------------- */

    /* Snapshot */
stream instr;
Body *btab=NULL;
cube *mtab=NULL;
int nbody=0;

    /* Special 3D {x,y,v} body to represent the image */
cube *ctab=NULL;
int ndata=0;

    /* placeholders for array of angles to test */
real theta1[MAXANG], theta2[MAXANG];
int ntheta1, ntheta2;
    

   /* local declarations */
void setparams(), read_data(), snap_fit();
void printvec(), yrotate(), zrotate();
real *mk_coords();
int read_model();


/* -------------------------------------------------------------------------- */

void nemo_main()
{ 
    setparams();
    read_data();              /* get the data  (an image, or table ...) */
    while(read_model()) {     /* while models (snapshot) available */
        snap_fit();
    }
}

void setparams()
{
    ntheta1 = nemoinpd(getparam("theta1"),theta1,MAXANG);
    ntheta2 = nemoinpd(getparam("theta2"),theta2,MAXANG);
    instr=stropen(getparam("in"),"r");
}

/*
 * read_model: read in the model (snapshot) into the 'btab' structure
 */

int read_model()
{
    real tsnap, x,y,z;
    int bits;
    Body *bp;
    byte *allocate();

    for(;;) {       /* infinite loop until some data found, or nothing */
        get_history(instr);             /* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
            return 0;
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0)
            continue;       /* skip diagnostics */
        break;
    }
    /* 
     * before we return sucessfully though, stuff the particle
     * weight into the Aux() array, 
     * The first time around, when mtab=NULL, also allocate the
     * mtab which will hold the 3 model coordinates to be directly
     * compared to the data (ctab)
     */

    for (bp=btab; bp<btab+nbody; bp++) {
        Aux(bp) = 1.0;          /* weight */
    }

    if (mtab==NULL) {
        mtab = (cube *) allocate(nbody*sizeof(cube));
    }    
    
    return 1;
}

/* 
 * read_data: read in data, and stuff them into the 'cube' structure
 *            each pixel get world coordinates from a lookup table 
 */

void read_data()
{
    imageptr  iptr=NULL;
    stream instr;
    bool Qall;
    int nx,ny,nz, ix,iy,iz, idx;
    real f, cutoff, sum;
    real *xw, *yw, *zw;         /* lookup table */
    vector tmpv, pos, data_com;
    matrix tmpm, qpole;
    byte *allocate();

    if (hasvalue("cube")) {
         Qall = !hasvalue("cutoff");
        if (!Qall) cutoff=getdparam("cutoff");
        instr =  stropen(getparam("cube"),"r");
        read_image(instr,&iptr);
	if (iptr==NULL) error("%s: Bad image cube",getparam("cube"));
        xw = mk_coords(Nx(iptr), Xmin(iptr), Dx(iptr));
        yw = mk_coords(Ny(iptr), Ymin(iptr), Dy(iptr));
        zw = mk_coords(Nz(iptr), Zmin(iptr), Dz(iptr));
        
        nx = Nx(iptr);                              /* cube dimensions */
        ny = Ny(iptr);
        nz = Nz(iptr);
	if (nz==1) warning("Input data %s has only one plane",
				getparam("cube"));

        /* 
         * first loop over all points, count how many, and get center
         * of mass 
         */
        CLRV(data_com);
        sum=0;
        for(ix=0; ix<nx; ix++) {
            for(iy=0; iy<ny; iy++) {
                for(iz=0; iz<nz; iz++) {
                    f=CubeValue(iptr,ix,iy,iz);
                    if(Qall || f>=cutoff) {
                        sum += f;
                        tmpv[0] = xw[ix];
                        tmpv[1] = xw[iy];
                        tmpv[2] = xw[iz];
                        ADDV(data_com, data_com, tmpv);
                        ndata++;
                    }
                }
            }
        }
        dprintf(0,"%d/%d data (%g%%) from the cube used in fit\n",
                ndata, nx*ny*nz, 100*(real)ndata / (real)(nx*ny*nz));
        if (sum==0.0 || ndata==0) error("No data found");
        DIVVS(data_com,data_com,sum);
        printvec("Data -- Center of Mass:",data_com);

        ctab = (cube *) allocate(ndata*sizeof(cube));
        idx = 0;
        CLRM(qpole);
        for(ix=0; ix<nx; ix++) {
            for(iy=0; iy<ny; iy++) {
                for(iz=0; iz<nz; iz++) {
                    f=CubeValue(iptr,ix,iy,iz);
                    if(Qall || f>=cutoff) {
                        ctab[idx].r[0] = xw[ix];
                        ctab[idx].r[1] = yw[iy];
                        ctab[idx].r[2] = zw[iz];
                        SUBV(pos,ctab[idx].r,data_com);
                        MULVS(tmpv,pos,f);
                        OUTVP(tmpm,tmpv,pos);
                        ADDM(qpole,qpole,tmpm);
                        idx++;
                    }
                }
            }
        }
        DIVMS(qpole,qpole,sum);
        printvec("       qpole[0]:   ", qpole[0]);
        printvec("       qpole[1]:   ", qpole[1]);
        printvec("       qpole[2]:   ", qpole[2]);
        free_image(iptr);
        strclose(instr);
    } else {
        error("Currently no other options then to supply an image; cube=");
    }

    
}

/*
 *  mk_coords:  build a lookup table of world coordinates
 */

real *mk_coords(n, start, incr)
int n;
real start, incr;
{
    byte *allocate();
    real *a;
    int i;

    if (n<=0) error("mk_coords: bad array length %d\n",n);
    a = (real *) allocate(n*sizeof(real));
#if 0    
    for (i=1, a[0]=start; i<n; i++)
        a[i] = a[i-1] + incr;
#else
    for (i=0; i<n; i++)
        a[i] = (real) i;
#endif
    return a;
}

void printvec(name, vec)
string name;
vector vec;
{
        printf("%s  %10.5f  %10.5f  %10.5f  %10.5f\n",
                   name, absv(vec), vec[0], vec[1], vec[2]);
                   
}

/*
 * snap_fit: the actual work horse 
 */

void snap_fit()
{
    int i, i1, i2;
    matrix mat1, mat2, rot;
    vector tmpv;
    Body tmpp, *bp, *qp=&tmpp;


    /* case: both theta's fixed */
    printf("Model fit: \n");

    for (i1=0; i1<ntheta1; i1++) {          
        printf("Theta1=%g\n", theta1[i1]);
        yrotate(mat1,theta1[i1]);
        for (i2=0; i2<ntheta2; i2++) {
            printf("Theta2=%g\n", theta2[i2]);
            zrotate(mat2,theta2[i2]);
            MULM(rot, mat1, mat2);      /* rotation matrix */
            for(i=0, bp=btab; i<nbody; i++, bp++) {
                tmpv[0] = Pos(bp)[0];    /* xvar */
                tmpv[1] = Pos(bp)[1];    /* yvar */
                tmpv[2] = Vel(bp)[2];    /* zvar */
                MULMV(ctab[i].r, rot, tmpv);
                /* accumulate to find center of mass */
                /* ... */
            }

            for(i=0, bp=btab; i<nbody; i++, bp++) {
                /* moment of inertia */
                /* ... */
            }
            
            /* compare moment of inertia of model with that of data */
        }
    }
}

void yrotate(mat,theta)
matrix mat;
real theta;
{
    SETMI(mat);         /* unit matrix */
    mat[2][2] =    mat[0][0] = cos(DEG2RAD * theta);
    mat[0][2] =  -(mat[2][0] = sin(DEG2RAD * theta));
}

void zrotate(mat,theta)
matrix mat;
real theta;
{
    SETMI(mat);         /* unit matrix */
    mat[1][1] =    mat[0][0] = cos(DEG2RAD * theta);
    mat[1][0] =  -(mat[0][1] = sin(DEG2RAD * theta));
}

