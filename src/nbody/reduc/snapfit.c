/*
 * SNAPFIT: various fitting methods to snapshots
 *          
 *  Originally designed to match an observed 3D datacube to 
 *  a 6D model. 
 *      28-may-92   V0.0 Initial design - fairly interactive
 *	 4-may-92    0.4 fixed matrix rotation order bug
 *	16-jun-92    0.6 rotation order -> math positive (see snaprotate(1))
 *      17-jul-03    0.7 make it compile under nemo 3.1.1
 *      29              a  : patch, since scratch files appear to not work
 *       2-aug-03       b  : patch back, stropen() has been fixed
 *	20-may-04       c  : declare local variables for times()
 *      15-jul-04       d  : SINGLE_PREC cleanup
 */

/* this code will only work in 3D */
#define THREEDIM

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#if 0
#include <bodytrans.h>    /* this defines some conflicting macro's, even in nrutil.h */
#else
typedef real (*rproc_body)(Body *, real, int);
extern rproc_body btrtrans(string expr);
#endif
#include <image.h>

typedef struct {    /* special fake cube particle to aid coding */
    vector r;		/* position in the cube */
    real w;		/* weight or intensity */
} cube, *cubeptr;

#define RAD2DEG   (180.0/PI)
#define DEG2RAD   (PI/180.0)

#define MAXANG    181

/* -------------------------------------------------------------------------- */

string defv[] = {
    "in=???\n           Input model snapshot(5NEMO)",
    "cube=???\n         Input data cube (image(5NEMO) or ascii table)",
    "cutoff=\n          Cutoff applied to datacube, if needed.",
    "weight=1\n         Weight applied to model data",
    "out=\n             Optional output",
    "minchi=\n          Min chi-quared to trigger output",
    "frame=model\n      Output in model|data units",
    "zmode=separ\n      zvar needs simultaneous|separate fit",
    "theta1=\n          Inclination Angles to test",
    "theta2=\n          Position Angle Angles to test",
    "contour=\n         Optional Output image file with chi-squared",
    "iter=0\n           number of more iterations after best on matrix",
    "times=all\n        Times selected from snapshot models",
    "VERSION=0.7d\n     15-jul-04 PJT",
    NULL,
};

string usage = "fit snapshots to some model";

/* -------------------------------------------------------------------------- */

    /* Model: Snapshot */
local stream instr, outstr=NULL;   /* file pointer to input and output */
local body *btab=NULL;             /* snapshot */
local cube *mtab=NULL;	     /* to be fitted coordinates out of the snapshot */
local int nmodel=0;		     /* number of model bodies */
local real tsnap;
local string times;
local rproc_body weight;

    /* Data: A special 3D {x,y,v} body */
local cube *dtab=NULL;	/* observed data */
local int ndata=0;		/* number of observed data */

    /* placeholders for array of angles to test */
local real theta1[MAXANG], theta2[MAXANG];
local int ntheta1, ntheta2;

local stream constr=NULL;   /* images of chi-squared country */

local int my_debug=2;
local bool Qsimul;
local real rscale, vscale;
local int maxiter=0;


/* forward declarations */

void setparams(void);
int read_model(void);
void read_data(void);
real *mk_coords(int n, real start, real incr);
void printvec(string name, vector vec);
void snap_fit(void);
int write_snapshot(stream outstr, int nbody, body *btab, real t1, real t2, real rscale, real vscale);
void yrotate(matrix mat, real theta);
void zrotate(matrix mat, real theta);
int eigenframe(vector frame[], matrix mat);
int invert(vector frame[]);

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
    string zmode;

    if (hasvalue("theta1")) {
        ntheta1 = nemoinpr(getparam("theta1"),theta1,MAXANG);
        if (ntheta1<1) error("Too many theta1's: maximum %d",MAXANG);
    }
    if (hasvalue("theta2")) {
        ntheta2 = nemoinpr(getparam("theta2"),theta2,MAXANG);
        if (ntheta2<1) error("Too many theta2's: maximum %d",MAXANG);
    }
    instr=stropen(getparam("in"),"r");
    if(hasvalue("out"))
        outstr=stropen(getparam("out"),"w");
    if(hasvalue("contour"))
        constr=stropen(getparam("contour"),"w");
    zmode = getparam("zmode");
    if (strncmp(zmode,"si",2)==0)
        Qsimul = TRUE;
    else if (strncmp(zmode,"se",2)==0)
        Qsimul = FALSE;
    else
        error("%s: zmode must be 'simultaneous' or 'separate'",zmode);

    weight = btrtrans(getparam("weight"));
    maxiter = getiparam("iter");    
    times = getparam("times");
}

/*
 * read_model: read in the model (snapshot) into the 'btab' structure
 */

int read_model()
{
    int i, bits, count;
    Body *bp, *bq;

    for(;;) {       /* infinite loop until some data found, or nothing */
        get_history(instr);             /* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
            return 0;			/* eof */
        get_snap(instr, &btab, &nmodel, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0)
            continue;       /* skip diagnostics */
        break;
    }
    /* 
     * The first time around, when mtab=NULL, allocate the
     * 'mtab' which can hold the 3 model coordinates to be directly
     * compared to the data 'dtab'
     */

    if (mtab==NULL)
        mtab = (cube *) allocate(nmodel*sizeof(cube));
    
    count=0;
    for (bp=btab, i=0; i<nmodel; bp++, i++) {
        mtab[i].w = (weight)(bp,tsnap,i);
        if (mtab[i].w > 0) count++;
    }
    dprintf(0,"Snapshot: time=%g Using %d/%d bodies\n",tsnap,count,nmodel);

    return 1;
}

/* 
 * read_data: read in data, and stuff them into the 'cube' structure
 *            each pixel get world coordinates from a lookup table 
 */

vector unit_frame[3] = {
    { 1.0, 0.0, 0.0, },
    { 0.0, 1.0, 0.0, },
    { 0.0, 0.0, 1.0, },
};


vector data_frame[3];

void read_data()
{
    imageptr  iptr=NULL;
    stream instr, scrstr;
    bool Qall;
    int nx,ny,nz, ix,iy,iz, i;
    real f, cutoff, sum, vals[4];
    real *xw, *yw, *zw;         /* lookup table */
    vector tmpv, tmpr, data_com, frame[3];
    matrix tmpm, qpole;
    char line[256];

    if (hasvalue("cube")) {
        Qall = !hasvalue("cutoff");
        if (!Qall) cutoff=getdparam("cutoff");
        instr =  stropen(getparam("cube"),"r");

        if (qsf(instr)) {       /* if binary; assume it's an image */
	   dprintf(1,"BINARY image cube\n");
           rewind(instr);
    
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
            * first loop over all points, count how many, 
            */
           for(ix=0; ix<nx; ix++) {
              for(iy=0; iy<ny; iy++) {
                 for(iz=0; iz<nz; iz++) {
                    f=CubeValue(iptr,ix,iy,iz);
                    if(Qall || f>=cutoff) ndata++;
                 }
              }
           }
           dprintf(0,"%d/%d data (%g%%) from the cube used in fit\n",
                   ndata, nx*ny*nz, 100*(real)ndata / (real)(nx*ny*nz));
           if (ndata==0) error("No data found; use different cutoff");


           /* 
            * allocate space, and stuff data into array 
            */
           dtab = (cube *) allocate(ndata*sizeof(cube));
           i = 0;
           for(ix=0; ix<nx; ix++) {
              for(iy=0; iy<ny; iy++) {
                 for(iz=0; iz<nz; iz++) {
                    f=CubeValue(iptr,ix,iy,iz);
                    if(Qall || f>=cutoff) {
                        dtab[i].r[0] = xw[ix];
                        dtab[i].r[1] = yw[iy];
                        dtab[i].r[2] = zw[iz];
                        dtab[i].w    = f;
                        i++;
                    }
                 }
              }
           }
           free_image(iptr);
        } else {                /* assume it's a table */
	   dprintf(1,"ASCII table\n");
           rewind(instr);
           scrstr = stropen("xxx","s");     /* scratch file for 'cube' */
           while ( fgets(line,256,instr) != NULL) {
              if(line[0]=='#') continue;
              line[strlen(line)-1] = '\0';
              switch (nemoinpr(line,vals,4)) {
                case 2:
                  vals[2] = 0.0;      /* allow rare 2D fit?? */
                case 3:
                  vals[3] = 1.0;      /* fill in default weight */
                case 4:
                  break;
                default:
                  warning("Parsing: %s",line);
                  break;
              }
	      if (Qall || vals[3]>cutoff) {
		put_data(scrstr,"XYVI",RealType,vals,4,0);
		ndata++;
	      }
           }
	   rewind(scrstr);
	   dprintf(0,"%d valid data found\n",ndata);
           if (ndata==0) error("No data found; use different cutoff");
         
           dtab = (cube *) allocate(ndata*sizeof(cube));

           for(i=0; i<ndata; i++) {
              get_data(scrstr,"XYVI",RealType,vals,4,0);
	      dprintf(1,"%g %g %g %g\n",vals[0],vals[1],vals[2],vals[3]);
              dtab[i].r[0] = vals[0];
              dtab[i].r[1] = vals[1];
              dtab[i].r[2] = vals[2];
              dtab[i].w    = vals[3];
           }
           strclose(scrstr);
        }
        strclose(instr);
    } else
        error("No cube or table input");

    /* 
     * since the data is never modified throughout the fitting process
     * it's convenient to compute all relevant scaling and rotation 
     * angles which need to be compared with the model later on.
     */

    sum=0;                          /* center of mass */
    CLRV(data_com);
    for(i=0; i<ndata; i++) {
        dprintf(3,"%g %g %g %g\n",
                dtab[i].r[0], dtab[i].r[1], dtab[i].r[2], dtab[i].w);
        sum += dtab[i].w;
	MULVS(tmpv, dtab[i].r, dtab[i].w);
        ADDV(data_com, data_com, tmpv);
    }
    if (sum==0.0) error("Total data weighs 0.0");
    DIVVS(data_com,data_com,sum);
    printvec("Data - C.O.M.  :   ", data_com);

    CLRM(qpole);                    /* moment of inertia */
    for(i=0; i<ndata; i++) {
        SUBV(tmpr, dtab[i].r, data_com);
        MULVS(tmpv, tmpr, dtab[i].w);
        OUTVP(tmpm, tmpv, tmpr);
        ADDM(qpole, qpole, tmpm);
    }
    DIVMS(qpole,qpole,sum);
    printvec("       qpole[0]:   ", qpole[0]);
    printvec("       qpole[1]:   ", qpole[1]);
    printvec("       qpole[2]:   ", qpole[2]);
    rscale = sqrt(qpole[0][0]+qpole[1][1]);
    vscale = sqrt(qpole[2][2]);

    eigenframe(frame, qpole);			/* get rot frame */
    if (dotvp(unit_frame[0], frame[0]) < 0.0)
        MULVS(frame[0], frame[0], -1.0);   
    if (dotvp(unit_frame[2], frame[2]) < 0.0)
        MULVS(frame[2], frame[2], -1.0);   
    CROSSVP(frame[1], frame[2], frame[0]); 	/* make it R.H. */
    printvec(" e_x:", frame[0]);		/* show it */
    printvec(" e_y:", frame[1]);
    printvec(" e_z:", frame[2]);
    for (i = 0; i < NDIM; i++)			/* save it */
        SETV(data_frame[i], frame[i]);

}

/*
 *  mk_coords:  build a lookup table of world coordinates
 */

real *mk_coords(n, start, incr)
int n;
real start, incr;
{
    real *a;
    int i;

    if (n<=0) error("mk_coords: bad array length %d",n);
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
        dprintf(my_debug,"%s  %10.5f  %10.5f  %10.5f  %10.5f\n",
                   name, absv(vec), vec[0], vec[1], vec[2]);
                   
}

/*
 * snap_fit: the actual work horse 
 */

vector model_frame[3];

void snap_fit()
{
    int i, i1, i2, i1_min=-1, i2_min=-1;
    matrix mat1, mat2, rot, tmpm, w_qpole;
    vector tmpv, tmpr, w_pos, frame[3], framet[3];
    real w_sum, sum, sum_min=1000.0, w_rscale=1, w_vscale=1, rscale_min, vscale_min;
    Body tmpp, *bp, *qp=&tmpp;
    imageptr iptr=NULL;


    if (ntheta1==0 || ntheta2==0) {
        error("No searching implemented yet, search must be manual");
        return;
    }
    if (ntheta1==1 || ntheta2==1) my_debug = 1;
    /* case 1: both theta's fixed: for those values fit is shown */
    dprintf(my_debug,"Model fit: \n");

    if (ntheta1>1 && ntheta2>2) {
        create_image(&iptr,ntheta1,ntheta2);
        Xmin(iptr) = theta1[0];     
        Ymin(iptr) = theta2[0];
        Dx(iptr) = theta1[1] - theta1[0];
        Dy(iptr) = theta2[1] - theta2[0];
    }


    printf("2\\1   ");
    for (i1=0; i1<ntheta1; i1++)
        printf(" %5.1f",theta1[i1]);
    printf("\n\n");


    for (i2=0; i2<ntheta2; i2++) {
        zrotate(mat2,theta2[i2]);
        printf("%5.1f :",theta2[i2]);
        for (i1=0; i1<ntheta1; i1++) {          
            dprintf(my_debug,"Theta1,2= %g %g\n", theta1[i1], theta2[i2]);
            yrotate(mat1,theta1[i1]);
            MULM(rot, mat2, mat1);      /* rotation matrix : order=yz */

#if 1
            printvec("  rot[0]: ", rot[0]);
            printvec("  rot[1]: ", rot[1]);
            printvec("  rot[2]: ", rot[2]);
#endif

            w_sum = 0.0;
            CLRV(w_pos);                /* C.O.M. of this model cube */
            for(i=0, bp=btab; i<nmodel; i++, bp++) {
                if (mtab[i].w <= 0) continue;
                MULMV(tmpr, rot, Pos(bp));
                MULMV(tmpv, rot, Vel(bp));
                mtab[i].r[0] = tmpr[0];     /* xvar, yvar, zvar */
                mtab[i].r[1] = tmpr[1];
                mtab[i].r[2] = tmpv[2];
                w_sum += mtab[i].w;
                SETV(tmpv,mtab[i].r);
                MULVS(tmpv,tmpv,mtab[i].w);
                ADDV(w_pos, w_pos, tmpv);
            }
            if(w_sum==0.0) error("weight is zero");
            DIVVS(w_pos,w_pos,w_sum);
            printvec("       C.O.M.  :   ", w_pos);

            CLRM(w_qpole);
            for(i=0, bp=btab; i<nmodel; i++, bp++) {
                if (mtab[i].w <= 0) continue;
                SUBV(tmpr, mtab[i].r, w_pos);
                MULVS(tmpv, tmpr, mtab[i].w);
                OUTVP(tmpm, tmpv, tmpr);
                ADDM(w_qpole, w_qpole, tmpm);
            }
            DIVMS(w_qpole, w_qpole, w_sum);
            if (!Qsimul) {
                w_rscale = sqrt(w_qpole[0][0]+w_qpole[1][1]);
                w_vscale = sqrt(w_qpole[2][2]);
                w_qpole[0][0] *= rscale*rscale/(w_rscale*w_rscale);
                w_qpole[0][1] *= rscale*rscale/(w_rscale*w_rscale);
                w_qpole[1][0] *= rscale*rscale/(w_rscale*w_rscale);
                w_qpole[1][1] *= rscale*rscale/(w_rscale*w_rscale);
                w_qpole[0][2] *= rscale*vscale/(w_rscale*w_vscale);
                w_qpole[1][2] *= rscale*vscale/(w_rscale*w_vscale);
                w_qpole[2][0] *= rscale*vscale/(w_rscale*w_vscale);
                w_qpole[2][1] *= rscale*vscale/(w_rscale*w_vscale);
                w_qpole[2][2] *= vscale*vscale/(w_vscale*w_vscale);
            }
            printvec("       qpole[0]:   ", w_qpole[0]);
            printvec("       qpole[1]:   ", w_qpole[1]);
            printvec("       qpole[2]:   ", w_qpole[2]);

            eigenframe(frame, w_qpole);			/* get rot frame */
            if (dotvp(unit_frame[0], frame[0]) < 0.0)
               MULVS(frame[0], frame[0], -1.0);   
            if (dotvp(unit_frame[2], frame[2]) < 0.0)
               MULVS(frame[2], frame[2], -1.0);   
            CROSSVP(frame[1], frame[2], frame[0]); 	/* make it R.H. */
            for (i = 0; i < NDIM; i++)
                SETV(model_frame[i], frame[i]);

            printvec(" frame e_x   :", frame[0]);
            printvec("       e_y   :", frame[1]);
            printvec("       e_z   :", frame[2]);
            invert(frame);
            TRANM(model_frame, frame);       /* model_frame is now inverse */
                                            /* and can be mult'd with model_ */
                                            /* and compared with unit I */
            printvec(" e_x^-1   :", model_frame[0]);
            printvec(" e_y^-1   :", model_frame[1]);
            printvec(" e_z^-1   :", model_frame[2]);



            sum = 0.0;
            for (i=0; i<3; i++) {
                MULMV(tmpv, data_frame, model_frame[i]);
                printvec(" data*model^-1:", tmpv);
                SUBV(tmpv, tmpv, unit_frame[i]);
                sum += absv(tmpv);
            }
            dprintf(my_debug,"Theta1,2,sum= %g %g %g\n",theta1[i1], theta2[i2], sum);
            sum = log10(sum);
            if (sum<sum_min) {
                sum_min = sum;
                i1_min = i1;
                i2_min = i2;
                rscale_min = w_rscale;
                vscale_min = w_vscale;
            }
            printf(" %5.2f",sum);
            MapValue(iptr,i1,i2) = sum;
        } /* i1 */
        printf("\n");
    } /* i2 */

    printf("\n");
    printf("       ");
    for (i1=0; i1<ntheta1; i1++)
        printf(" %5.1f",theta1[i1]);
    printf("\nMin: theta1= %g theta2= %g log10(sum)= %g time= %g scale= %g %g %g %g\n",
            theta1[i1_min], theta2[i2_min], sum_min, tsnap,
            rscale, vscale, rscale_min, vscale_min);
    if (constr) write_image(constr,iptr);
    if (outstr) write_snapshot( outstr,
				nmodel,btab,
				theta1[i1_min], theta2[i2_min], 
				rscale/w_rscale, vscale/w_vscale);

}

write_snapshot( outstr, nbody, btab, t1, t2, rscale, vscale)
stream outstr;
int nbody;
Body *btab;
real t1, t2, rscale,vscale;
{
    warning("output snapshot not supported yet");
}

/* rotating: is positive when rotating counter-clockwise, or math
 *           positive 
 */

/*
 * yrotate: construct a rotation matrix for a rotation of 'theta' 
 *          around the y-axis
 */

void yrotate(mat,theta)
matrix mat;
real theta;
{
    SETMI(mat);         /* unit matrix */
    mat[2][2] =    mat[0][0] = cos(DEG2RAD * theta);
    mat[2][0] =  -(mat[0][2] = sin(DEG2RAD * theta)); /* PJT */
/*  mat[0][2] =  -(mat[2][0] = sin(DEG2RAD * theta)); /* JEB */

}

/*
 * yrotate: construct a rotation matrix for a rotation of 'theta' 
 *          around the y-axis
 */

void zrotate(mat,theta)
matrix mat;
real theta;
{
    SETMI(mat);         /* unit matrix */
    mat[1][1] =    mat[0][0] = cos(DEG2RAD * theta);
    mat[0][1] =  -(mat[1][0] = sin(DEG2RAD * theta)); /* PJT */
/*  mat[1][0] =  -(mat[0][1] = sin(DEG2RAD * theta)); /* JEB */
}


#include "nrutil.h"

#if 1

eigenframe(frame, mat)          /* float version */
vector frame[];
matrix mat;
{   
    float **q, *d, **v;
    int i, j, nrot;
        
    q = fmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
        for (j = 1; j <= 3; j++)
            q[i][j] = mat[i-1][j-1];
    d = fvector(1, 3);
    v = fmatrix(1, 3, 1, 3);
    jacobi(q, 3, d, v, &nrot);
    eigsrt(d, v, 3);
    for (i = 1; i <= 3; i++)
        for (j = 1; j <= 3; j++)
            frame[i-1][j-1] = v[j][i];
} 

#else

eigenframe(frame, mat)          /* double version */
vector frame[];
matrix mat;
{   
    double **q, *d, **v;
    int i, j, nrot;
        
    q = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
        for (j = 1; j <= 3; j++)
            q[i][j] = mat[i-1][j-1];
    d = dvector(1, 3);
    v = dmatrix(1, 3, 1, 3);
    jacobi_d(q, 3, d, v, &nrot);
    eigsrt_d(d, v, 3);
    for (i = 1; i <= 3; i++)
        for (j = 1; j <= 3; j++)
            frame[i-1][j-1] = v[j][i];
} 
#endif 

invert(frame)
vector frame[];
{
    real mat[NDIM*NDIM], det;

    matinv(frame, 3, 3, &det);
}
