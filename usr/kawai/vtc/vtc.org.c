#pragma global noalias
#define CALCPOTENTIAL (1)

int md2_mpi_rank; /* dummy */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "vtc.h"
#include "vtclocal.h"

/* parameter for file I/O etc */
static int snapin_flag = 0;
static int snapout_flag = 0;
static char snapinfile[255]={'?'};
static char snapoutfile[255]={'?'};
static int snapout_acc_flag = 0;
static int snapio_binary_flag = 0;
static int ionemo_flag = 0;
static int command_exec_flag = 0;
static char command_name[255]={'?'};
static int outlogstep = 1;
static int cputime_flag = 0;
static int test_id = -1;
#ifdef ANIM
static int animation_flag = 0;
#endif /* ANIM */
static int initcond_flag = 0;
static int initn;

/* parameter for time integration */
static double  dt = 0.01;
static double  dtsnapout = 1;
static double tstart = 0.0;
static double tend = 10.0;
static double  dtimageout = 1;
static int imageout_flag = 0;
static double image_scale = 0.3;

/* parameter for force calculation */
static double mmin = 0.0;
static double xmin = 0.0, xmax = 0.0;
static int autoscale = 1;
static int node_div_crit = 8;
static double eps = 0.02;
static double theta = 0.75;
static int ncrit = 8192;
static int grape_flag = 0;
static double pprad = 0.2;
static int full_dof_flag = TRUE;
static int me_order = 1;
static int negativemass = FALSE;

static void
show_usage(char *command)
{
    fprintf(stderr, "usage: %s -i fname [options]\n", command);
#ifdef ANIM
    fprintf(stderr, "    -A            animate\n");
#endif /* ANIM */
    fprintf(stderr, "    -a            output acceleration in addition to m, x, v, p\n");
#ifdef IONEMO
    fprintf(stderr, "    -b            use NEMO binary format I/O (use ASCII otherwise)\n");
#endif /* IONEMO */
    fprintf(stderr, "    -c            output timing info\n");
    fprintf(stderr, "    -D num        snapshot output interval [1]\n");
    fprintf(stderr, "    -d num        timestep [0.01]\n");
    fprintf(stderr, "    -e num        softening [0.02]\n");
    fprintf(stderr, "    -F            never use 'modified' P2M2 (i.e. use regular P2M2 for p=2)\n");
#ifndef NOGRAPE
    fprintf(stderr, "    -G            use GRAPE if available\n");
#endif /* NOGRAPE */
    fprintf(stderr, "    -i fname      input file [?]\n");
    fprintf(stderr, "    -I num        use num-particle homo sphere as initial distribution\n");
    fprintf(stderr, "    -l num        log output interval [1]\n");
    fprintf(stderr, "    -N num        max num of particles each leaf cell can have [8]\n");
    fprintf(stderr, "    -n num        Barnes' vectorization parameter [8192]\n");
    fprintf(stderr, "    -o fname      output file in printf format [none] ex) -o out%%03d\n");
    fprintf(stderr, "    -p num        multipole-expansion order [1] (0<p<7)\n");
    fprintf(stderr, "    -r num        sphere radius on which pseudoparticles are distributed [0.2]\n");
    fprintf(stderr, "    -s num        scaling of image to animate [0.3]\n");
    fprintf(stderr, "    -T num        end time [10.0]\n");
    fprintf(stderr, "    -t num        opening angle [0.75]\n");
#if USEX11
    fprintf(stderr, "    -V num        test mode [-1] (try 0, 1, 2 and see what will happen\n");
    fprintf(stderr, "                              control the view with keys +-hjklq spc)\n");
#endif /* USEX11 */
    fprintf(stderr, "    -X command    execute command immediately after snapshot output\n");
    fprintf(stderr, "    -W num        GIF image output interval [none]\n");
    fprintf(stderr, "    -z            mass can be negative (for MD simulation)\n");
    fprintf(stderr, "    -h            print this help\n");
    exit(1);
}

#ifdef NOGRAPE
#define GRAPEOPTS ""
#else
#define GRAPEOPTS "G"
#endif /* NOGRAPE */

#ifdef IONEMO
#define IONEMOOPTS "b"
#else
#define IONEMOOPTS ""
#endif /* IONEMO */

#ifdef ANIM
#define ANIMOPTS "A"
#else
#define ANIMOPTS ""
#endif /* ANIM */

#if USEX11
#define TESTOPTS "V:"
#else
#define TESTOPTS ""
#endif /* USEX11 */

extern char *optarg;
extern int optind;
static void
parse_argv(int argc, char **argv)
{
    int c;
    char *infile = NULL;
    char* param = "acD:d:e:FI:i:l:N:n:o:p:r:s:T:t:W:X:zh" GRAPEOPTS IONEMOOPTS ANIMOPTS TESTOPTS;

    while ((c = getopt(argc, argv, param)) != EOF) {
	switch (c) {
#ifdef ANIM
	case 'A':
            animation_flag = 1;
	    break;
#endif /* ANIM */
	case 'a':
	    snapout_acc_flag = 1;
	    break;
#ifdef IONEMO
	case 'b':
	    ionemo_flag = 1;
	    break;
#endif /* IONEMO */
	case 'c':
	    cputime_flag = 1;
	    break;
	case 'D':
	    dtsnapout = atof(optarg);
	    break;
	case 'd':
	    dt = atof(optarg);
	    break;
	case 'e':
	    eps = atof(optarg);
	    break;
	case 'F':
	    full_dof_flag = FALSE;
	    break;
#ifndef NOGRAPE
	case 'G':
	    grape_flag = 1;
	    break;
#endif /* NOGRAPE */
	case 'i':
	    strcpy(snapinfile, optarg);
	    snapin_flag = 1;
	    break;
	case 'I':
	    initcond_flag = 1;
	    initn = atoi(optarg);
	    break;
	case 'l':
	    outlogstep = atoi(optarg);
	    break;
	case 'N':
	    node_div_crit = atoi(optarg);
	    if (node_div_crit < 1) {
		fprintf(stderr, "-N num must be positive\n");
		exit(1);
	    }
	    break;
	case 'n':
	    ncrit = atof(optarg);
	    if (ncrit < 1) {
		fprintf(stderr, "-n num must be positive\n");
		exit(1);
	    }
	    break;
	case 'o':
	    strcpy(snapoutfile, optarg);
	    snapout_flag = 1;
	    break;
	case 'p':
	    me_order = atoi(optarg);
	    break;
	case 'r':
	    pprad = atof(optarg);
	    break;
	case 's':
	    image_scale = atof(optarg);
	    break;
	case 'T':
	    tend = atof(optarg);
	    break;
	case 't':
	    theta = atof(optarg);
	    break;
#if USEX11
	case 'V':
	    test_id = atoi(optarg);
	    break;
#endif /* USEX11 */
	case 'X':
	    strcpy(command_name, optarg);
	    command_exec_flag = 1;
	    break;
	case 'W':
	    dtimageout = atof(optarg);
	    imageout_flag = 1;
	    break;
	case 'z':
	    negativemass = TRUE;
	    break;
	case 'h':
	default:
	    show_usage(argv[0]);
	    break;
	}
    }
}

static void
write_ascii(char *fname, double time, Nbodyinfo *nb)
{
    FILE *fp;
    int i, k;
    int n = nb->n;

    fprintf(stderr, "will write %s...", fname);
    fp = fopen(fname, "w");
    if (fp == NULL) {
	perror("weite_ascii");
	exit(1);
    }
    fprintf(fp, "%d\n", n);
    fprintf(fp, "%d\n", 3);
    fprintf(fp, "% 15.13E\n", time);
    for (i = 0; i < n; i++) {
	fprintf(fp, " % 15.13E\n", nb->m[i]);
    }
    for (i = 0; i < n; i++) {
	fprintf(fp, " % 15.13E % 15.13E % 15.13E\n",
		nb->x[i][0], nb->x[i][1], nb->x[i][2]);
    }
    for (i = 0; i < n; i++) {
	fprintf(fp, " % 15.13E % 15.13E % 15.13E\n",
		nb->v[i][0], nb->v[i][1], nb->v[i][2]);
    }

    if (snapout_acc_flag == 1) {
	for (i = 0; i < n; i++) {
	    fprintf(fp, " % 15.13E\n", nb->p[i]);
	}
	for (i = 0; i < n; i++) {
	    fprintf(fp, " % 15.13E % 15.13E % 15.13E\n",
		    nb->a[i][0], nb->a[i][1], nb->a[i][2]);
	}
    }
    fprintf(stderr, "done\n");
    fclose(fp);
}

static void
create_initcond(Nbodyinfo *nb, int n)
{
    int i, k;
    double m = 1.0/n;

    nb->n = n;
    nb->m = (double *)malloc(sizeof(double)*n);
    if (NULL == nb->m) {
	perror("create_initcond");
	exit(1);
    }
    nb->x = (double (*)[3])malloc(sizeof(double)*3*n);
    nb->v = (double (*)[3])malloc(sizeof(double)*3*n);
    if (NULL == nb->x || NULL == nb->v) {
	perror("create_initcond");
	exit(1);
    }
    for (i = 0; i < n; i++) {
	nb->m[i] = m;
    }
    for (i = 0; i < n; i++) {
	nb->v[i][0] = 0.0;
	nb->v[i][1] = 0.0;
	nb->v[i][2] = 0.0;
    }
    i = 0;
    while (i < n) {
	double x, y, z;
	x = 2.0*drand48()-1.0;
	y = 2.0*drand48()-1.0;
	z = 2.0*drand48()-1.0;
	if (x*x+y*y+z*z<1.0) {
	    nb->x[i][0] = x;
	    nb->x[i][1] = y;
	    nb->x[i][2] = z;
	    i++;
	}
    }
}

static void
read_binary(char *fname, double *timep, Nbodyinfo *nb)
{
#ifdef IONEMO
    int *np = NULL;
    double *tp = NULL;

    io_nemo(fname, "double, n, t, x, v, m, read, info",
	    &np, &tp, &(nb->x), &(nb->v), &(nb->m));
    io_nemo(fname,"close");
    nb->n = *np;
    if (tp == NULL) {
	fprintf(stderr, "read_binary: no time tag in input file %s\n", fname);
	fprintf(stderr, "set time to 0.0\n");
	*timep = 0.0;
    }
    else {
	*timep = *tp;
    }
    fprintf(stderr, "done.\n");

#else /* IONEMO */
    fprintf(stderr, "binary I/O is not supported\n");
    exit(1);
#endif /* IONEMO */
}

static void
write_binary(char *fname, double time, Nbodyinfo *nb)
{
#ifdef IONEMO
    static int *np = NULL;
    static double *tp = NULL;

    if (!np) {
	np = (int *)malloc(sizeof(int));
	tp = (double *)malloc(sizeof(double));
    }
    np[0] = nb->n;
    tp[0] = time;
    if (snapout_acc_flag == 1) {
	io_nemo(fname, "double, n, t, x, v, m, p, a, save, info",
		&np, &tp, &(nb->x), &(nb->v), &(nb->m), &(nb->p), &(nb->a));
    }
    else {
	io_nemo(fname, "double, n, t, x, v, m, p, save, info",
		&np, &tp, &(nb->x), &(nb->v), &(nb->m), &(nb->p));
    }
    io_nemo(fname,"close");

#else /* IONEMO */
    fprintf(stderr, "binary I/O is not supported\n");
    exit(1);
#endif /* IONEMO */
}

static void
read_ascii(char *fname, double *timep, Nbodyinfo *nb)
{
    FILE *fp;
    int i, k;
    int n;
    double dummy;

    fprintf(stderr, "will read %s...", fname);
    fp = fopen(fname, "r");
    if (fp == NULL) {
	perror("read_ascii");
	exit(1);
    }
    fscanf(fp, "%d\n", &(nb->n));
    fscanf(fp, "%d\n", &dummy);
    fscanf(fp, "%lf\n", timep);
    n = nb->n;
    nb->m = (double *)malloc(sizeof(double)*n);
    if (NULL == nb->m) {
	perror("read_ascii");
	exit(1);
    }
    nb->x = (double (*)[3])malloc(sizeof(double)*3*n);
    nb->v = (double (*)[3])malloc(sizeof(double)*3*n);
    if (NULL == nb->m || NULL == nb->x || NULL == nb->v) {
	perror("read_ascii");
	exit(1);
    }
    for (i = 0; i < n; i++) {
	fscanf(fp, "%lf\n", nb->m+i);
    }
    for (i = 0; i < n; i++) {
	fscanf(fp, "%lf %lf %lf\n",
	       nb->x[i]+0, nb->x[i]+1, nb->x[i]+2);
    }
    for (i = 0; i < n; i++) {
	fscanf(fp, "%lf %lf %lf\n",
	       nb->v[i]+0, nb->v[i]+1, nb->v[i]+2);
    }
    fclose(fp);
    fprintf(stderr, "done.\n");
}

static double
get_potential_energy(Nbodyinfo *nb)
{
    int i;
    double e = 0.0;

    for (i = 0; i < nb->n; i++) {
	e += nb->m[i]*nb->p[i];
	Cfprintf(stderr, "p[%d] = %6.5f\n", i, nb->p[i]);
    }
    e *= 0.5;
    return e;
}

static double
get_kinetic_energy(Nbodyinfo *nb)
{
    int i, k;
    double e = 0.0;
    double v2;

    for (i = 0; i < nb->n; i++) {
	v2 = 0.0;
	for (k = 0; k < 3; k++) {
	    v2 += nb->v[i][k]*nb->v[i][k];
	}
	e += nb->m[i]*v2;
    }
    e *= 0.5;
    return e;
}

void
plot_nbody(double time, Nbodyinfo *nb)
{
#ifdef ANIM
    if (animation_flag) {
	int i;
	int err = 0;
	plot_star(nb->x, nb->n, time, image_scale, nb->m, nb->m[0]);
    }
#endif /* ANIM */
}

static void
push_velocity(double dt, Nbodyinfo *nb)
{
    int i;
    int n = nb->n;
    double *nbv0, *nbv1, *nbv2;
    double *nba0, *nba1, *nba2;

    nbv0 = (double *)&(nb->v[0][0]);
    nbv1 = (double *)&(nb->v[0][1]);
    nbv2 = (double *)&(nb->v[0][2]);
    nba0 = (double *)&(nb->a[0][0]);
    nba1 = (double *)&(nb->a[0][1]);
    nba2 = (double *)&(nb->a[0][2]);
    for (i = 0; i < n; i++) {
	nbv0[i*3] += dt*nba0[i*3];
	nbv1[i*3] += dt*nba1[i*3];
	nbv2[i*3] += dt*nba2[i*3];
    }
}

static void
push_position(double dt, Nbodyinfo *nb)
{
    int i;
    int n = nb->n;
    double *nbx0, *nbx1, *nbx2;
    double *nbv0, *nbv1, *nbv2;

    nbx0 = (double *)&(nb->x[0][0]);
    nbx1 = (double *)&(nb->x[0][1]);
    nbx2 = (double *)&(nb->x[0][2]);
    nbv0 = (double *)&(nb->v[0][0]);
    nbv1 = (double *)&(nb->v[0][1]);
    nbv2 = (double *)&(nb->v[0][2]);
    for (i = 0; i < n; i++) {
	nbx0[i*3] += dt*nbv0[i*3];
	nbx1[i*3] += dt*nbv1[i*3];
	nbx2[i*3] += dt*nbv2[i*3];
    }
}
void
get_cmterms(Nbodyinfo *nb, double cm[3], double cmv[3])
{
    int i, k;
    int n = nb->n;
    double totm = 0.0;

    for (k = 0; k < 3; k++) {
	cm[k] = 0.0;
	cmv[k] = 0.0;
    }
    for (k = 0; k < 3; k++) {
	for (i = 0; i < n; i++) {
	    totm += nb->m[i];
	    cm[k] += nb->x[i][k]*nb->m[i];
	    cmv[k] += nb->v[i][k]*nb->m[i];
	}
    }
    for (k = 0; k < 3; k++) {
	cm[k] /= totm;
	cmv[k] /= totm;
    }
}

static void
time_integration_loop(Forceinfo *fi, Nbodyinfo *nb)
{
    int i, nstep;
    int outsnapstep, outimagestep;
    double W, K, E0, E, Q;
    double cm[3], cmv[3];
    double t = tstart;

    nstep = (int)((tend-tstart)/dt+0.1);
    dt = (tend-tstart)/nstep;
    outsnapstep = (int) (dtsnapout/dt+0.1);
    outimagestep = (int) (dtimageout/dt+0.1);
    PR(tstart,f); PR(tend,f); PR(dt,f); PRL(dtsnapout,f); PRL(dtimageout,f);

    if (cputime_flag) {
	vtc_turnon_cputime();
    }
    vtc_init_cputime();
    vtc_print_cputime("initial vtc_get_force start at");
#if 0 /* direct sum */
    vtc_get_default_direct_params(fi);
    if (grape_flag) {
	fi->calculator = GRAPE;
	vtc_set_scale(32.0, nb->m[0]);
    }
    vtc_get_force_direct(fi, nb);
#else /* tree */
    if (grape_flag) {
	fi->calculator = GRAPE_FORCEONLY;
	vtc_get_force_tree(fi, nb);
#if CALCPOTENTIAL
        fi->calculator = GRAPE_POTENTIALONLY;
        vtc_get_force_tree(fi, nb);
	fi->calculator = GRAPE_FORCEONLY;
#endif
    }
    else {
	vtc_get_force_tree(fi, nb);
    }
#endif
    vtc_print_cputime("initial vtc_get_force end at");
    fprintf(stderr, "ninteraction: %e list_len_avg: %6.5f nwalk: %d\n",
	    fi->ninteraction, (double)fi->ninteraction/nb->n, fi->nwalk);
    K = get_kinetic_energy(nb);
#if CALCPOTENTIAL
    W = get_potential_energy(nb);
    E = E0 = W+K;
    PR(E0, f); PR(W, f);
#endif
    PRL(K, f);
    plot_nbody(t, nb);

    /* calculate force at T=0 and quit */
    if (snapout_flag && dtsnapout == 0.0) {
	static int nout = 0;
	char fn[255];
	snapout_acc_flag = 1;
	sprintf(fn, snapoutfile, nout);
	PRL(fn, s);
	if (ionemo_flag) {
	    write_binary(fn, t, nb);
	}
	else {
	    write_ascii(fn, t, nb);
	}
	vtc_print_cputime("exit at");
	exit(0);
    }
    for (i = 0; i < nstep; i++) {
 	if (i%outlogstep == outlogstep - 1 && cputime_flag) {
	    vtc_turnon_cputime();
	}
	else {
	    vtc_turnoff_cputime();
	}

        push_velocity(0.5*dt, nb);
        push_position(dt, nb);
	t = tstart+(i+1)*dt;
	vtc_init_cputime();
#if 0 /* direct sum */
	if (grape_flag) {
	    fi->calculator = GRAPE;
	}
	vtc_get_force_direct(fi, nb);
#else /* tree */
	vtc_get_force_tree(fi, nb);
#endif
	vtc_print_cputime("vtc_get_force end at");
        push_velocity(0.5*dt, nb);
	vtc_print_cputime("integration end at");

	plot_nbody(t, nb);
 	if (i%outlogstep == outlogstep - 1) {
	    fprintf(stderr, "ninteraction: %e list_len_avg: %6.5f nwalk: %d\n",
		    fi->ninteraction, fi->ninteraction/nb->n, fi->nwalk);
	    K = get_kinetic_energy(nb);
#if CALCPOTENTIAL
	    if (grape_flag) {
		fi->calculator = GRAPE_POTENTIALONLY;
		vtc_get_force_tree(fi, nb);
		fi->calculator = GRAPE_FORCEONLY;
	    }
	    W = get_potential_energy(nb);
	    E = W+K;
	    Q = K/(K-E);
	    fprintf(stderr, "\nT= %f K= %f Q= %f E= %f Eerr= %f Eabserr= %f\n",
		    t, K, Q, E, (E0-E)/E0, E0-E);
#else
	    fprintf(stderr, "\nT= %f K= %f\n", t, K);
#endif
	    fprintf(stderr, "step: %d\n", i);
	    get_cmterms(nb, cm, cmv);
	    fprintf(stderr, "CM  %15.13e %15.13e %15.13e\n",
		    cm[0], cm[1], cm[2]);
	    fprintf(stderr, "CMV %15.13e %15.13e %15.13e\n",
		    cmv[0], cmv[1], cmv[2]);
 	}
	if (i%outsnapstep == outsnapstep - 1 && snapout_flag) {
	    static int nout = 0;
	    char fn[255];
	    char cmd[255];

	    sprintf(fn, snapoutfile, nout);
	    PRL(fn, s);
	    if (ionemo_flag) {
		write_binary(fn, t, nb);
	    }
	    else {
		write_ascii(fn, t, nb);
	    }

	    if (command_exec_flag) {
		sprintf(cmd, "%s %s %03d", command_name, fn, nout);
		fprintf(stderr, "execute %s\n", cmd);
		system(cmd);
	    }
	    nout++;
	    PR(t, f); PRL(fn, s);
	}
	if (i%outimagestep == outimagestep - 1 && imageout_flag) {
	    static int nout = 0;
	    char fn[255];
	    char cmd[255];

	    sprintf(fn, snapoutfile, nout);
	    strcat(fn, ".gif");
	    PRL(fn, s);
	    vtc_plotstar(nb->x, nb->n, t, image_scale, fn);
	    nout++;
	    PR(t, f); PRL(fn, s);
	}
	vtc_print_cputime("exit at");
    } /* i loop */
}

static void
init_nbodyinfo(Nbodyinfo *nb)
{
    int k;

    nb->m = NULL;
    nb->x = NULL;
    nb->v = NULL;
    nb->a = NULL;
    nb->p = NULL;
}

int
main(int argc, char**argv)
{
    int i, k;
    int n;
    Forceinfo force;
    Nbodyinfo nbody;

    init_nbodyinfo(&nbody);

    vtc_init_cputime();
    parse_argv(argc, argv);
    if (snapin_flag == 0 && !initcond_flag) {
	fprintf(stderr, "snapshot input file required (-i)\n");
	show_usage(argv[0]);
    }
    if (initcond_flag) {
	create_initcond(&nbody, initn);
	tstart = 0.0;
    }
    else {
	if (ionemo_flag) {
	    read_binary(snapinfile, &tstart, &nbody);
	}
	else {
	    read_ascii(snapinfile, &tstart, &nbody);
	}
    }
    n = nbody.n;
    nbody.a = (double (*)[3])malloc(sizeof(double)*3*n);
    nbody.p = (double *)malloc(sizeof(double)*n);
    if (NULL == nbody.a || NULL == nbody.p) {
	perror("main");
	exit(1);
    }
    if (tstart > tend) {
	fprintf(stderr, "start time %f larger than end time %f. abort.\n",
		tstart, tend);
	exit(1);
    }
    vtc_get_default_tree_params(&force); /* get current values */
    /* modify some of them */
    force.theta = theta;
    force.eps = eps;
    force.ncrit = ncrit;
    force.node_div_crit = node_div_crit;
    force.p = me_order;
    force.negativemass = negativemass;
    force.pprad = pprad;
    force.full_dof = full_dof_flag;
    if (grape_flag) {
	force.calculator = GRAPE_FORCEONLY;
	force.calculator = GRAPE;
    }
    force.test_id = test_id;

    PR(snapinfile,s); PRL(snapoutfile,s);
    PRL(n,d); 
    PR(eps, f); PR(theta, f); PR(ncrit, d); PRL(node_div_crit, d); 
    time_integration_loop(&force, &nbody);
#if USE_GM_API
    fprintf(stderr, "will m2_gm_finalize...\n");
    m2_gm_finalize();
    fprintf(stderr, "done m2_gm_finalize.\n");
#endif /* USE_GM_API */
    exit(0);
}
