#ifndef __VTC_H__
#define __VTC_H__

typedef enum {
  NOCALCULATOR = 0,
  HOST,
  HOST_FORCEONLY,
  HOST_POTENTIALONLY,
  GRAPE,
  GRAPE_FORCEONLY,
  GRAPE_POTENTIALONLY,
} Calculator;

typedef struct {
  int n;
  double *m;
  double (*x)[3];
  double (*v)[3];
  double (*a)[3];
  double *p;
} Nbodyinfo;

typedef struct {
    double theta;          /* opening parameter */
    double eps;            /* softening parameter */
    int ncrit;             /* Barnes' vectorization parameter */
    int node_div_crit;     /* max # of particles a leaf cell can contain */
    int p;                 /* multipole-expansion order */
    int negativemass;      /* if TRUE, particles may have negative mass */
    double pprad;          /* sphere radius on which pseudo particles are distributed */
    int npp;               /* # of pseudo particles */
    int tdesign;               /* spherical t-design */
    double (*pppos)[3];    /* position of pseudo particles */
    int full_dof;          /* use modified P2M2 if p<=2 */
    Calculator calculator; /* HOST, GRAPE, GRAPE_FORCEONLT, etc. */
    double ninteraction;   /* # of paiwise interactions */
    int nwalk;             /* # of tree traversals */
    int test_id;           /* -1 for normal operation */
} Forceinfo;

/* for direct-summation algorithm */
void vtc_get_default_direct_params(Forceinfo *fi);
void vtc_get_force_direct(Forceinfo *tr, Nbodyinfo *nb);

/* for Barnes-Hut tree algorithm */
void vtc_get_default_tree_params(Forceinfo *fi);
void vtc_get_force_tree(Forceinfo *tr, Nbodyinfo *nb);
void vtc_set_scale(double xmax, double mmin); /* fixed-point format hardware specific */
int vtc_set_rscale(double rscale); /* not used for normal operation */
void vtc_set_mac_margin(double margin);
double vtc_get_mac_margin(void);

/* timing analysis tools */
void vtc_init_cputime(void);
double vtc_get_cputime(void);
void vtc_turnon_cputime(void);
void vtc_turnoff_cputime(void);
void vtc_print_cputime(char *comment);

/* visualization tools */
void vtc_viewtree(int n, double (*bpos)[3],
		  int nc, double (*cpos)[3], double *csize);
void vtc_viewlist(int nall, double (*ballpos)[3], int n, double (*bpos)[3],
		  int nc, double (*cpos)[3], double *csize, int ni);
void vtc_plotstar(char *fname, Nbodyinfo *nb, char *msg, double scale, double center[2], double cmass, double *xheavy, double *xmin);

/* GRAPE relevant funcs */
int vtc_does_include_self_interaction(int calculator);
void vtc_close_grape(void);

#endif /* __VTC_H__ */
