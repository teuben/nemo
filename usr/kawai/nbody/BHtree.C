// 
// BHtree.C
// Version 1999/1/1 --- cleaned up and some comments  added.
//
// The tree construction/handling package designed for TREE
// implementation of SPH/NBODY programs
//
//
// This program uses a new tree construction method, based on
// Morton ordering and top-down tree construction.
//
// non-local functions (not complete yet...)
//
// setup_tree()
//    allocate the memory for tree if necessary
//    create the tree structure
// set_cm_quantities_for_default_tree()
//    calculate mass and position for tree nodes
//
// set_hmax_for_sph() (SPH use only)
// check_and_set_nbl() (SPH use only)
//
// calculate_gravity_using_tree() (NO GRAPE)
//
// evaluate_gravity_using_tree_and_list()(GRAPE)

#define NOFORCE (0)

#define LISTLENCK (0)
// count length of the interaction list with ncrit=1.
// should be 0 for practical use.

extern "C" double cpusec();


#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

#include  <string.h>
#include  <stdlib.h>
#include  <errno.h>
#include  <math.h>
#include  <stdiostream.h>
#ifdef P2M2
#include "p2m2.h"
#endif

#ifdef GRAPE5
#include <gp5util.h>
static double current_xmin;
static double current_xmax;
static double current_mmin;
#endif /* GRAPE5 */

#define real double
#include "BHtree.h"

void calculate_force_from_interaction_list_using_grape(
    vector * pos_list, real * mass_list,
    int list_length, int first_leaf, int ni,
    real eps2,
    vector * acc_list, real * phi_list);

#ifdef PAR_TRAV
#include  <pthread.h>
#include  <sched.h>
#include  <sys/time.h>
#include  <unistd.h>

#define NWORKERS (1)

#if 0
#define check(status,string) //
#elif 1
#define check(status,string) if (status != 0) {         \
    errno = status;                                             \
    fprintf(stderr, "%s status %d: %s\n", string, status, strerror(status));}
#else
#define check(status,string) { \
    if (status != 0) { \
      errno = status; \
      fprintf(stderr, "%s status %d: %s\n", status, string, strerror(status)); \
    } else { \
      fprintf(stderr, "status %d  %s  trav %d calc %d\n", status, string, list_curr, calc_curr); \
    } \
    fflush(stderr);\
}
#endif										 

pthread_mutex_t trav_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  trav_cv = PTHREAD_COND_INITIALIZER;
pthread_t threads[NWORKERS];
static int list_curr;
static int calc_curr;
static int trav_done;
static int add_list_nlist = 0;

static bhparticle * fc_bp[2];
static real (*fc_pos_list[2])[3];
static real * fc_mass_list[2];
static int fc_list_length[2];
static int fc_first_leaf[2];
static int fc_nparticle[2];
static real fc_eps2[2];
static int fc_ncrit[2];

static inline void 
trav_lock(void) 
{ 
    int status = pthread_mutex_lock(&trav_mutex); 
    check(status, "lock"); 
}

static inline void 
trav_unlock(void *) 
{ 
    int status = pthread_mutex_unlock(&trav_mutex); 
    check(status, "unlock"); 
} 

static inline void
trav_signal(char *msg)
{
//    printf(">>> signal %s\n", msg);
    int status = pthread_cond_signal(&trav_cv); 
    check(status, "signal"); 
}

static inline void
trav_sleep(void)
{
    pthread_cleanup_push(trav_unlock, NULL); 

#ifdef __DECCXX
    struct timespec delta, abstime;
    delta.tv_sec = 0;
    delta.tv_nsec = (long)1e8;

    pthread_get_expiration_np(&delta, &abstime);
    int status = pthread_cond_timedwait(&trav_cv, &trav_mutex, &abstime);
    check(status, "awake"); 
#else
    struct timeval tval;
    struct timespec tspec;
    struct timezone dummy;

    gettimeofday(&tval, &dummy);
    TIMEVAL_TO_TIMESPEC(&tval, &tspec);
    tspec.tv_nsec += (long)1e8;
    int status = pthread_cond_timedwait(&trav_cv, &trav_mutex, &tspec);
    check(status, "awake"); 
#endif

    pthread_cleanup_pop(0);
}

#if 1

#define TRAV_WAIT(cond, msg) \
{ \
    pthread_cleanup_push(trav_unlock, NULL); \
    while (cond) \
    {  \
	trav_signal(msg); \
	int status = pthread_cond_wait (&trav_cv, &trav_mutex); \
	check(status, "wait");  \
    } \
    pthread_cleanup_pop(1); \
} \

#else

#define TRAV_WAIT(cond, msg) \
{ \
    pthread_cleanup_push(trav_unlock, NULL); \
    while (cond) \
    {  \
	trav_signal(msg); \
	printf(">>> wait " msg " l: %d c: %d\n", list_curr, calc_curr); \
	int status = pthread_cond_wait (&trav_cv, &trav_mutex); \
	check(status, "wait");  \
    } \
    pthread_cleanup_pop (1); \
} \

#endif

static void*
force_calculator(void* arg)
{
    void *dummy;
    static vector * acc_list = NULL;
    static real * phi_list = NULL;

    while (1)
    {
	// wait until list for next box is prepared
	trav_lock();
	TRAV_WAIT(list_curr == calc_curr, "fc l==c");

	// perform force calculation for the (calc_curr)th
	// acceleration box
	if (acc_list == NULL){
	    acc_list = new vector[fc_ncrit[calc_curr%2] + 100];
	    phi_list = new real[fc_ncrit[calc_curr%2] + 100];
	}
	calculate_force_from_interaction_list_using_grape(
	    (vector *)fc_pos_list[calc_curr%2], fc_mass_list[calc_curr%2],
	    fc_list_length[calc_curr%2], fc_first_leaf[calc_curr%2],
	    fc_nparticle[calc_curr%2], fc_eps2[calc_curr%2], acc_list, phi_list);

	real epsinv = 0.0;
	if (fc_eps2[calc_curr%2] != 0.0) {
	    epsinv = 1.0/sqrt(fc_eps2[calc_curr%2]);
	}
	for(int i = 0; i < fc_nparticle[calc_curr%2]; i++){
	    real_particle * p = (fc_bp[calc_curr%2]+i)->get_rp();
	    p->set_acc_gravity(acc_list[i]);
	    p->set_phi_gravity(phi_list[i] + p->get_mass()*epsinv);
	}

	// increment calc_curr
	trav_lock();
	calc_curr++;
	trav_signal("fc");
	trav_unlock(dummy);
    }
}

#ifdef GRAPE5

void
local_g5_accel(int *ni, double (*xi)[3],   /* ip */
	       int *nj, double (*xj)[3], double *mj, /* jp */
	       double (*a)[3], double *p,  /* val to be returned */
	       double *eps)               /* scalar */
{
    void * dummy;
    int iout = 0;
    int offs, offr, nii, c, c0;
    int i, ic, np, nc;
/*
    static double epsold = -1;
    static double xminold = -1;
    static double xmaxold = -1;
    static double mminold = -1;
    */

    g5_set_mj(0, *nj, mj);
    g5_set_xj(0, *nj, xj);
    g5_set_n(*nj);

#if 0
    if (*eps != epsold ||
	current_xmin != xminold ||
	current_xmax != xmaxold ||
	current_mmin != mminold)
    {
//	cerr << "will g5_set_eps_to_all()" << endl << endl;
	g5_set_eps_to_all(*eps);
	epsold = *eps;
	xminold = current_xmin;
	xmaxold = current_xmax;
	mminold = current_mmin;
    }
#else
	g5_set_eps_to_all(*eps);
#endif


    np = g5_get_number_of_pipelines_per_board();
    nc = g5_get_number_of_boards();
    c0 = g5_get_firstcluster();

    offs = 0;
    for (c = c0; c < nc+c0 && offs < *ni; c++)
    {
	nii = np;
	if (offs+nii > *ni)
	{
	    nii = *ni - offs;
	}
	g5_set_xiMC(c, nii, (double (*)[3])xi[offs]);
	g5_runMC(c);
	offs += nii;
    }

    for (offr = 0; offr < *ni;)
    {
	for (c = c0; c < nc+c0; c++)
	{
	    if (offr < *ni)
	    {
		nii = np;
		if (offr+nii > *ni)
		{
		    nii = *ni - offr;
		}
		g5_get_forceMC(c, nii, (double (*)[3])a[offr], &p[offr]);
		offr += nii;
	    }
	    if (offs < *ni)
	    {
		nii = np;
		if (offs+nii > *ni)
		{
		    nii = *ni - offs;
		}
		g5_set_xiMC(c, nii, (double (*)[3])xi[offs]);
		g5_runMC(c);
		offs += nii;

		if (list_curr != calc_curr+2 && !trav_done)
		{
		    trav_lock();
		    trav_signal("ag5");
		    trav_sleep();
		    trav_unlock(dummy);
		}
	    }
	}
    }

    for (i = 0; i < *ni; i++)
    {
	p[i] *= -1;
    }

/*
    if (*eps != 0.0) {
	for (i = 0; i < *ni; i++) {
	    p[i] += mj[i]/(*eps);
	}
    }
    */
}

#endif // GRAPE5
#endif // PAR_TRAV

void dump_octal(BHlong x)
{
    char st[256];
    sprintf(st," %lo ",x);
    cerr <<  st ;
}
    

BHlong conv_to_morton(int ix, int keybits)
{
    BHlong dum = 0;
    //    char st[256];
    //    cerr << "conv_to_morton "; PR(ix); PRL(keybits);
    //    sprintf(st,"%lo",ix);
    //    cerr << "ix = " << st << endl;
    int i, j;
    for(i = j= 0; i<keybits; i++,j+=3){
	if (ix & (1<<i)){
	    dum |= ((BHlong) 1)<<j;
	}
    }
    //sprintf(st,"%lo",dum);
    //    cerr << "dum = " << st << endl;
    return dum;

}
    
inline BHlong  construct_key(const vector & pos, real rscale, int ioffset, int keybits)
{
    long ix[3];
    vector p = pos;
    for(int i = 0; i<3; i++){
	ix[i] = (long) (p[i]*rscale+ioffset);
    }
    return (conv_to_morton(ix[0],keybits)<<2 )
	|(conv_to_morton(ix[1],keybits)<<1 )
	|(conv_to_morton(ix[2],keybits));
}
    

void bhparticle::set_key(real rscale, int ioffset, int keybits)
{
    key =  construct_key(rp->get_pos(), rscale, ioffset, keybits);
    //    PRL(key);
    
}

int compare_key(bhparticle * p1, bhparticle * p2)
{
    long comp = ((long) p1->get_key()) - ((long) p2->get_key());
    if (comp > 0L){
	return 1;
    }else if (comp == 0L){
	return 0;
    }else{
	return -1;
    }
}

void sort_bh_array( bhparticle * r, int lo, int up )
{
    int i, j;
    bhparticle tempr;
    while ( up>lo ) {
	i = lo;
	j = up;
	tempr = r[lo];
	/*** Split file in two ***/
	while ( i<j ) {
	    for ( ; r[j].key > tempr.key; j-- );
	    for ( r[i]=r[j]; i<j && r[i].key<=tempr.key; i++ );
	    r[j] = r[i];
	}
	r[i] = tempr;
	/*** Sort recursively, the smallest first ***/
	if ( i-lo < up-i ) { sort_bh_array(r,lo,i-1);  lo = i+1; }
	else    { sort_bh_array(r,i+1,up);  up = i-1; }
    }
}
void check_bh_array( bhparticle * r, int size )
{
    for(int i = 0; i<size-1;i++){
	if(r[i].get_key() > r[i+1].get_key()){
	    PR(i); PR(r[i].get_key()); PRL(r[i+1].get_key());
	    cerr << "Sort failed ... \n";
	    exit (1);
	}
    }
}

real initialize_key(int nbody,
		    real_particle * rp,
		    int & nkeysize,
		    bhparticle * &bhp)
{
    if (nbody > nkeysize || bhp == NULL){
	if (bhp != NULL){
	    delete [] bhp;
	}
	nkeysize = nbody+100;
	bhp = new bhparticle[nkeysize];
#ifdef REUSE_PREVIOS_DATA
	// With present quicksort routine, the pre-sorted data
	// actuallt DEGRADE its performance. So DO NOT ACTIVATE
	// THIS PART --- JM 1998/12/22
	for(int i = 0; i<nbody; i++){
	    bhparticle * p = bhp + i;
	    p->set_rp(rp+i);
	}
#endif	
    }
    real rmax = 1.0/64.0;
    for(int i = 0; i<nbody; i++){
	vector p = (rp+i)->get_pos();
	for (int k = 0; k<3; k++){
	    while (fabs(p[k])>=rmax) rmax *= 2;
	}
    }

    real rscale = 1.0/rmax*default_ix_offset;
	
    for(int i = 0; i<nbody; i++){
	bhparticle * p = bhp + i;
#ifndef REUSE_PREVIOS_DATA	
	p->set_rp(rp+i);
#endif	
	p->set_key(rscale, default_ix_offset, default_key_length);
	//	PR(i); PRL(p->get_key());
    }
//    cerr << "Call quicksort, cpu = " << cpusec() << endl;
    sort_bh_array(bhp,0,nbody-1);
    //    qsort(bhp, nbody, sizeof(bhparticle), compare_key);
    // The private sort routine is for some unknow reason
    // much faster than qsort of the system for large N
//    cerr << "Exit quicksort, cpu = " <<cpusec() << endl;
    for(int i = 0; i<nbody; i++){
	bhparticle * p = bhp + i;
	//	PR(i); PR(p->get_key()); PRL(p->get_rp()->get_index());
    }
    return rmax;
}

#ifdef GRAPE5
void bhnode::assign_root(vector root_pos, real length, bhparticle * bp, int np,
			 double mmin, int autoscale)
#else
void bhnode::assign_root(vector root_pos, real length, bhparticle * bp, int np)
#endif
{
    pos = root_pos;
    l = length;
    bpfirst = bp;
    nparticle = np;

#ifdef GRAPE5
    double xmin, xmax;

    xmax = root_pos[0]+length/2.0;
    xmin = root_pos[0]-length/2.0;

    for (int k = 1; k < 3; k++)
    {
	if (xmin > root_pos[k]-length/2.0)
	{
	    xmin = root_pos[k]-length/2.0;
	}
	if (xmax < root_pos[k]+length/2.0)
	{
	    xmax = root_pos[k]+length/2.0;
	}
    }

    if (autoscale) {
	current_xmin = xmin*2.0;
	current_xmax = xmax*2.0;
	current_mmin = mmin;
	g5_set_range(current_xmin, current_xmax, current_mmin);
	cerr << "set_range "; PR(current_xmin); PR(current_xmax); PR(current_mmin);
    }

#endif /* GRAPE5 */

}


    


void bhnode::create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
				   BHlong current_key,
				   int current_level,
				   int n_critical)
{
//    cerr << "create tree called "; PRC(nparticle); PRL(n_critical);
//    PRL(heap_remainder);
    if (heap_remainder <= 0){
	cerr << "create_tree: no more free node... exit\n";
	exit(1);
    }
    if (nparticle <= n_critical) return;
    if (current_level == 0) return;
    //    cerr << "Enter recursion\n";
    //    dump();
    BHlong keyscale = ((BHlong) 1)<<((current_level-1)*3);
    bhparticle * bptmp = bpfirst;
    int npremain = nparticle;
    for(int i=0; i<8;i++)child[i] = NULL;
    isleaf = 1;
    for(int i=0; i<8;i++){
	BHlong new_key = current_key + keyscale * i;
	vector new_pos = pos + vector( ((i&4)*0.5-1)*l/4,
				       ((i&2)    -1)*l/4,
				       ((i&1)*2  -1)*l/4);
	
	if(bptmp->get_key() - new_key <keyscale){
	    // current bptmp is actually in the current subnode
	    // search for the end location
	    int p0 = 0;
	    int p1 = npremain-1;
	    if ((bptmp+p1)->get_key() - new_key >=keyscale){
		while (p1 - p0 > 1){
		    int pnew = (p0+p1)/2;
		    if ((bptmp+pnew)->get_key() - new_key <keyscale){
			p0 = pnew;
		    }else{
			p1 = pnew;
		    }
		}
		p1 = p0;
	    }
	    p1 ++;
	    isleaf = 0;
	    child[i] = heap_top;
	    heap_top ++;
	    heap_remainder -- ;
	    child[i]->bpfirst = bptmp;
	    child[i]->pos = new_pos;
	    child[i]->l = l*0.5;
	    child[i]->nparticle = p1;
	    child[i]->isleaf = 1;
	    child[i]->create_tree_recursive(heap_top, heap_remainder,
					    new_key, current_level-1, n_critical);
	    bptmp += p1;
	    npremain -= p1;
	    if (npremain <= 0) return;
				  
	}
    }

    //dump();
}


void spc(int indent)
{
    for(int i=0;i<indent;i++)cerr << " ";
}

void bhnode::dump(int indent = 0)
{
    int i;
    spc(indent); cerr << "node pos " << pos ;
#ifdef SPH    
    cerr << " h " << hmax_for_sph;
#endif
    cerr << endl;
#ifndef P2M2
    spc(indent); cerr << "node cm  " << cmpos << " m " << cmmass ;
#endif // P2M2
    if (isleaf){
	cerr << " IS LEAF" ;PRL(nparticle);
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    for(int j=0;j<indent+2;j++)cerr << " ";
	    real_particle * p = (bp+i)->get_rp();
	    PR(p->get_index()); PRL(p->get_pos());
	}
    }else{
	cerr << " IS _not_ LEAF ";PRL(nparticle);
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->dump(indent + 2);
	    }
	}
    }
}

// inbox: returns 0 if the particle is in the box;

int inbox(vector  & cpos, // center of the box
	  vector  & pos,  // position of the particle
	  real l)         // length of one side of the box
    
{
    for(int  i = 0; i< ndim; i++){
	if (fabs(pos[i]-cpos[i]) > l*0.5) return 1;
    }
    return 0;
}
	
	
int bhnode::sanity_check()
{
    int i;
    int iret = 0;
    if (isleaf){
	// this is the lowest level node. Things to check:
	// all particles are in the cell
	bhparticle * bp = bpfirst;
	cout << "Leaf np="<< nparticle <<endl;
	for(i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    vector ppos = p->get_pos();
	    if(inbox(pos,ppos,l)){
		cerr << "Error, particle out of box ... \n";
		dump();
		return 1;
	    }
	}
    }else{

	// This is the non-leaf node. Check the position and side
	// length of the child cells and then check recursively..
	cout << "Non Leaf " << pos  <<endl;
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		int err = 0;
	        err = child[i]->sanity_check();
		if (l*0.5 != child[i]->get_length()) err += 2;
		vector relpos = pos-child[i]->get_pos();
		for (int k = 0 ; k<ndim;k++){
		    if (fabs(relpos[k]) !=l*0.25)err += 4;
		}
		if (err){
		    cerr << "Child " << i << " Error type = " << err << endl;
		    dump();
		}
		iret += err;
	    }
	}
    }
    return iret;
}

#ifdef SPH	
void  bhnode::set_hmax_for_sph()
{
    int i;
    hmax_for_sph = 0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real hp = (bp+i)->get_rp()->get_h();
	    if(hmax_for_sph < hp) hmax_for_sph = hp;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_hmax_for_sph();
		if (hmax_for_sph < child[i]->hmax_for_sph)
		    hmax_for_sph = child[i]->hmax_for_sph;
	    }
	}
    }
}
#endif

#ifdef P2M2
// expand potential at dst node

void bhnode::expand_by_pp(P2m2 *p2m2, bhnode *dst, real *dm)
{
    real rr, cosg;
    vector x1;
    real r1, m1;
    static vector x0[NPPMAX];
    static real r0[NPPMAX];
    static real pln[TDESIGNMAX/2+1];

    // x0: pseudo particle pos
    // x1: physical particle pos

    for (int j = 0; j < p2m2->npp; j++) {
	dm[j] = 0.0;
	x0[j] = p2m2->pppos0[j] * dst->get_length();
	r0[j] = abs(x0[j]);
    }

#ifndef TESTCODE
    if (isleaf || nparticle < p2m2->npp) { // this cell does not have pseudo particle
	for (int i = 0; i < nparticle; i++) {
	    x1 = (bpfirst+i)->get_rp()->get_pos() - dst->get_cmpos();
#else
    if (isleaf) { // this cell does not have pseudo particle
	for (int i = 0; i < nparticle; i++) {
	    x1 = (bpfirst+i)->get_rp()->get_pos() - dst->get_pos();
#endif
	    r1 = abs(x1);
	    m1 = (bpfirst+i)->get_rp()->get_mass();
	    for (int j = 0; j < p2m2->npp; j++) {
		rr = r1/r0[j];
		cosg = x0[j]*x1/r0[j]/r1;
		plgndr0(p2m2->spherical_design/2+1, cosg, pln);
		for (int q = 0; q < p2m2->spherical_design/2+1; q++) {
		    dm[j] = dm[j] + pow(rr, (real)q)*m1*(2*q+1)*pln[q];
		}
	    }
	}
    } else {
	for (int i = 0; i < p2m2->npp; i++) {
#ifndef TESTCODE
	    x1 = pppos[i] - dst->get_cmpos();
#else
	    x1 = pppos[i] - dst->get_pos();
#endif
	    r1 = abs(x1);
	    m1 = ppmass[i];
	    for (int j = 0; j < p2m2->npp; j++) {
		rr = r1/r0[j];
		cosg = x0[j]*x1/r0[j]/r1;
		plgndr0(p2m2->spherical_design/2+1, cosg, pln);
		for (int q = 0; q < p2m2->spherical_design/2+1; q++) {
		    dm[j] = dm[j] + pow(rr, (real)q)*m1*(2*q+1)*pln[q];
		}
	    }
	}
    }
    for (int j = 0; j < p2m2->npp; j++) {
	dm[j] /= p2m2->npp;
    }

#if 0
    // check if pseudo-particles give the same expansion coefficient as
    // that given by physical-particles
    if (isleaf) {
	cout << "nparticle: " << nparticle << endl;
	cout << "npp      : " << p2m2->npp << endl << endl;
	for (int l = 0; l < p2m2->order+2; l++) {
	    cout << "l: " << l << endl;

	    cout << "physical: ";
	    for (int m = -l; m <= l; m++) {
		real are = 0.0;
		real aim = 0.0;
		for (int i = 0; i < nparticle; i++) {
		    real re, im;
		    vector xi = (bpfirst+i)->get_rp()->get_pos() - dst->get_pos();
		    real ri = abs(xi);
		    real mi = (bpfirst+i)->get_rp()->get_mass();
		    ylm(l, m, acos(xi[2]/ri), atan2(xi[1], xi[0]),
			&re, &im);
		    are += mi*pow(ri, l)*re;
		    aim += mi*pow(ri, l)*(-im);
		}
		are /= 4.0*M_PI/(2.0*l+1.0);
		aim /= 4.0*M_PI/(2.0*l+1.0);
		cout << are << " " << aim << " ";
	    }

	    cout << endl << "pseudo  : ";
	    for (int m = -l; m <= l; m++) {
		real bre = 0.0;
		real bim = 0.0;
		for (int i = 0; i < p2m2->npp; i++) {
		    real re, im;
		    vector xi = p2m2->pppos0[i] * dst->get_length();
		    real ri = abs(xi);
		    real mi = dm[i];
		    ylm(l, m, acos(xi[2]/ri), atan2(xi[1], xi[0]),
			&re, &im);
		    bre += mi*pow(ri, l)*re;
		    bim += mi*pow(ri, l)*(-im);
		}
		bre /= 4.0*M_PI/(2.0*l+1.0);
		bim /= 4.0*M_PI/(2.0*l+1.0);

		cout << bre << " " << bim << " ";
	    }
	    cout << endl;
	    cout << endl;
	}
    } // isleaf

    static int cnt = 0;
    cnt++;
    if (nparticle > 100)
    {
	exit(1);
    }
#endif // coeff check

}
#endif // P2M2

#ifdef P2M2

void bhnode::set_cm_quantities(P2m2 *p2m2)
{
    if (p2m2->full_dof == 0 || p2m2->order > 2) {
	set_cm_quantities_anyorder(p2m2);
    }
    else if (p2m2->order == 1) {
	set_cm_quantities_1storder(p2m2);
    }
    else if (p2m2->order == 2) {
	set_cm_quantities_2ndorder(p2m2);
    }
    else {
	cerr << "incorrect p2m2->order: " << p2m2->order << endl;
    }
}

void bhnode::set_cm_quantities_anyorder(P2m2 *p2m2)
{
    real dm[NPPMAX];

    // set cm
    cmpos = 0.0;
    cmmass = 0.0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(int i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    cmpos += mchild*(bp+i)->get_rp()->get_pos();
	    cmmass += mchild;
	}
    }else{
	for (int i = 0; i < 8; i++){
	    if (child[i] != NULL){
		child[i]->set_cm_quantities_anyorder(p2m2);
		real mchild = child[i]->cmmass;
		cmpos += mchild*child[i]->cmpos;
		cmmass += mchild;
	    }
	}
    }
    cmpos /= cmmass;

    // set pseudo particles
#ifndef TESTCODE
    if (nparticle >= p2m2->npp)
#endif
    {
	// p2m2 expansion is applied only when the cell contains
	// more than p2m2->npp particle

	if (isleaf) {
	    for (int i = 0; i < p2m2->npp; i++) {
#ifndef TESTCODE
		pppos[i] = cmpos + p2m2->pppos0[i] * l;
#else
		pppos[i] = pos + p2m2->pppos0[i] * l;
#endif
		ppmass[i] = 0.0;
	    }
	    expand_by_pp(p2m2, this, dm);
	    for (int j = 0; j < p2m2->npp; j++) {
		ppmass[j] += dm[j];
	    }
	}
	else {
	    for (int i = 0; i < p2m2->npp; i++) {
#ifndef TESTCODE
		pppos[i] = cmpos + p2m2->pppos0[i] * l;
#else
		pppos[i] = pos + p2m2->pppos0[i] * l;
#endif
		ppmass[i] = 0.0;
	    }
	    for (int i = 0; i < 8; i++) {
		if (child[i] != NULL) {
		    child[i]->expand_by_pp(p2m2, this, dm);
		    for (int j = 0; j < p2m2->npp; j++) {
			ppmass[j] += dm[j];
		    }
		}
	    }
	}
    }
}

void bhnode::set_cm_quantities_1storder(P2m2 *p2m2)
{
    int i;
    cmpos = 0.0;
    cmmass = 0.0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    cmpos += mchild*(bp+i)->get_rp()->get_pos();
	    cmmass += mchild;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_cm_quantities_1storder(p2m2);
		real mchild = child[i]->cmmass;
		cmpos += mchild*child[i]->cmpos;
		cmmass += mchild;
	    }
	}
    }
    cmpos /= cmmass;
    ppmass[0] = cmmass;
    pppos[0] = cmpos;
}

void
printmatrix(real a[3][3])
{
    int i, j;

    cout.precision(10.6);
    for (i = 0; i < 3; i++)
    {
	for (j = 0; j < 3; j++)
	{
	    cout << a[i][j] << " ";
	}
	cout << endl;
    }
}

// #define DEBUG1

void bhnode::set_cm_quantities_2ndorder(P2m2 *p2m2)
{
    int nrot;
    real eval[3];
    real q[3][3];
    real p[3][3];
    real x[3][3];

    // set cm
    cmpos = 0.0;
    cmmass = 0.0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(int i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    cmpos += mchild*(bp+i)->get_rp()->get_pos();
	    cmmass += mchild;
	}
    }else{
	for (int i = 0; i < 8; i++){
	    if (child[i] != NULL){
		child[i]->set_cm_quantities_2ndorder(p2m2);
		real mchild = child[i]->cmmass;
		cmpos += mchild*child[i]->cmpos;
		cmmass += mchild;
	    }
	}
    }
    cmpos /= cmmass;

    // set ppmass
    ppmass[0] = ppmass[1] = ppmass[2] = cmmass / 3.0;

    // set q[][];
    double xx = 0.0, yy = 0.0, zz = 0.0;
    double xy = 0.0, yz = 0.0, zx = 0.0;
    if (isleaf) {
	bhparticle * bp = bpfirst;
	for(int i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    vector pchild = (bp+i)->get_rp()->get_pos() - cmpos;

	    // cout << "nparticle: " << nparticle
	    // << " pchild: " << pchild << endl;

	    xx += mchild*pchild[0]*pchild[0];
	    yy += mchild*pchild[1]*pchild[1];
	    zz += mchild*pchild[2]*pchild[2];
	    xy += mchild*pchild[0]*pchild[1];
	    yz += mchild*pchild[1]*pchild[2];
	    zx += mchild*pchild[2]*pchild[0];
	}
    }
    else {
	for (int i = 0; i < 8; i++){
	    if (child[i] != NULL){
		for (int j = 0; j < p2m2->npp; j++) {
		    real mchild = child[i]->ppmass[j];
		    vector pchild = child[i]->pppos[j] - cmpos;
		    xx += mchild*pchild[0]*pchild[0];
		    yy += mchild*pchild[1]*pchild[1];
		    zz += mchild*pchild[2]*pchild[2];
		    xy += mchild*pchild[0]*pchild[1];
		    yz += mchild*pchild[1]*pchild[2];
		    zx += mchild*pchild[2]*pchild[0];
		}
	    }
	}
    }
    q[0][0] = xx-0.5*yy-0.5*zz;
    q[1][1] = yy-0.5*zz-0.5*xx;
    q[2][2] = zz-0.5*xx-0.5*yy;
    q[0][1] = q[1][0] = 1.5 * xy;
    q[1][2] = q[2][1] = 1.5 * yz;
    q[2][0] = q[0][2] = 1.5 * zx;

#ifdef DEBUG1
    real q0[3][3];
    cout << "q:" << endl;
    for (int i = 0; i < DIM; i++) {
	for (int j = 0; j < DIM; j++) {
	    q0[i][j] = q[i][j];
	}
    }
    printmatrix(q0);
#endif

    jacobi(q, eval, p, &nrot);
    eigenvalsort(eval, p);

#ifdef DEBUG1
    cout << "nrot: " << nrot << endl;
    cout << "eval: " << eval[0] << " "
	 << eval[1] << " " << eval[2] << endl;
#endif

#if 1
    // pseudo-particle 0
    x[0][0] = 0.0;
    x[1][0] = 2.0*sqrt((eval[0]+2.0*eval[1])/3.0/cmmass);
    x[2][0] = 0.0;
    // pseudo-particle 1
    x[0][1] = sqrt((2.0*eval[0]+eval[1])/cmmass);
    x[1][1] = -x[1][0]/2.0;
    x[2][1] = 0.0;
    // pseudo-particle 2
    x[0][2] = -x[0][1];
    x[1][2] = x[1][1];
    x[2][2] = 0.0;
#else
    // pseudo-particle 0
    x[0][0] = 2.0*sqrt((2.0*eval[0]+eval[1])/3.0/cmmass);
    x[1][0] = 0.0;
    x[2][0] = 0.0;
    // pseudo-particle 1
    x[0][1] = -x[0][0]/2.0;
    x[1][1] = sqrt((eval[0]+2.0*eval[1])/cmmass);
    x[2][1] = 0.0;
    // pseudo-particle 2
    x[0][2] = x[0][1];
    x[1][2] = -x[1][1];
    x[2][2] = 0.0;
#endif


    // pppos[0..2] = transpose(p*x)
    for (int i = 0; i < DIM; i++) {
	for (int j = 0; j < DIM; j++) {
	    pppos[j][i] = 0.0;
	    for (int k = 0; k < DIM; k++) {
		pppos[j][i] += p[i][k] * x[k][j];
	    }
	}
    }

#ifdef DEBUG1
    if (isleaf) {
	PRL(pppos[0]);
	PR(l); PRL(nparticle);
	PRL(pppos[0]/l);
    }


    vector hoe;
    hoe = pppos[0] + pppos[1] + pppos[2];
    hoe /= 3.0;
    if (abs(hoe)>1e-10) {
	cout << "-------------------------------------------" << endl;
	cout << "hoe: " << hoe << endl << "cm: " << cmpos << endl;
	cout << "-------------------------------------------" << endl;
    }
    if (cmmass != ppmass[0]+ppmass[1]+ppmass[2]) {
	cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "cmmass: "<< cmmass << endl;
	cout << "ppmass: "<< ppmass[0]+ppmass[1]+ppmass[2] << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
    }


    real q1[3][3];

    xx = yy = zz = xy = yz = zx= 0.0;
    for (int k = 0; k < 3; k++) {
	xx += pppos[k][0]*pppos[k][0];
	yy += pppos[k][1]*pppos[k][1];
	zz += pppos[k][2]*pppos[k][2];
	xy += pppos[k][0]*pppos[k][1];
	yz += pppos[k][1]*pppos[k][2];
	zx += pppos[k][2]*pppos[k][0];
    }
    q1[0][0] = (xx-0.5*yy-0.5*zz)*cmmass/3.0;
    q1[1][1] = (yy-0.5*zz-0.5*xx)*cmmass/3.0;
    q1[2][2] = (zz-0.5*xx-0.5*yy)*cmmass/3.0;
    q1[0][1] = q1[1][0] = (1.5 * xy)*cmmass/3.0;
    q1[1][2] = q1[2][1] = (1.5 * yz)*cmmass/3.0;
    q1[2][0] = q1[0][2] = (1.5 * zx)*cmmass/3.0;

    cout << "q0:" << endl;
    printmatrix(q0);
    cout << "q1:" << endl;
    printmatrix(q1);
#endif


    for (int k = 0; k < 3; k++) {
	pppos[k] += cmpos;
    }

#ifdef DEBUG1
    hoe = ppmass[0]*pppos[0] + ppmass[1]*pppos[1] + ppmass[2]*pppos[2];

    cout << "cmmass: " << cmmass << endl;
    cout << ppmass[0]/cmmass << " " << ppmass[1]/cmmass << " " << ppmass[2]/cmmass << endl;
    cout << "dipole(cm): " << cmmass*cmpos << endl;
    cout << "dipole(pp): " << hoe << endl;
    cout << endl;
#endif    

}

#else // no P2M2
void bhnode::set_cm_quantities()
{
    int i;
    cmpos = 0.0;
    cmmass = 0.0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    cmpos += mchild*(bp+i)->get_rp()->get_pos();
	    cmmass += mchild;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_cm_quantities();
//		cout << "length " << l << " child " << i << endl;
		real mchild = child[i]->cmmass;
		cmpos += mchild*child[i]->cmpos;
		cmmass += mchild;
	    }
	}
    }

    if (cmmass != 0.0) {
	cmpos /= cmmass;
    }else{
	cmpos = pos;
    }
}
#endif // P2M2



inline real separation_squared(bhnode * p1,bhnode * p2)
{
    real r2 = 0;
    real xmin = (p1->get_length()+ p2->get_length())*0.5;
    vector dx = p1->get_pos() - p2->get_pos();
    for (int k = 0; k<ndim; k++){
	real adx = fabs(dx[k]);
	if (adx > xmin){
	    adx -= xmin;
	}else{
	    adx = 0;
	}
	r2 += adx*adx;
    }
    return r2;
}


#if 1

inline real separation_squared0(bhnode * p1, vector & pos2)
{
    real r2 = 0;
    real xmin = p1->get_length()*0.5;
    real adx, ady, adz;
    vector * p = p1->get_posp();
    adx = fabs((*p)[0]-pos2[0]);
    ady = fabs((*p)[1]-pos2[1]);
    adz = fabs((*p)[2]-pos2[2]);
    if (adx>xmin){
	adx -= xmin;
	r2 += adx*adx;
    }
    if (ady>xmin){
	ady -= xmin;
	r2 += ady*ady;
    }
    if (adz>xmin){
	adz -= xmin;
	r2 += adz*adz;
    }
    return r2;
}

inline real separation_squared(bhnode * p1, vector & pos2)
{
    real r2 = 0;
    real xmin = p1->get_length()*0.5;
    real adx, ady, adz;
    vector * p = p1->get_posp();
    adx = fabs((*p)[0]-pos2[0])-xmin;
    ady = fabs((*p)[1]-pos2[1])-xmin;
    adz = fabs((*p)[2]-pos2[2])-xmin;
    if (adx>0.0){
	r2 += adx*adx;
    }
    if (ady>0.0){
	r2 += ady*ady;
    }
    if (adz>0.0){
	r2 += adz*adz;
    }
    
    return r2;
}

#else

real separation_squared(bhnode * p1, vector & pos2)
{
    real r2 = 0;
    real xmin = p1->get_length()*0.5;
    vector dx = p1->get_pos() - pos2;
    for (int k = 0; k<ndim; k++){
	real adx = fabs(dx[k]);
	if (adx > xmin){
	    adx -= xmin;
	}else{
	    adx = 0;
	}
	r2 += adx*adx;
    }
    return r2;
}

#endif

inline int  are_overlapped(bhnode * p1,bhnode * p2)
{
    real xmin = (p1->get_length()+ p2->get_length())*0.4999999999999999;
    vector dx = p1->get_pos() - p2->get_pos();
#if 0
    for (int k = 0; k<ndim; k++){
	if(fabs(dx[k]) > xmin) return 0;
    }
#else
    if((fabs(dx[0]) > xmin) ||(fabs(dx[1]) > xmin) ||(fabs(dx[2]) > xmin)) {
	return 0;
    }
	
#endif
    return 1;
}

#ifdef SPH
int check_and_set_nbl(bhnode * p1,bhnode * p2);
int check_and_set_nbl(bhnode * p1,bhnode * p2)
{
    int iret = 0;
    if(p1 == NULL) return iret;
    if(p2 == NULL) return iret;
    real rcrit = p1->hmax_for_sph;
    if (rcrit <  p2->hmax_for_sph)rcrit = p2->hmax_for_sph;
    rcrit *= 2;
    real rcrit2 = rcrit*rcrit;
    if(separation_squared(p1,p2) > rcrit2) return iret;
    if (p1->isleaf == 0 || p2->isleaf == 0){
	//
	// either node is not leaf. Go down
	//
	if (p1->isleaf || (p2->isleaf == 0) && p2->l > p1->l){
	    //
	    // p1 is already leaf, or both are node but p2 is larger.
	    // go down p2 by 
	    for(int i = 0; i<8;i++)
		iret |= check_and_set_nbl(p1, p2->child[i]);
	}else{
	    //
	    // now, we can go down p1 ...
	    for(int i = 0; i<8;i++)
		iret |= check_and_set_nbl(p1->child[i], p2);
	}
    }else{
	//
	// both are leaves; can directly evaluate particles
	bhparticle * bpi = p1->bpfirst;
	bhparticle * bpj = p2->bpfirst;
	for(int i = 0; i < p1->nparticle; i++)
	    for(int j = 0; j < p2->nparticle; j++){
		if((bpi+i)->get_key() <(bpj+j)->get_key() ){
		    //
		    // here, comparison of key guarantee that a pair
		    // is evaluated only once
	  
		    iret |=check_and_set_nbl(*((bpi+i)->get_rp()),
					     *((bpj+j)->get_rp()));
		}
	    }
    }
    return iret;
}

#endif
static bhparticle * bp = NULL;
static int bhpsize = 0;
static int bnsize = 0;
static bhnode * bn;

#ifdef TESTCODE
// only for test code
bhnode * get_bn(void)
{
    return(bn);
}
#endif

#ifdef P2M2
void set_cm_quantities_for_default_tree(P2m2 *p2m2)
{
    bn->set_cm_quantities(p2m2);
}
#else // no P2M2
void set_cm_quantities_for_default_tree()
{
    bn->set_cm_quantities();
}
#endif // P2M2

#ifdef P2M2

void real_system::setup_tree(P2m2 *p2m2)
{
    real rsize = initialize_key(n,get_particle_pointer(),bhpsize,bp);
//    cerr << "Setup tree: called\n";

#if LISTLENCK
    int expected_bnsize =  (int)(bhpsize*1.6+100);
#else

#ifdef USE_P2M2_ALWAYS
    int expected_bnsize =  (int)(bhpsize*1.5+100);
#else // for practical use
    int expected_bnsize;


#if 0
    if (bhpsize > 1e6) {
	if (p2m2->npp==1) {
	    expected_bnsize=  (int)(bhpsize*0.6);
	}
	else if (p2m2->npp==4) {
	    expected_bnsize=  (int)(bhpsize*0.7/p2m2->npp);
	    PRL(expected_bnsize);
	}
	else if (p2m2->npp==12) {
	    expected_bnsize=  (int)(bhpsize*2.0/p2m2->npp);
	    PRL(expected_bnsize);
	}
	else {
	    expected_bnsize=  (int)(bhpsize*0.6/p2m2->npp);
	}
    }
    else
#endif
{
//    cout << "BHNODE size: " << sizeof(bhnode) << endl;
//    cout << "BHPARTICLE size: " << sizeof(bhparticle) << endl;
//    cout << "NBODY_PARTICLE size: " << sizeof(nbody_particle) << endl;
//    exit(0);

	if (p2m2->npp==1) {
	    expected_bnsize=  (int)(bhpsize*0.8);
	}
	else if (p2m2->npp==4) {
	    expected_bnsize=  (int)(bhpsize*1.5/p2m2->npp);
	    PRL(expected_bnsize);
	}
	else if (p2m2->npp==12) {
	    expected_bnsize=  (int)(bhpsize*2.0/p2m2->npp);
	    PRL(expected_bnsize);
	}
	else {
//	    expected_bnsize=  (int)(bhpsize*1.5/p2m2->npp);
	    expected_bnsize=  (int)(bhpsize*1.5/p2m2->npp);
	}
    }

#endif // USE_P2M2_ALWAYS

#endif

    if (bnsize < expected_bnsize){
	if (bnsize != 0){
	    delete [] bn;
	}
	bnsize = expected_bnsize;
	bn = new bhnode[bnsize];

	static vector * pbuf = new vector[bnsize*p2m2->npp];
	static real *mbuf = new real[bnsize*p2m2->npp];
	for(int j = 0; j<bnsize;j++) {
	    (bn+j)->alloc_pp(pbuf+j*p2m2->npp, mbuf+j*p2m2->npp);
	}
    }
    for(int j = 0; j<bnsize;j++) (bn+j)->clear(p2m2->npp);
#ifdef GRAPE5
    bn->assign_root(vector(0.0), rsize*2, bp, n, massmin, autoscale);
#else
    bn->assign_root(vector(0.0), rsize*2, bp, n);
#endif
    bhnode * btmp = bn+1;
    int heap_remainder = bnsize-1;
    BHlong key = 0;

#if LISTLENCK
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, 1);

#else

#ifdef USE_P2M2_ALWAYS
#ifdef TESTCODE
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, n*2);
#else
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, 1);
#endif

#else // for practical use
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, p2m2->npp*node_div_crit);
#endif // USE_P2M2_ALWAYS

#endif
    PR(bnsize);    PRL(heap_remainder);
    //    PRL(bn->sanity_check());
}

#else // no P2M2

void real_system::setup_tree()
{
    real rsize = initialize_key(n,get_particle_pointer(),bhpsize,bp);
//    cerr << "Setup tree: called\n";

#if LISTLENCK
    int expected_bnsize =  (int)(bhpsize*1.6+100);
#else
//    int expected_bnsize =  (int)(bhpsize*0.6+100);
    int expected_bnsize =  (int)(bhpsize*0.8+100);
//    for GB99
//    int expected_bnsize =  (int)(1260000);
#endif

    if (bnsize < expected_bnsize){
	if (bnsize != 0){
	    delete [] bn;
	}
	bnsize = expected_bnsize;
	bn = new bhnode[bnsize];
    }
    for(int j = 0; j<bnsize;j++) (bn+j)->clear();
#ifdef GRAPE5
    bn->assign_root(vector(0.0), rsize*2, bp, n, massmin, autoscale);
#else
    bn->assign_root(vector(0.0), rsize*2, bp, n);
#endif
    bhnode * btmp = bn+1;
    int heap_remainder = bnsize-1;
    BHlong key = 0;

#if LISTLENCK
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, 1);
#else
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, node_div_crit);
#endif

    PR(bnsize);    PRL(heap_remainder);
    //    PRL(bn->sanity_check());
}

#endif // P2M2


	
#ifdef SPH	
int sph_system::set_nnb_using_tree()
{
    setup_tree();
    real_particle * psph = get_particle_pointer();
    apply_vf(real_particle::clear_nnb);
    bn->set_hmax_for_sph();
    int iret = check_and_set_nbl(bn, bn);
    apply_vf(real_particle::sort_nblist);
    return iret;
}
#endif

void accumulate_force_from_point(vector dx, real r2, real eps2, 
				 vector & acc,
				 real & phi,
				 real jmass)
{
#if 0
    static int firstcall = 1;
    if (!firstcall) return;
    firstcall = 0;
    cout << "hoe" << endl;
    phi = 0;
    acc = 0;
#else
    double r2inv;
    if ((r2+eps2) == 0.0) {
	r2inv = 0;
    }
    else {
	r2inv = 1/(r2+eps2);
    }
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;

//    cout << "phi " << phi << " jmass " << jmass << " rinv " << rinv << endl;

    phi -= jmass*rinv;
    acc += jmass*r3inv*dx;
#endif
}

static real direct_interactions = 0;
static real leaf_cells = 1;

static real total_interactions;
static int tree_walks;
static int nisum;
void clear_tree_counters()
{
    total_interactions = 0;
    tree_walks = 0;
    nisum = 0;
    direct_interactions = 0;
    leaf_cells = 1;
}
void print_tree_counters()
{
#ifdef TREE
    real avg = total_interactions/nisum;
    PRC(nisum); PRC(tree_walks); PRC(total_interactions); PRL(avg);
    PR(direct_interactions); PR(leaf_cells); PRL(direct_interactions/leaf_cells);
//    cout <<"tree_walks = " <<tree_walks << " ntaverage = " << avg << endl;
#endif // TREE
}

void calculate_force_from_interaction_list(const vector & pos,
					   real eps2, 
					    vector & acc,
					    real & phi,
					    vector * poslist,
					    real * masslist,
					    int list_length)
{
    acc = 0.0;
    phi = 0.0;

#if !LISTLENCK

    for(int i = 0; i<list_length; i++){
	vector dx = *(poslist+i)-pos;
	real r2 = dx*dx;
	accumulate_force_from_point(dx, r2, eps2, acc, phi,*(masslist+i));
    }
#endif
}



#ifdef P2M2
void bhnode::accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
					vector & acc, real & phi, P2m2 *p2m2)
{
#ifndef TESTCODE
    vector dx = cmpos - ipos;
#else
    vector dx = pos - ipos;
#endif
    real r2 = dx*dx;

    if (r2*theta2 > l*l){

	// node and position is well separated;

#if defined(USE_P2M2_ALWAYS)
	if (0) {
#else
	if (nparticle <= p2m2->npp) {
#endif
	    bhparticle * bp = bpfirst;
	    for(int i = 0; i < nparticle; i++){
		vector dx = (bp+i)->get_rp()->get_pos()-ipos;
		real r2 = dx*dx;
		accumulate_force_from_point(dx, r2, eps2, acc, phi,
					    (bp+i)->get_rp()->get_mass());
		total_interactions++;
		direct_interactions++;
	    }
	    leaf_cells++;
	} else {
	    for (int i = 0; i < p2m2->npp; i++) {
		vector dx = pppos[i] - ipos;
		real r2 = dx*dx;

		accumulate_force_from_point(dx, r2, eps2, acc, phi, ppmass[i]);
//		cout << "acc " << acc << " phi " << phi
//		     << " ppmass[" << i << "] " << ppmass[i] << endl;
//		PR(pppos[i]); PRL(abs(pppos[i]));
//		PRL(pos);
		total_interactions++;
	    }
	}
    }else{
	int i;
	if (isleaf){
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		vector dx = (bp+i)->get_rp()->get_pos()-ipos;
		real r2 = dx*dx;
		accumulate_force_from_point(dx, r2, eps2, acc, phi,
					    (bp+i)->get_rp()->get_mass());
		total_interactions++;
		direct_interactions++;
	    }
	    leaf_cells++;
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		    child[i]->accumulate_force_from_tree(ipos,eps2,theta2, acc, phi, p2m2);
		}
	    }
	}
    }
}

#else // no P2M2

void bhnode::accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
					vector & acc,
					real & phi)
{
    vector dx = cmpos - ipos;
    real r2 = dx*dx;

    if (r2*theta2 > l*l){
	// node and position is well separated;
	accumulate_force_from_point(dx, r2, eps2, acc, phi, cmmass);
	total_interactions++;
    }else{
	int i;
	if (isleaf){
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		vector dx = (bp+i)->get_rp()->get_pos()-ipos;
		real r2 = dx*dx;
		accumulate_force_from_point(dx, r2, eps2, acc, phi,
					    (bp+i)->get_rp()->get_mass());
		total_interactions++;
	    }
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		    child[i]->accumulate_force_from_tree(ipos,eps2,theta2, acc, phi);
		}
	    }
	}
    }
}

#endif // P2M2

void bhnode::add_to_interaction_list(bhnode & dest_node, real theta2,
				     vector * pos_list,
				     real * mass_list,
				     int & nlist,
				     int list_max,
				     int & first_leaf
#ifdef P2M2
, P2m2 *p2m2
#endif // P2M2
)
{
    void *dummy;

#ifdef PAR_TRAV

#if 0
    static int add_list_nlist_max = 500;
    add_list_nlist++;
    if (add_list_nlist > add_list_nlist_max && list_curr > calc_curr)
    {
	add_list_nlist = 0;
	add_list_nlist_max = (10 + fc_list_length[calc_curr%2]*0.03);
	trav_lock();
	trav_signal("add_to");
	trav_sleep();
	trav_unlock();
    }
#else
    if (add_list_nlist < nlist && list_curr > calc_curr)
    {
//	add_list_nlist += 300;
	add_list_nlist += (int)(10 + fc_list_length[calc_curr%2]*0.03);
	trav_lock();
	trav_signal("add_to");
	trav_sleep();
	trav_unlock(dummy);
    }
#endif
#endif // PAR_TRAV

    // check if node and position are well separated
#ifdef P2M2
//    if(!are_overlapped(this,&dest_node) && (separation_squared(&dest_node,pos)*theta2 > l*l)){
    if(!are_overlapped(this,&dest_node) && (separation_squared(&dest_node,cmpos)*theta2 > l*l)){
#else // no P2M2
    if(!are_overlapped(this,&dest_node) && (separation_squared(&dest_node,cmpos)*theta2 > l*l)){
#endif // P2M2
#ifdef P2M2
	if (nparticle <= p2m2->npp) {
	    for (int j = 0; j < nparticle; j++) {
		*(pos_list+nlist) = (bpfirst+j)->get_rp()->get_pos();
		*(mass_list+nlist) = (bpfirst+j)->get_rp()->get_mass();
		nlist++;
		direct_interactions++;
	    }
	    leaf_cells++;
	} else {
	    for (int j = 0; j < p2m2->npp; j++) {
		*(pos_list+nlist) = pppos[j];
		*(mass_list+nlist) = ppmass[j];
		nlist++;
	    }
	}
#else // no P2M2
	*(pos_list+nlist) = cmpos;
	*(mass_list+nlist) = cmmass;
	nlist ++;
#endif // P2M2
    }else{
	int i;
	if (isleaf || this == (&dest_node)) {
	    if (this == (&dest_node)) {
		// adding the particles in the node itself
		first_leaf = nlist;
	    }
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		*(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
		*(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
		nlist ++;
	    }
	} else {
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		    child[i]->add_to_interaction_list(dest_node, theta2,
						      pos_list, mass_list,
						      nlist, list_max,
						      first_leaf
#ifdef P2M2
						      , p2m2
#endif // P2M2
);
		}
	    }
	}
	if (nlist > list_max){
	    cerr << "List length exceeded " << nlist << endl;
	    exit(1);
	}
    }
}



#if defined(HARP3) || defined(GRAPE5)

#ifdef HARP3
extern "C" void h3open_();
extern "C" void h3close_();
extern "C" int get_nboards();
extern "C" void accel_by_harp3_separate_noopen_(int * ni, vector * xi,
						int * nj, vector * xj, real *m,
						vector *  a, real *p, real * eps2);
void
my_h3_accel(int *ni, vector *xi,
	    int *nj, vector *xj, real *m,
	    vector *a, real *p, real *eps2)
{
    int memsize = 43000*get_nboards();

    if (*nj <= memsize) {
	accel_by_harp3_separate_noopen_(ni, xi,	nj, xj, m, a, p, eps2);
	return;
    }

    int j0 = 0, nj0 = memsize;
    static int old_size = 0;
    static vector *atmp = NULL;
    static real *ptmp = NULL;

    if (old_size < *ni) {
	old_size = ((*ni)/10000+1)*10000;
	if (atmp != NULL) {
	    delete atmp;
	    delete ptmp;
	}
	atmp = new vector[old_size];
	ptmp = new real[old_size];
	if (atmp == NULL || ptmp == NULL) {
	    cerr << "BHtree.C line " << __LINE__
		 << ": memory allocation failed." << endl;
	    exit(1);
	}
    }

    for (int i = 0; i < *ni; i++) {
	a[i] = 0.0;
	p[i] = 0.0;
    }
    while (j0 < *nj) {
	if (j0+nj0 > *nj) {
	    nj0 = *nj-j0;
	}
//	cerr << endl;
//	PR(*ni); PR(j0); PR(nj0); PR(old_size); PR(*nj);
//	cerr << endl;

	accel_by_harp3_separate_noopen_(ni, xi,	&nj0, xj+j0, m+j0,
					atmp, ptmp, eps2);

	for (int i = 0; i < *ni; i++) {
	    a[i] += atmp[i];
	    p[i] += ptmp[i];
	}
	j0 += nj0;
    }

}

#endif /* HARP3 */


#ifdef GRAPE5
void
my_g5_accel(int *ni, double (*xi)[3],   /* ip */
		 int *nj, double (*xj)[3], double *mj, /* jp */
		 double (*a)[3], double *p,  /* val to be returned */
		 double *eps)               /* scalar */
{
    int iout = 0;
    int offs, offr, nii, c, c0;
    int i, ic, np, nc;
/*
    static double epsold = -1;
    static double xminold = -1;
    static double xmaxold = -1;
    static double mminold = -1;
    */

    if (JMEMSIZE < *nj) {
	cerr << "nj " << *nj << "exceeded GRAPE-5 JMEMSIZE(" << JMEMSIZE << ")" << endl;
	exit(1);
    }

    g5_set_mj(0, *nj, mj);
    g5_set_xj(0, *nj, xj);
    g5_set_n(*nj);

#if 0
    if (*eps != epsold ||
	current_xmin != xminold ||
	current_xmax != xmaxold ||
	current_mmin != mminold)
    {
	cerr << "will g5_set_eps_to_all()" << endl << endl;
	g5_set_eps_to_all(*eps);
	epsold = *eps;
	xminold = current_xmin;
	xmaxold = current_xmax;
	mminold = current_mmin;
    }
#else
	g5_set_eps_to_all(*eps);
#endif

    np = g5_get_number_of_pipelines_per_board();
    nc = g5_get_number_of_boards();
    c0 = g5_get_firstcluster();

    offs = 0;
    for (c = c0; c < nc+c0 && offs < *ni; c++)
    {
	nii = np;
	if (offs+nii > *ni)
	{
	    nii = *ni - offs;
	}
	g5_set_xiMC(c, nii, (double (*)[3])xi[offs]);
	g5_runMC(c);
	offs += nii;
    }

    for (offr = 0; offr < *ni;)
    {
	for (c = c0; c < nc+c0; c++)
	{
	    if (offr < *ni)
	    {
		nii = np;
		if (offr+nii > *ni)
		{
		    nii = *ni - offr;
		}
		g5_get_forceMC(c, nii, (double (*)[3])a[offr], &p[offr]);
		offr += nii;
	    }
	    if (offs < *ni)
	    {
		nii = np;
		if (offs+nii > *ni)
		{
		    nii = *ni - offs;
		}
		g5_set_xiMC(c, nii, (double (*)[3])xi[offs]);
		g5_runMC(c);
		offs += nii;
	    }
	}
    }
    for (i = 0; i < *ni; i++)
    {
	p[i] *= -1;
    }
}

#endif // GRAPE5

void calculate_force_from_interaction_list_using_grape(
    vector * pos_list, real * mass_list,
    int list_length, int first_leaf, int ni,
    real eps2,
    vector * acc_list, real * phi_list)
{
    static int call_count = 0;
    static real holdtime = 0.0;
    static int h3_open_state = 0;
    if (h3_open_state == 0){
#if defined(HARP3)
	h3open_();
#elif defined(GRAPE5)
	g5_open();
#endif
	holdtime = cpusec();
	h3_open_state = 1;
    }
    //    PR(ni);PRL(list_length);
    nisum += ni;
    tree_walks += 1;
    total_interactions += ((real)ni)*list_length;
#if defined(HARP3)

    // split interaction list if its length exceeds total size of the
    // particle memory on GRAPE-4 PBs
    my_h3_accel(&ni,pos_list+first_leaf,
		&list_length,pos_list, mass_list,
		acc_list, phi_list, &eps2);
/*
    accel_by_harp3_separate_noopen_(&ni,pos_list+first_leaf,
				    &list_length,pos_list, mass_list,
				    acc_list, phi_list, &eps2);
				    */
#elif defined(GRAPE5)

#if LISTLENCK  // performs no force calculation
    if (tree_walks%100000 == 0)
    {
	cerr << "done " << tree_walks << endl;
    }
    for (int i = 0; i < ni; i++)
    {
	for (int k = 0; k < 3; k++)
	{
	    acc_list[i][k] = 0.0;
	}
	phi_list[i] = 0.0;
    }
#else

#ifdef PAR_TRAV
    double eps = sqrt(eps2);
    local_g5_accel(&ni,(double (*)[3])(pos_list+first_leaf),
		     &list_length, (double (*)[3])pos_list, mass_list,
		     (double (*)[3])acc_list, phi_list, &eps);
#else // !PAR_TRAV
    double eps = sqrt(eps2);
    my_g5_accel(&ni,(double (*)[3])(pos_list+first_leaf),
		     &list_length, (double (*)[3])pos_list, mass_list,
		     (double (*)[3])acc_list, phi_list, &eps);
#endif // PAR_TRAV

#endif

#endif
    call_count += ni;
    if (call_count > 50000) {
	real now = cpusec();
	if (now - holdtime > 60.0) {
#if defined(HARP3)
	    cerr << "Close and release GRAPE-4\n";
	    h3close_();
#elif defined(GRAPE5)
	    cerr << "Close and release GRAPE-5\n";
	    g5_close();
#endif
	    h3_open_state = 0;
	}
	call_count = 0;
    }
}

static int h3_open_state = 0;

#endif /* HARP3 || GRAPE5 */

void bhnode::evaluate_gravity_using_tree_and_list(bhnode & source_node,
						  real theta2,
						  real eps2,
						  int ncrit
#ifdef P2M2
						  , P2m2 *p2m2
#endif // P2M2
)
{
  void *dummy;

  // memory allocation for the interaction list buffer

  // cal necessary size
#ifdef GRAPE5
    static int nc = 0;
    static int firstcallg5 = 1;
    if (firstcallg5) {
	firstcallg5 = 0;
	nc = g5_get_number_of_boards();
    }
    int list_max = JMEMSIZE*nc;
#else
    int list_max = 30000;
#endif

#ifdef PAR_TRAV
      //	static real pos_list[2][list_max][3];
      //	static real mass_list[2][list_max];
    static real (*pos_list[2])[3];
    static real *mass_list[2];

    if (mass_list[0] == NULL) {
	pos_list[0] = new real[list_max][3];
	pos_list[1] = new real[list_max][3];
	mass_list[0] = new real[list_max];
	mass_list[1] = new real[list_max];
    }
#else // !PAR_TRAV
    static real *mass_list = NULL;
    static vector *pos_list = NULL;
    static vector *acc_list = NULL;
    static real *phi_list = NULL;

    if (mass_list == NULL) {
	mass_list = new real[list_max];
	pos_list = new vector[list_max];
	acc_list = new vector[ncrit + 100];
	phi_list = new real[ncrit + 100];
    }
#endif // PAR_TRAV

    real epsinv = 0.0;
    if (eps2 != 0.0) {
	1.0/sqrt(eps2);
    }

#ifdef PAR_TRAV

    //    PR(pos); PR(nparticle); PRL(isleaf);
    if((nparticle > ncrit) && (isleaf==0)){
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->evaluate_gravity_using_tree_and_list(source_node,
							       theta2,
							       eps2,
							       ncrit
#ifdef P2M2
							       , p2m2
#endif // P2M2
);
	    }
	}
    }else{

	trav_lock();
	TRAV_WAIT(list_curr > calc_curr+1, "ev_grav l>c+1");

	//
	// node is below critical ... first create list
	//
	int list_length = 0;
	int first_leaf = -1;

//	fprintf(stderr, ">>> will add_to...   trav: %d, calc: %d\n", list_curr, calc_curr);

	add_list_nlist = 0;
	source_node.add_to_interaction_list(*this,  theta2,
					    (vector *)pos_list[list_curr%2],
					    mass_list[list_curr%2],
					    list_length,
					    list_max,
					    first_leaf
#ifdef P2M2
					    , p2m2
#endif // P2M2
);

//	fprintf(stderr, ">>> done add_to... list_length: %d\n", list_length);



	if (first_leaf == -1){
	    cerr << "evaluate_gravity: impossible error \n";
	    cerr << "failed to find the node in the tree \n";
	    exit(1);
	}

	trav_lock();

	fc_bp[list_curr%2] = bpfirst;
	fc_mass_list[list_curr%2] =mass_list[list_curr%2];
	fc_pos_list[list_curr%2] = pos_list[list_curr%2];
	fc_list_length[list_curr%2] = list_length;
	fc_first_leaf[list_curr%2] = first_leaf;
	fc_nparticle[list_curr%2] = nparticle;
	fc_eps2[list_curr%2] = eps2;
	fc_ncrit[list_curr%2] = ncrit;
	list_curr++;

	trav_signal("ev_grav");
	trav_unlock(dummy);

	// this thread does not perform any force calculation.
    }

#else // !PAR_TRAV

    //    PR(pos); PR(nparticle); PRL(isleaf);
    if((nparticle > ncrit) && (isleaf==0)){
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->evaluate_gravity_using_tree_and_list(source_node,
							       theta2,
							       eps2,
							       ncrit
#ifdef P2M2
							       , p2m2
#endif // P2M2
);
	    }
	}
    }else{

	//
	// node is below critical ... first create list
	//
	int list_length = 0;
	int first_leaf = -1;
	source_node.add_to_interaction_list(*this,  theta2,
					    pos_list,
					    mass_list,
					    list_length,
					    list_max,
					    first_leaf
#ifdef P2M2
					    , p2m2
#endif // P2M2
					    );
	if (first_leaf == -1){
	    cerr << "evaluate_gravity: impossible error \n";
	    cerr << "failed to find the node in the tree \n";
	    exit(1);
	}

#if defined(HARP3) || defined(GRAPE5)
	bhparticle *bp = bpfirst;

#if NOFORCE // only for performance checking
        for(int i = 0; i < nparticle; i++){
          acc_list[i] = 0.0;
          phi_list[i] = 0.0;
          nisum = 1;
        }
#else
	calculate_force_from_interaction_list_using_grape(
	    pos_list, mass_list,list_length, first_leaf,
	    nparticle, eps2, acc_list, phi_list);
#endif

	for(int i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    p->set_acc_gravity(acc_list[i]);
	    p->set_phi_gravity(phi_list[i] + p->get_mass()*epsinv);
	}

#else /* NO GRAPE */
	for(int i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    vector acc;
	    real phi;
	    calculate_force_from_interaction_list(pos_list[i+first_leaf],eps2, acc, phi,
					  pos_list,mass_list,list_length);
	    p->set_acc_gravity(acc);
	    p->set_phi_gravity(phi + p->get_mass()*epsinv);
	}
#endif /* HARP3 || GRAPE5 */
    }
#endif // PAR_TRAV
}

#ifdef PAR_TRAV

void create_force_calculator(void)
{
    int status = 0;
    static int firstcall = 1;
    struct sched_param param;

    list_curr = 0;
    calc_curr = 0;
    trav_done = 0;
    if (firstcall)
    {
	firstcall = 0;
	for (int w = 0; w < NWORKERS; w++)
	{ 
	    fprintf(stderr, "will pthread_create: %016x\n", threads[w]);

#ifdef __DECCXX
	    param.sched_priority = PRI_FIFO_MAX;
#else
	    param.sched_priority = sched_get_priority_max(SCHED_FIFO);
#endif
	    pthread_setschedparam(pthread_self(), SCHED_FIFO, &param);
//	    pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
	    status = pthread_create(&threads[w], NULL, force_calculator, (void*)w); 
	    check(status, "pthread_create"); 
	    fprintf(stderr, "pthread_create: %016x\n\n", threads[w]);
	} 
    }
}

#endif // PAR_TRAV

void evaluate_gravity_using_default_tree_and_list(real theta2,
					  real eps2,
					  int ncrit
#ifdef P2M2
						  , P2m2 *p2m2
#endif // P2M2
)
{
#ifdef PAR_TRAV
    void *dummy;

    trav_lock();
    list_curr = 0;
    calc_curr = 0;
    trav_done = 0;
    trav_unlock(dummy);

    bn->evaluate_gravity_using_tree_and_list(*bn, theta2,eps2, ncrit
#ifdef P2M2
					     , p2m2
#endif // P2M2
);

    trav_lock();
    trav_done = 1;
    TRAV_WAIT(list_curr > calc_curr, "ev_grav_def l>c");

#else // !PAR_TRAV

    bn->evaluate_gravity_using_tree_and_list(*bn, theta2,eps2, ncrit
#ifdef P2M2
					     , p2m2
#endif // P2M2
);

#endif // PAR_TRAV
}

#ifdef P2M2
void real_particle::calculate_gravity_using_tree(real eps2, real theta2, P2m2 *p2m2)
#else // no P2M2
void real_particle::calculate_gravity_using_tree(real eps2, real theta2)
#endif // P2M2
{
    acc_gravity = 0;

    if (eps2 == 0.0) {
	phi_gravity = 0.0;
    }
    else
    {
	phi_gravity = mass/sqrt(eps2);
    }

#ifdef P2M2
    bn->accumulate_force_from_tree(pos,eps2,theta2,
				  acc_gravity, phi_gravity, p2m2);
#else // no P2M2
    bn->accumulate_force_from_tree(pos,eps2,theta2,
				  acc_gravity, phi_gravity);
#endif // P2M2

    nisum += 1;

#if LISTLENCK
    if (nisum%10000 == 0)
    {
	cerr << "nisum " << nisum << " ntavg "
	     << total_interactions/nisum << endl;
    }
#endif
}


#ifdef TESTXXX
//
// do the serious test of
// construction of tree
// consistency of tree
// validity of the neighbour list (by comparing with the result
// of direct calculation

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    pb.initialize_h_and_nbl(pow(1.0/n,0.33333));
    static sph_system pbcopy = pb;
    copy_sph_particles(&pb, &pbcopy);
    real_particle * psph = pb.get_particle_pointer();
    real_particle * psphcopy = pbcopy.get_particle_pointer();
    pb.set_nnb_using_tree();
    cerr << "Dumping copy ... \n";
    cerr << "checking NB \n";
    int error = 0;
    for(int i = 0; i<n; i++){
	(psph+i)->sort_nblist();
	int err = 0;
	if((psph+i)->get_nnb() != (psphcopy+i)->get_nnb()){
	    cerr << "Neighbour count differs for "; PRL(i);
	    err = 1;
	    
	}
	if (err == 0){
	    for(int j = 0; (j< (psph+i)->get_nnb()) && (err == 0); j++){
		if ((psph+i)->get_neighbor(j)->get_index()!=
		    (psphcopy+i)->get_neighbor(j)->get_index()) err = 1;
	    }
	}
	if(err){
	    (psph+i)->dump();
	    (psphcopy+i)->dump();
	    error ++;
	}
    }
    PRL(error);
}
#endif /* TESTXXX */

#ifdef TEST
//
// do the serious test of
// construction of tree
// consistency of tree
// validity of the neighbour list (by comparing with the result
// of direct calculation

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    bhparticle * bp = NULL;

    bn = new bhnode[n];
    for(int i = 0; i<1; i++){
	real rsize = initialize_key(n,pb.get_particle_pointer(),nkey,bp);
	for(int j = 0; j<n;j++) (bn+j)->clear();
#ifdef GRAPE5
	bn->assign_root(vector(0.0), rsize*2, bp, n, massmin, autoscale);
#else
	bn->assign_root(vector(0.0), rsize*2, bp, n);
#endif
	bhnode * btmp = bn+1;
	int heap_remainder = n-1;
	BHlong key = 0;
        bn->create_tree_recursive(btmp,  heap_remainder,key, default_key_length, 4);
    }
    PRL(bn->sanity_check());
    pb.initialize_h_and_nbl(pow(1.0/n,0.33333));
    bn->set_hmax_for_sph();
    //    bn->dump();
    bn->set_cm_quantities();
    //    bn->dump();
    static sph_system pbcopy = pb;
    copy_sph_particles(&pb, &pbcopy);
    real_particle * psph = pb.get_particle_pointer();
    real_particle * psphcopy = pbcopy.get_particle_pointer();
    for(int i = 0; i<n; i++){
	(psph+i)->clear_nnb();
    }
    PRL(check_and_set_nbl(bn, bn));
    cerr << "Dumping copy ... \n";
    cerr << "checking NB \n";
    int error = 0;
    for(int i = 0; i<n; i++){
	(psph+i)->sort_nblist();
	int err = 0;
	if((psph+i)->get_nnb() != (psphcopy+i)->get_nnb()){
	    cerr << "Neighbour count differs for "; PRL(i);
	    err = 1;
	    
	}
	if (err == 0){
	    for(int j = 0; (j< (psph+i)->get_nnb()) && (err == 0); j++){
		if ((psph+i)->get_neighbor(j)->get_index()!=
		    (psphcopy+i)->get_neighbor(j)->get_index()) err = 1;
	    }
	}
	if(err){
	    (psph+i)->dump();
	    (psphcopy+i)->dump();
	    error ++;
	}
    }
    PRL(error);
    pb.use_self_gravity = 1;
    pb.eps2_for_gravity = 0.01;
#define COMPARISON_WITH_DIRECT    
#ifdef COMPARISON_WITH_DIRECT    
    pb.calculate_uncorrected_gravity_direct();
    copy_sph_particles(&pb, &pbcopy);
    psphcopy = pbcopy.get_particle_pointer();
    cerr << "Direct force \n";
    for(int i = 0; i<n; i++){
	real phi = (psphcopy+i)->get_phi_gravity();
	vector acc  = (psphcopy+i)->get_acc_gravity();
	PR(i); PR(phi); PRL(acc);
    }
#endif


    
    cerr << "Tree   force \n";
    for(int j = 0; j<10; j++){
	PRL(j);
	pb.apply_vf(real_particle::clear_acc_phi_gravity);
	for(int i = 0; i<n; i++){
	    (psph+i)->calculate_gravity_using_tree(pb.eps2_for_gravity, 0.4);
	}
    }
    pb.apply_vf(real_particle::clear_acc_phi_gravity);
    bn->evaluate_gravity_using_tree_and_list(*bn,0.4,pb.eps2_for_gravity,1);
#ifdef COMPARISON_WITH_DIRECT    
    real perrmax = 0;
    real ferrmax = 0;
    for(int i = 0; i<n; i++){
	real phi = (psph+i)->get_phi_gravity();
	real phierr = (psphcopy+i)->get_phi_gravity()-phi;
	vector acc  = (psph+i)->get_acc_gravity();
	vector accerr  = (psphcopy+i)->get_acc_gravity()-acc;
	PR(i); PR(phi); PRC(acc); PRC(phierr); PRL(accerr);
	real prelerr = fabs(phierr/phi);
	real frelerr = abs(accerr)/abs(acc);
	if(perrmax < prelerr) perrmax = prelerr;
	if(ferrmax < frelerr) ferrmax = frelerr;
    }
    PR(perrmax);    PRL(ferrmax);
#else
    for(int i = 0; i<n; i++){
	real phi = (psph+i)->get_phi_gravity();
	vector acc  = (psph+i)->get_acc_gravity();
	PR(i); PR(phi); PRL(acc); 
    }
	
#endif    
    
}
#endif /* TEST */

	
#ifdef TESTXX
//
// Sample test for timing purpose...

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    bhparticle * bp = NULL;

    bn = new bhnode[n];
    for(int i = 0; i<10; i++){
	real rsize = initialize_key(n,pb.get_particle_pointer(),nkey,bp);
	for(int j = 0; j<n;j++) (bn+j)->clear();
#ifdef GRAPE5
	bn->assign_root(vector(0.0), rsize*2, bp, n, massmin, autoscale);
#else
	bn->assign_root(vector(0.0), rsize*2, bp, n);
#endif
	bhnode * btmp = bn+1;
	int heap_remainder = n-1;
	BHlong key = 0;
        bn->create_tree_recursive(btmp,
				  heap_remainder,key,
				  default_key_length, 8 );
	PRL(heap_remainder);
    }
    PRL(bn->sanity_check());
    real_particle * psph = pb.get_particle_pointer();
    real h0 = pow(1.0/n,0.33333);
    for(int i = 0; i<10; i++){
	
	pb.apply_vf(real_particle::set_h, h0);
	pb.apply_vf(real_particle::clear_nnb);
	bn->set_hmax_for_sph();
	//    bn->dump();
	PRL(check_and_set_nbl(bn, bn));
	pb.apply_vf(real_particle::sort_nblist);
		
    }
}
#endif /* TESTXX */
