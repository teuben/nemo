//
// gravity.C
//
// Version 1999/2/9 Jun Makino
//
// Change in the calling format for apply_vf, from function name to
// add explicit address operator. This was necessary to pass g++ 2.8.1
// 
// Version 1999/1/1 Jun Makino  --- cleaned up and some comments  added.
//
// The gravitational force calculation package designed for TREE
// implementation of SPH/NBODY programs
//
//
// Structure:
//
//  calculate_gravity (top level routine)
//      calculate_uncorrcted_gravity (calculate gravity with const softening)
//          setup_tree
//	    set_cm_quantities_for_default_tree
//	    calculate_gravity_using_tree,eps2_for_gravity(NON-GRAPE)
//          evaluate_gravity_using_default_tree_and_list(GRAPE)
//      corrcted_gravity (apply SPH form-factor correction)
//

#ifndef NOGRAPHICS
#define GRAPHICS
#endif

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

#include  <stdlib.h>
#include  <math.h>
#include  <iostream>
#include  <fstream>

using namespace std;

#define real double
#include "vector.h"
#include "nbody_particle.h"
#define NBODY

#define NBODY
#ifdef NBODY
typedef nbody_particle real_particle;
typedef nbody_system real_system;
#endif
#ifdef SPH
typedef sph_particle real_particle;
typedef sph_system real_system;
#endif

#ifdef TREE
#ifdef P2M2
void set_cm_quantities_for_default_tree(P2m2 *p2m2);
#else // no P2M2
void set_cm_quantities_for_default_tree();
#endif // P2M2
#endif // TREE   
extern "C" double cpusec();
void clear_tree_counters();
void print_tree_counters();

#ifndef REAL_GRAVITY
void real_system::calculate_gravity()
{
    // This is peudo gravity: just harmonic force
    apply_vf(real_particle::set_harmonic_gravity,1.0);
    
}
#endif
#ifdef REAL_GRAVITY
void real_system::calculate_gravity()
{
    if (use_self_gravity){
	calculate_uncorrected_gravity();
	correct_gravity();
    }
}
#endif

void evaluate_gravity_using_default_tree_and_list(real theta2,
					  real eps2,
					  int ncrit
#ifdef P2M2
						  , P2m2 *p2m2
#endif // P2M2
);

void real_system::calculate_uncorrected_gravity()
{
#ifndef TREE    

#  if defined(HARP3) || defined(GRAPE5)
    calculate_uncorrected_gravity_using_grape();
#  else /* NO GRAPE */
    calculate_uncorrected_gravity_direct();
#  endif /* HARP3 || GRAPE5 */

#else /* TREE */

    cerr << "Call setup tree, cpu = " <<cpusec() << endl;

#ifdef P2M2
    setup_tree(&p2m2);
    cerr << "Call set cm, cpu = " <<cpusec() << endl;
    set_cm_quantities_for_default_tree(&p2m2);
    cerr << "Call evaluate_gravity, cpu = " <<cpusec() << endl;
    clear_tree_counters();
#else // no P2M2
    setup_tree();
    cerr << "Call set cm, cpu = " <<cpusec() << endl;
    set_cm_quantities_for_default_tree();
    cerr << "Call evaluate_gravity, cpu = " <<cpusec() << endl;
    clear_tree_counters();
#endif // P2M2


#  if defined(HARP3) || defined(GRAPE5)
#ifdef P2M2
    evaluate_gravity_using_default_tree_and_list(theta_for_tree*theta_for_tree,
						 eps2_for_gravity,ncrit_for_tree,
						 &p2m2);
#else // no P2M2
    evaluate_gravity_using_default_tree_and_list(theta_for_tree*theta_for_tree,
						 eps2_for_gravity,ncrit_for_tree);
#endif // P2M2

#  else /* NO GRAPE */
    apply_vf(&real_particle::clear_acc_phi_gravity);
#ifdef P2M2
    for(int i = 0;i<n;i++){
	(pb+i)->calculate_gravity_using_tree(eps2_for_gravity,
					    theta_for_tree*theta_for_tree,
					    &p2m2);
    }
#else // no P2M2
    apply_vf(&real_particle::calculate_gravity_using_tree,
	     eps2_for_gravity, theta_for_tree*theta_for_tree);
#endif // P2M2
#  endif /* HARP3 || GRAPE5 */

#endif /* TREE */
    
    print_tree_counters();
    cerr << "Exit evaluate_gravity, cpu = " <<cpusec() << endl;    
}

#if defined(HARP3) || defined(GRAPE5)

#ifdef HARP3
extern "C" void h3open_();
extern "C" void h3close_();
extern "C" void accel_by_harp3_separate_noopen_(int * ni, vector * xi,
						int * nj, vector * xj, real *m,
						vector *  a, real *p, real * eps2);
#endif /* HARP3 */

#ifdef GRAPE5
#include <gp5util.h>
void my_g5_accel(int *ni, double (*xi)[3], int *nj, double (*xj)[3], double *mj,
		 double (*a)[3], double *p, double *eps);
#endif /* GRAPE5 */

void real_system::calculate_uncorrected_gravity_using_grape()
{
    int i;
#if defined(HARP3)
    int jmemsize = 40000;
#elif defined(GRAPE5)
    int jmemsize = JMEMSIZE;
#endif
    static int nharp3 = 0;
    static int call_count = 0;
    static real holdtime = 0.0;
    static int h3_open_state = 0;
    static vector *xgrape = NULL;
    static real *mgrape = NULL;
    static vector *agrape = NULL;
    static real *pgrape = NULL;
    static vector *atmp = NULL;
    static real *ptmp = NULL;

    if (nharp3 < nsize){
	if(nharp3 != 0){
	    delete [] xgrape;
	    delete [] mgrape;
	    delete [] agrape;
	    delete [] pgrape;
	    delete [] atmp;
	    delete [] ptmp;
	}
	nharp3 = nsize;

	xgrape = new vector[nsize+100];
	mgrape = new real[nsize+100];
	agrape = new vector[nsize+100];
	pgrape = new real[nsize+100];
	atmp = new vector[nsize+100];
	ptmp = new real[nsize+100];
    }
    real_particle * p;

    for (i = 0; i < n; i++) {
	agrape[i] = 0.0;
	pgrape[i] = 0.0;
    }
    for (p = pb, i = 0; i < n; p++, i++){
	xgrape[i] = p->get_pos();
	mgrape[i] = p->get_mass();
    }

    int j0 = 0, nj0 = jmemsize;
    while (j0 < n) {
	if (j0+nj0 > n) {
	    nj0 = n-j0;
	}

	if (h3_open_state == 0){
#if defined(HARP3)
	    h3open_();
#elif defined(GRAPE5)
	    g5_open();
#endif
	    holdtime = cpusec();
	    h3_open_state = 1;
	}
#if defined(HARP3)
	accel_by_harp3_separate_noopen_(&n, xgrape, &nj0, xgrape+j0, mgrape+j0,
					atmp, ptmp, &eps2_for_gravity);
#elif defined(GRAPE5)
	double eps = sqrt(eps2_for_gravity);
	my_g5_accel(&n, (double (*)[3])xgrape, &nj0, (double (*)[3])xgrape+j0, mgrape+j0,
		  (double (*)[3])atmp, ptmp, &eps);
#endif
	call_count ++;
	if (call_count * n > 100000){
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

	for (i = 0; i < n; i++) {
	    agrape[i] += atmp[i];
	    pgrape[i] += ptmp[i];
	}
	j0 += nj0;
    }

    real epsinv = 0.0;
    if (eps2_for_gravity != 0.0) {
	epsinv = 1.0/sqrt(eps2_for_gravity);
    }
    for(p = pb, i = 0; i < n; p++, i++){
	p->set_acc_gravity(agrape[i]);
	p->set_phi_gravity(pgrape[i]);
	p->correct_phi_self_gravity(epsinv);
    }
}
#endif /* HARP3 || GRAPE5 */

void real_system::calculate_uncorrected_gravity_direct()
{
    apply_vf(&real_particle::clear_acc_phi_gravity);
    int i, j;
    real_particle * pi;
    real_particle * pj;

    for(i = 0,  pi = &(pb[0]); i<n-1; i++,pi++){
	for(j = i+1,  pj = pi+1; j<n; j++,pj++){
	    accumulate_mutual_gravity(*pi, *pj, eps2_for_gravity);
	}
	if ((i+1)%10000 == 0) {
	    PRL(i);
	}
    }
}

#if 0
template <class T>
void accumulate_mutual_gravity(T & p1,
			       T & p2,
			       real eps2)
#endif    
void accumulate_mutual_gravity(real_particle & p1,
			       real_particle & p2,
			       real eps2)
{
    vector dx = p1.pos-p2.pos;
    double r2inv = 1/(dx*dx+eps2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    p1.phi_gravity -= p2.mass*rinv;
    p2.phi_gravity -= p1.mass*rinv;
    p1.acc_gravity -= p2.mass*r3inv*dx;
    p2.acc_gravity += p1.mass*r3inv*dx;
}

#ifdef SPH
#if 1
// The next function perform the gravity correction with the simplest possible
// method, namely assuming that one particle is uniform density sphere and the
// other a point mass.
//
void apply_sph_correction_for_mutual_gravity(real_particle & p1,
					     real_particle & p2,
					     real eps2,
					     int correct_both)
{
    vector dx = p1.pos-p2.pos;
    real r2 = dx*dx;
    real scaledist = (p1.h + p2.h)*1;
    real scaledist2 = scaledist * scaledist ;

    if (r2 > scaledist2) return;
    real hsinv2 = 1.0/(scaledist2+eps2);
    real hsinv = sqrt(hsinv2);
    real phi_h = hsinv;
    real fcoef = hsinv2*hsinv;
    
    double r2inv = 1/(r2+eps2);
    double r = sqrt(r2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    real phifact = rinv-hsinv+0.5*fcoef*(r2-scaledist2);
    real ffact = r3inv - fcoef;

    p1.phi_gravity += p2.mass*phifact;
    p1.acc_gravity += p2.mass*ffact*dx;
    if(correct_both){
	p2.phi_gravity += p1.mass*phifact;
	p2.acc_gravity -= p1.mass*ffact*dx;
    }
}
#else


// Correction using spline kernel (assuming eps2 = 0)
//
void apply_real_correction_for_mutual_gravity(real_particle & p1,
					     real_particle & p2,
					     real eps2,
					     int correct_both)
{
    vector dx = p1.pos-p2.pos;
    real r2 = dx*dx;
    real scaledist = (p1.h + p2.h);
    real scaledist2 = scaledist * scaledist ;
    //    PRC(r2); PRL(scaledist2);
    if (r2 > scaledist2) return;
    real u2 = r2/scaledist2;
    real u = sqrt(u);
    real f, g;
    real scale3 = scaledist2*scaledist;
    double r2inv = 1/(r2+eps2);
    double r = sqrt(r2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    if  (u < 1){
	f = (1.4 - ((0.05*u-0.015)*u2+1.0/3.0)*u2*2)/scaledist;
	g = ((0.5*u-1.2)*u2+4.0/3.0)/scale3;
    }else{
	f = -rinv/15.0+(1.6-(((-u/30+0.3)*u-1)+4.0/3.0)*u2)/scaledist;
	g = r3inv*((((-u/6+1.2)*u-3)+8.0/3.0)*u2*u-1.0/15.0);
    }
    real phifact = rinv-f;
    real ffact = r3inv - g;
    //    PRC(phifact/rinv); PRC(ffact/r3inv); PRL(dx);
    p1.phi_gravity += p2.mass*phifact;
    p1.acc_gravity += p2.mass*ffact*dx;
    if(correct_both){
	p2.phi_gravity += p1.mass*phifact;
	p2.acc_gravity -= p1.mass*ffact*dx;
    }
}
#endif
#endif

#ifdef SPH
void real_particle::correct_gravity(real eps2)
{
    //    cerr << "correct_gravity for "; PRC(index); PRL(eps2);
    //    PRC(acc_gravity); PRL(phi_gravity);
    for(int j = 0; j<nnb; j++){
	real_particle * pj = pnb[j];
	if (pj->index > index){
	    apply_real_correction_for_mutual_gravity(*this, *pj, eps2, 1);
	}
    }
    //    PRC(acc_gravity); PRL(phi_gravity);
}
#else
void real_particle::correct_gravity(real eps2)
{
}
#endif
void real_system::correct_gravity()
{
    apply_vf(&real_particle::correct_gravity, eps2_for_gravity);
}
