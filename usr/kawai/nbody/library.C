#include  <stdlib.h>
#include  <math.h>
#include  <stdiostream.h>

#ifdef GRAPE5
#include <gp5util.h>
#endif /* GRAPE5 */

#define real double
#include "vector.h"
#include "nbody_particle.h"

extern "C" void
nbody_force_gravity(int n, real eps,
		    real *m, real (*x)[3], real (*a)[3], real *p);
extern "C" void
nbody_force_gravity_init(real theta, int ncrit, int order);

extern "C" double cpusec();

void evaluate_gravity_using_default_tree_and_list(real theta2, real eps2, int ncrit
#ifdef P2M2
,P2m2 *p2m2
#endif // P2M2
);
void set_cm_quantities_for_default_tree(
#ifdef P2M2
P2m2 *p2m2
#endif // P2M2
);
void clear_tree_counters();
void print_tree_counters();

static nbody_system pb;

void nbody_system::set_particle_distribution(int n, real *m, real (*x)[3])
{
    this->n = n;
	
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }

    if (pb == NULL){
	pb = new nbody_particle[n];
	nsize = n;
    }
    for (int i = 0; i< n; i++){
	(pb+i)->set_mass(m[i]);
	(pb+i)->set_index(i);
    }
    for (int i = 0; i < n; i++){
	vector vtmp;
	for (int k = 0; k < 3; k++) {
	    vtmp[k] = x[i][k];
	}
	(pb+i)->set_pos(vtmp);
    }

#ifdef GRAPE5
    static int firstcall = 1;
    if (firstcall) {
	firstcall = 0;
	g5_open();
	cerr << "opened" << endl;
	set_massmin();
	cerr << "massmin" << endl;
	g5_set_range(-1.0, 1.0, get_massmin());
	cerr << "setrange" << endl;
	g5_close();
	cerr << "closed" << endl;
    }
#endif // GRAEP5
}

void 
nbody_system::get_particle_force(real (*a)[3], real *p)
{
    for (int i = 0; i < n; i++){
	vector vtmp;
	vtmp = (pb+i)->get_acc_gravity();
	for (int k = 0; k < 3; k++) {
	    a[i][k] = vtmp[k];
	}
    }    
    for (int i = 0; i< n; i++){
	p[i] = (pb+i)->get_phi_gravity();
    }
}

void
nbody_system::calculate_gravity_tree(int n, real eps,
				     real *m, real (*x)[3], real (*a)[3], real *p)
{
    set_particle_distribution(n, m, x);

#ifdef P2M2
    setup_tree(&p2m2);
    set_cm_quantities_for_default_tree(&p2m2);
#else
    setup_tree();
    set_cm_quantities_for_default_tree();
#endif // P2M2
    clear_tree_counters();

#  if defined(HARP3) || defined(GRAPE5)
    evaluate_gravity_using_default_tree_and_list(theta_for_tree*theta_for_tree,
						 eps*eps, ncrit_for_tree
#ifdef P2M2
						 , &p2m2
#endif // P2M2
);
#  else /* NO GRAPE */
    for(int i = 0;i<n;i++){
	(pb+i)->calculate_gravity_using_tree(eps*eps, theta_for_tree*theta_for_tree
#ifdef P2M2
					     ,&p2m2
#endif // P2M2
);
    }
#  endif /* HARP3 || GRAPE5 */

    get_particle_force(a, p);
}

void
nbody_system::calculate_gravity_tree_init(real theta, int ncrit, int order)
{
    ncrit_for_tree = ncrit;
    theta_for_tree = theta;
#ifdef P2M2
    load_design(order, 0.2, 1);
#endif
    node_div_crit = 8;
}

void
nbody_force_gravity_init(real theta, int ncrit, int order)
{
    pb.calculate_gravity_tree_init(theta, ncrit, order);
}

void
nbody_force_gravity(int n, real eps, real *m, real (*x)[3], real (*a)[3], real *p)
{
    pb.calculate_gravity_tree(n, eps, m, x, a, p);
}
