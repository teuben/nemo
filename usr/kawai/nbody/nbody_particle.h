#ifndef  NBODY_PARTICLE_H
#  define   NBODY_PARTICLE_H
/*-----------------------------------------------------------------------------
 *  nbody-particle : basic class for simple nbody implementation
 *  J. Makino 1998/11/29
 *-----------------------------------------------------------------------------
 */
#include <float.h>
#include "vector.h"

#ifndef ONED
#define THREED
#endif
#define REAL_GRAVITY

#ifdef P2M2

#define NPPMAX (500)
#define TDESIGNMAX (100)
typedef struct p2m2_t {
    int order;
    int spherical_design;
    int npp;
    real ppscale;
    vector pppos0[NPPMAX];
    int full_dof;
} P2m2;

#endif // P2M2

class nbody_particle
{
    private:
        vector pos;
	vector vel;
	vector acc_gravity;
	real phi_gravity;
	real mass;
	int index;
	
    public:
	nbody_particle(){
	    pos = 0.0;
	    vel = 0.0;
	    acc_gravity = 0.0;
	    phi_gravity =  mass =  0.0;
	    index = 0;
	}
        void  set_pos(const vector& new_pos)      {pos = new_pos;}
        void  set_vel(const vector& new_vel)      {vel = new_vel;}
        void  set_acc_gravity(const vector& new_acc)      {acc_gravity = new_acc;}
        void  set_phi_gravity(real new_phi)      {phi_gravity = new_phi;}

	void  clear_pos()                         {pos = 0.0;}
	void  clear_vel()                         {vel = 0.0;}
	void  clear_acc_phi_gravity()      {acc_gravity = 0.0;phi_gravity = 0.0;}
	void  correct_phi_self_gravity(real epsinv)      {phi_gravity += mass*epsinv;}

	void  inc_pos(const vector& d_pos)        {pos += d_pos; }
	void  inc_vel(const vector& d_vel)        {vel += d_vel;}
	void  update_vel(real dt)        {vel += acc_gravity*dt;}
	void  update_pos(real dt)        {pos = (pos+vel*dt).readjust();}

	void  scale_pos(const real scale_factor)  {pos *= scale_factor; }
	void  scale_vel(const real scale_factor)  {vel *= scale_factor; }

	vector  get_pos()                         {return pos;}
	vector*  get_posp()                         {return &pos;}
	vector  get_vel()                         {return vel;}
	real  get_phi_gravity()                         {return phi_gravity;}
	vector  get_acc_gravity()                         {return acc_gravity;}

	real get_mass()			{return mass;}
	void set_mass(real m)		{mass = m;}
	void set_index(int i){index = i;}
	int get_index(){return index;}
	void predict(real dt){
	    real dt2 = dt*dt*0.5;
	    pos = (pos + vel*dt + acc_gravity*dt2).readjust();
	    vel +=  acc_gravity*(dt*0.5);
	}
	void correct(real dt){
	    vel +=  acc_gravity*(dt*0.5);
	}

	void read(istream & );
	void write(ostream & );
	void dump();

	void plot(real parm);

	real kinetic_energy();
	real energy();
	real get_ke(){ return 0.5*mass*vel*vel;}

	void friend accumulate_mutual_gravity(nbody_particle & p1,
					      nbody_particle & p2,
					      real eps2);
#ifdef P2M2
	void calculate_gravity_using_tree(real eps2, real theta2, P2m2 *p2m2);
#else // no P2M2
	void calculate_gravity_using_tree(real eps2, real theta2);
#endif // P2M2
	void correct_gravity(real);
};

typedef vector (nbody_particle::*nbody_VMF_ptr)(void);     // vector member function pointer
typedef void (nbody_particle::*nbody_MF_ptr)(const vector &);     // member function pointer

typedef void (nbody_particle::*nbody_VF_ptr)(void);     // void member function pointer
typedef void (nbody_particle::*nbody_RF_ptr)(real);     // void member function
						    // pointer with real arg
typedef void (nbody_particle::*nbody_RRF_ptr)(real,real);     // void member function
						    // pointer with two real args
class nbody_system
    {
    private:
	int nsize;
	nbody_particle * pb;
#ifdef GRAPE5
	real massmin;
	int autoscale;
#endif /* GRAPE5 */

    public:
	int n;
	real time;
	real time0;
	real timestep;
	real eps2_for_gravity;
	int use_self_gravity;
	vector pos;
	vector vel;
	real   mass;
	real plot_xmax;
	real theta_for_tree;
	int ncrit_for_tree;
	int node_div_crit;
#ifdef P2M2
	P2m2 p2m2;
#endif // P2M2
	nbody_system(){
	    n = 0;
	    nsize = 0;
	    time = 0;
	    time0 = 0;
	    timestep = 0;
	    int node_div_crit = 8;
	    pb = NULL;
#ifdef GRAPE5
	    massmin = 0.0;
	    autoscale = 1;
#endif /* GRAPE5 */
	}
 
	void calculate_uncorrected_gravity();
	void calculate_uncorrected_gravity_using_grape();
	void calculate_uncorrected_gravity_direct();

	void read(istream & );
	void write(ostream & );
	void atos(istream & );
	void stoa(ostream & );
#ifdef IONEMO
	void btos(char *fname);
	void stob(char *fo, char *fi);
#endif // IONEMO
	void dump();

#ifdef P2M2
	void  setup_tree(P2m2 *p2m2);
#else // no P2M2
	void  setup_tree();
#endif // P2M2

	void generate_cube(int nx);
	void apply_vf(nbody_VF_ptr);
	void apply_vf(nbody_RF_ptr, real);
	void apply_vf(nbody_RRF_ptr, real, real);
	void plot(real time);
	real kinetic_energy();
	real energy();
	void show_energy();
	void calculate_gravity();
#ifdef TESTCODE
	void calculate_gravity_at(vector evalpos, vector &acc, real &phi); // only for test
	void calculate_gravity_direct_at(vector evalpos, vector &acc, real &phi); // only for test
	void calculate_gravity_classicalME_at(vector evalpos, real &phi); // only for test
	void load_homo_particle(int n);
#endif
	void create_uniform_sphere(int nbody, real power_index, real r0);

	void evolve( real dt, real tend);
	void evolve_onestep( real dt);
	void integrate( real dt);

	void calculate_cmterms();
	void make_collision(vector relpos, vector relv);
	nbody_particle * get_particle_pointer(){return pb;}
	void friend copy_nbody_particles(nbody_system * source,
				       nbody_system * desitination);
	void correct_gravity();

	// library interface
	void set_particle_distribution(int n, real *m, real (*x)[3]);
	void get_particle_force(real (*a)[3], real *p);
	void calculate_gravity_tree(int n, real eps,
				    real *m, real (*x)[3], real (*a)[3], real *p);
	void calculate_gravity_tree_init(real theta, int ncrit, int order);

#ifdef GRAPE5
	real get_massmin(void) {return massmin;}
	void set_massmin(real m) {
	    massmin = m;
	}
	void set_massmin(void){
	    int i;
	    massmin = DBL_MAX;
	    for (i = 0; i < n; i++)
	    {
		real m = (pb+i)->get_mass();
		if (m < massmin && m > 0.0)
		{
		    massmin = m;
		}
	    }
	}
	void set_autoscale(int as){
	    autoscale = as;
	}
	int get_autoscale(void){
	    return (autoscale);
	}
#endif /* GRAPE5 */
#ifdef P2M2
	void load_design(int order, real ss, int full_dof);
#endif // P2M2
};

#endif
