#ifndef  BHTREE_H
#  define   BHTREE_H
/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *  J. Makino 1998/12/14
 *-----------------------------------------------------------------------------
 */

const int default_key_length = 15;
const int default_ix_offset = 1<<(default_key_length-1);
// First include some particle class
// This need, someday, to be changed to template...
#include "vector.h"
#include "nbody_particle.h"
#define NBODY
#ifdef NBODY
typedef nbody_particle real_particle;
typedef nbody_system real_system;
#endif
#ifdef SPH
typedef sph_particle real_particle;
typedef sph_system real_system;
#endif

#ifdef __linux__
typedef  long long BHlong;
#else
typedef  long BHlong;
#endif

class bhparticle
{
private:
    real_particle * rp;
    BHlong key;
    
public:
    bhparticle(){
	rp = NULL;
	key = 0;
    }
    void set_key(real rscale, int ioffset, int keybits);
    int friend compare_key( bhparticle * p1,  bhparticle * p2);
    BHlong get_key(){return key;}
    void set_rp(real_particle * p){rp = p;}
    real_particle * get_rp(){return rp ;}
    void friend sort_bh_array( bhparticle * r, int lo, int up );
    
};


class bhnode
{
private:
    vector pos;
    real l;
    bhnode * child[8];
    bhparticle * bpfirst;
    int nparticle;
    int isleaf;
#ifdef SPH    
    real hmax_for_sph;
#endif    

#ifdef P2M2
    vector *pppos; // pos of pp
    real *ppmass; // mass of pp
#endif // P2M2

    vector cmpos;
    real cmmass;

public:

    bhnode(){
	pos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	bpfirst = NULL;
	nparticle = 0;
	isleaf = 1;
#ifdef SPH	
	hmax_for_sph = 0;
#endif	
#ifdef P2M2
	pppos = NULL;
	ppmass = NULL;
#else // no P2M2
	cmpos = 0.0;
	cmmass = 0.0;
#endif // P2M2
    }
#ifdef P2M2
    void clear(int num_of_pp){
#else // no P2M2
    void clear(){
#endif // P2M2
	pos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	bpfirst = NULL;
	nparticle = 0;
	isleaf = 1;
#ifdef SPH	
	hmax_for_sph = 0;
#endif	

#ifdef P2M2
	for(int i = 0; i<num_of_pp; i++)ppmass[i] = 0.0;
#else // no P2M2
	cmpos = 0.0;
	cmmass = 0.0;
#endif // P2M2
    }
    void set_pos(vector newpos){pos = newpos;}
    vector get_pos(){return pos;}
    vector get_cmpos(){return cmpos;}
    vector*  get_posp(){return &pos;}
    void set_length(real newl){l = newl;}
    real get_length(){return l;}
    void create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
			       BHlong current_key,
			       int current_level,
			       int n_critical);
#ifdef GRAPE5
    void assign_root(vector root_pos, real length, bhparticle * bp, int nparticle,
		     double mmin, int autoscale);
#else
    void assign_root(vector root_pos, real length, bhparticle * bp, int nparticle);
#endif
    void dump(int indent);
    int sanity_check();
#ifdef SPH    
    void set_hmax_for_sph();
    real get_hmax_for_sph(){return hmax_for_sph;}
#endif    
    int friend check_and_set_nbl(bhnode * p1,bhnode * p2);
#ifdef P2M2
    void accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
				   vector & acc, real & phi, P2m2 *p2m2);
    void add_to_interaction_list(bhnode & dest_node, real theta2,
				 vector * pos_list,
				 real * mass_list,
				 int & nlist,
				 int list_max,
				 int & first_leaf,
				 P2m2 *p2m2);
    void evaluate_gravity_using_tree_and_list(bhnode & source_node,
					      real theta2,
					      real eps2,
					      int ncrit,
					      P2m2 *p2m2);
#else // no P2M2
    void accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
				   vector & acc, real & phi);
    void add_to_interaction_list(bhnode & dest_node, real theta2,
				 vector * pos_list,
				 real * mass_list,
				 int & nlist,
				 int list_max,
				 int & first_leaf);
    void evaluate_gravity_using_tree_and_list(bhnode & source_node,
					      real theta2,
					      real eps2,
					      int ncrit);
#endif // P2M2

#ifdef P2M2
    void set_cm_quantities(P2m2 *p2m2);
    void set_cm_quantities_anyorder(P2m2 *p2m2);
    void set_cm_quantities_1storder(P2m2 *p2m2);
    void set_cm_quantities_2ndorder(P2m2 *p2m2);
    void alloc_pp(vector *pbuf, real *mbuf) {
	pppos = pbuf;
	ppmass = mbuf;
    };
    void expand_by_pp(P2m2 *p2m2, bhnode *dst, real *dm);
#else // no P2M2
    void set_cm_quantities();
#endif // P2M2

};

void clear_tree_counters();
void print_tree_counters();
#endif
