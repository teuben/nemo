/**
 * rebound integrate, NEMO style
 *
 * Originally:
 * A self-gravitating Plummer sphere is integrated using
 * the leap frog integrator. Collisions are not resolved. Note that the
 * fixed timestep might not allow you to resolve individual two-body
 * encounters. An alternative integrator is IAS15 which
 * comes with adaptive timestepping.
 */
#include <nemo.h>
#include "rebound.h"

string defv[] = {
    "in=???\n                     input snapshot",
    "out=???\n                    output snapshot",
    "tstop=10\n                   Set, or let go forever",
    "dtout=1\n                    Optional interval to save simulation",
    "dt=0.01\n                    Integration Step",
    "eps=0.01\n                   Gravitational Softening Length",
    "integrator=leapfrog\n        rebound integrator",
    "gravity=basic\n              none, basic, compensated or tree",
    "tolerance=0.25\n             tree opening angle - in radians?",
    "box=\n                       Optional Simulation Box size - in terms or r0",
    "G=1\n                        Gravitational Constant",
    "server=-1\n                  Use web based server on this port (if > 0)",
    "VERSION=0.3\n                17-jan-2024 PJT",
    NULL,
};

string usage="Rebound nbody integrator for NEMO";

void heartbeat(struct reb_simulation* r);
    
int integrator_type(string ri)
{
  if (streq(ri,"ias15"))     return REB_INTEGRATOR_IAS15;      // 0 - IAS15 integrator, 15th order, non-symplectic (default)
  if (streq(ri,"whfast"))    return REB_INTEGRATOR_WHFAST;     // 1 - WHFast integrator, symplectic, 2nd order, up to 11th order correctors
  if (streq(ri,"sei"))       return REB_INTEGRATOR_SEI;        // 2 - SEI integrator for shearing sheet simulations, symplectic, needs OMEGA variable
  if (streq(ri,"leapfrog"))  return REB_INTEGRATOR_LEAPFROG;   // 4 - LEAPFROG integrator, simple, 2nd order, symplectic
  if (streq(ri,"none"))      return REB_INTEGRATOR_NONE;       // 7 - Do not integrate anything
  if (streq(ri,"janus"))     return REB_INTEGRATOR_JANUS;      // 8 - Bit-wise reversible JANUS integrator.
  if (streq(ri,"mercurius")) return REB_INTEGRATOR_MERCURIUS;  // 9 - MERCURIUS integrator 
  if (streq(ri,"saba"))      return REB_INTEGRATOR_SABA;       // 10 - SABA integrator family (Laskar and Robutel 2001)
  if (streq(ri,"eos"))       return REB_INTEGRATOR_EOS;        // 11 - Embedded Operator Splitting (EOS) integrator family (Rein 2019)
  if (streq(ri,"bs"))        return REB_INTEGRATOR_BS;         // 12 - Gragg-Bulirsch-Stoer 
  if (streq(ri,"whfast512")) return REB_INTEGRATOR_WHFAST512;  // 21 - WHFast integrator, optimized for AVX512
  error("Invalid integrate=%s", ri);
  return -1;
}

int gravity_type(string ti)
{
  if (streq(ti,"none"))         return REB_GRAVITY_NONE;        // 0
  if (streq(ti,"basic"))        return REB_GRAVITY_BASIC;       // 1
  if (streq(ti,"compensated"))  return REB_GRAVITY_COMPENSATED; // 2
  if (streq(ti,"tree"))         return REB_GRAVITY_TREE;        // 3
  error("Invalid gravity=%s", ti);  
  return -1;
}


void nemo_main()
{
    string infile = getparam("in");
    string outfile = getparam("out");
    int port = getiparam("server");
    double eps = getdparam("eps");
    double box = getdparam("box");
    double dt = getdparam("dt");
    double dtout = getdparam("dtout");
    double G = getdparam("G");
    double tolerance = getdparam("tolerance");
    int reb_integrator = integrator_type(getparam("integrator"));
    int reb_gravity = gravity_type(getparam("gravity"));
    double tstop = INFINITY;
    double tsnap;
    int nbody;
    if (hasvalue("tstop")) tstop = getdparam("tstop");      
    
    struct reb_simulationarchive *sa = reb_simulationarchive_create_from_file(infile);
    if (!sa) error("bad simulation archive %s", infile);
    struct reb_simulation *r = reb_simulation_create_from_simulationarchive(sa, -1);
    if (!r) error("cannot find snapshot");
    reb_simulationarchive_free(sa);
    nbody = r->N;
    tsnap = r->t;
    dprintf(0,"Read snapshot N=%d at t=%g\n", nbody, tsnap);
    dprintf(0,"Integrator:     %d\n", reb_integrator);
    dprintf(0,"Gravity Solver: %d\n", reb_gravity);
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:<port>
    if (port > 0)
      reb_simulation_start_server(r, port);

    // Setup constants
    r->G              = G;        
    r->integrator     = reb_integrator;
    r->gravity        = reb_gravity;
    r->dt             = dt;
    r->softening      = eps;
    r->opening_angle2 = tolerance;
    r->heartbeat      = heartbeat;

    if (hasvalue("box"))
      reb_simulation_configure_box(r, box, 1, 1, 1);
    reb_simulation_save_to_file_interval(r, outfile, dtout);
    reb_simulation_synchronize(r);   // note docs are wrong on name
    reb_simulation_integrate(r, tstop);
    //reb_simulation_save_to_file(r, outfile);
    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.0*r->dt)){
        reb_simulation_output_timing(r, 0);
    }
}
