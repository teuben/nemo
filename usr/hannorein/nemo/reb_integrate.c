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
    "server=1234\n                Use web based server on this port (if > 0)",
    "eps=0.01\n                   Softening in terms of r0",
    "dt=2e-5\n                    Integration step, in terms of t0",
    "dtout=\n                     Optional interval to save simulation",
    "tstop=\n                     Set, or let go forever",
    "box=20\n                     Simulation box in terms or r0",
    "integrator=leapfrog\n        rebound integrator",
    "VERSION=0.1\n                16-jan-2024 PJT",
    NULL,
};

string usage="rebound+NEMO nbody integrator (toy model)";

void heartbeat(struct reb_simulation* r);
    
int integrator(string ri)
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
  error("Invalid integrator %s", ri);
  return -1;
}


void nemo_main()
{
    string infile = getparam("in");
    int port = getiparam("server");
    int nbody = 100;
    double eps = getdparam("eps");
    double box = getdparam("box");
    double dt = getdparam("dt");
    int reb_integrator = integrator(getparam("integrator"));
    double tstop = INFINITY;
    if (hasvalue("tstop")) tstop = getdparam("tstop");      
    dprintf(0,"Integrating with %d\n",reb_integrator);
    
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    if (port > 0)
      reb_simulation_start_server(r, 1234);

    if (hasvalue("in"))
      error("NEMO snapshot input from %s not yet supported",infile);

    // Setup system characteristics : a Plummer sphere in nbody (virial) units
    int _N = nbody;                 // Number of particles
    double G = 1;                   // Gravitational constant
    double M = 1;                   // Total mass of the cluster
    double R = 1;                   // Radius of the cluster
    double E = 3./64.*M_PI*M*M/R;   // Energy of the cluster
    double r0 = 16./(3.*M_PI)*R;    // Characteristic length scale
    double t0 = r->G*pow(M,5./2.)*pow(4.*E,-3./2.)*(double)_N/log(0.4*(double)_N); // Relaxation time
    printf("Characteristic size:              %f\n", r0);
    printf("Characteristic time (relaxation): %f\n", t0);

    // Setup constants
    r->G             = G;        
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->dt            = dt*t0;       // timestep
    r->softening     = eps*r0;      // Softening parameter
    r->heartbeat     = heartbeat;
    
    reb_simulation_configure_box(r, box*r0, 1, 1, 1);
    reb_simulation_add_plummer(r, _N, M, R);             // Adds particles
    reb_simulation_move_to_com(r);                       // force center of mass at 0
    if (hasvalue("dtout"))
      reb_simulation_save_to_file_interval(r, getparam("out"), getdparam("dtout"));
    reb_simulation_integrate(r, tstop);
    if (hasvalue("out"))
      reb_simulation_save_to_file(r, getparam("out"));
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.0*r->dt)){
        reb_simulation_output_timing(r, 0);
    }
}
