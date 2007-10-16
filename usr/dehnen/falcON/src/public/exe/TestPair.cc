#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

#define falcON_NONEMO

#include <main.h>
#include <body.h>
#include <forces.h>

using std::cerr;
using std::cout;
using std::endl;
using std::setw;
using std::ofstream;
using std::setprecision;

using namespace falcON;

const double twothird   = 2./3.;
const double Pi         = 3.14159265358979323846264338328;
const double Pih        = 0.5 * Pi;
const double TPi        = 2.0 * Pi;

void dmintmin(const bodies *BODS,
	      bodies::index A,
	      bodies::index B,
	      falcON::real  T,
	      falcON::real &dmin,
	      falcON::real &tmin)
{
  if(T==zero) {
    tmin = zero;
    dmin = sqrt(dist_sq(BODS->pos(A),BODS->pos(B)));
  } else {
    register vect R = BODS->pos(A)-BODS->pos(B), V=BODS->vel(A)-BODS->vel(B);
    tmin = max(zero, min(T, -(R*V)/norm(V)));
    dmin = sqrt(norm(R+tmin*V));
  }
}

void falcON::main(int argc, char* argv[]) falcON_THROWING
{
  if(argc < 6 || argc > 10) {
    cerr<<"\n Testing tree sticky/SPH particle support with a Plummer sphere\n\n"
	<<" usage: \"TestPair N S Ns s tau [M C Nc T]\"\n with \n"
	<<" N               : number of particles\n"
	<<" S = long int    : seed for RNG\n"
	<<" Ns              : first Nsub bodies to make sticky\n"
        <<" s               : size of sticky/SPH particles\n"
        <<" tau             : time step; if < 0: SPH\n"
	<<" M  (default 1)  : 1/0: max/sum (h_i,h_j)\n"
	<<" C  (default 0)  : 1/0: count partners?\n"
	<<" Nc (default "<<std::setw(2)<<Default::Ncrit
	<<") :  don't split cell with <=Ncrit bodies\n"
        <<" T  (default  0) : dump nodes to tree.cells and tree.leafs\n";
      
    exit(1);
  }

  register clock_t      cpu0 = clock(), cpu1;
  register real         time_prep, time_tree, time_subt;
  int                   p=1, T=0, Mx=1, Count=0;
  long                  Seed;
  unsigned              N, NS, Ncrit=Default::Ncrit;
  real                  S,TAU;

  N     = atoi(argv[p++]);
  Seed  = atoi(argv[p++]);
  NS    = atoi(argv[p++]); if(NS>N) NS=N;
  S     = atof(argv[p++]);
  TAU   = atof(argv[p++]);
  if(argc>p) Mx    = atoi(argv[p++]);
  if(argc>p) Count = atoi(argv[p++]);
  if(argc>p) Ncrit = atoi(argv[p++]);
  if(argc>p) T     = atoi(argv[p++]);

  srand48(Seed);

  bodies         *BODIES = new bodies(N,fieldset(bodies::DefBits),NS);
  const    real   mass=1./real(N);
  register double M,r,rq,P,Pot,veq,f0,v,vq,fs,R,cth,phi;
  LoopAllBodies(BODIES,Bi) {
    M   = pow(drand48(),twothird);
    rq  = M/(1-M);
    Pot =-1./sqrt(rq+1.);
    r   = sqrt(rq);
    veq =-2*Pot;
    do {
      P   = veq*pow(drand48(),1./4.5);
      vq  = veq-P;
      v   = sqrt(vq);
      f0  = v*pow(P,3.5);
      fs  = veq*f0;
      f0 *= v;
    } while (fs* drand48()>f0);
    cth = 2*drand48()-1.;
    phi = TPi*drand48();
    R   = r * sqrt(1.-cth*cth);
    Bi.pos()[0] =  R * sin(phi);
    Bi.pos()[1] =  R * cos(phi);
    Bi.pos()[2] =  r * cth;
    cth = 2*drand48()-1.;
    phi = TPi*drand48();
    R   = v * sqrt(1.-cth*cth);
    Bi.vel()[0] =  R * sin(phi);
    Bi.vel()[1] =  R * cos(phi);
    Bi.vel()[2] =  v * cth;
    Bi.mass() = mass;
    Bi.flag_as_active();
    if(bodyindex(Bi)<NS) {
      Bi.size() = S;
      Bi.flag_as_sticky();
      Bi.flag_as_active();
    }
    else {
      Bi.unflag_sticky();
      Bi.unflag_active();
    }
  }
  cpu1 = clock();
  time_prep = (cpu1 - cpu0) / real(CLOCKS_PER_SEC);
  cpu0 = cpu1;
  cerr<<" setup of (x,v) used "<<time_prep<<" seconds \n";

  forces FALCON(BODIES,0.);
  FALCON.grow(Ncrit);
  cpu1 = clock();
  time_tree = (cpu1 - cpu0) / real(CLOCKS_PER_SEC);
  cpu0 = cpu1;

  cerr<<" tree.grow() used   "<<time_tree<<" seconds \n";
  if(T) FALCON.dump_nodes("tree.cells","tree.leafs");

  unsigned Na=0;
  forces::indx_pair* BP = new forces::indx_pair[NS];
  cpu0 = clock();
  FALCON.make_iaction_list(BP,NS,Na,Mx,TAU,Count);
  cpu1 = clock();
  time_subt = (cpu1 - cpu0) / real(CLOCKS_PER_SEC);
  cerr<<" tree.make_iaction_list() used "<<time_subt<<" seconds \n";

  register unsigned i;
  if(TAU<zero) {
    ofstream out("pairs.c++");
    for(i=0; i!=min(NS,Na); ++i) {
      out<<setw(6) <<BODIES->bodyindex(BP[i][0])<<" "
	 <<setw(6) <<BODIES->bodyindex(BP[i][1])<<"  "
	 <<setw(10)<<setprecision(4)
	 <<dist(BODIES->pos(BP[i][0]),BODIES->pos(BP[i][1]))<<endl;
    }
  } else {
    ofstream out("pairs.c++");
    real dmin=0,tmin=0;
    for(i=0; i!=min(NS,Na); ++i) {
      dmintmin(BODIES,BP[i][0],BP[i][1],TAU,dmin,tmin);
      out<<setw(6) <<BODIES->bodyindex(BP[i][0])<<" "
	 <<setw(6) <<BODIES->bodyindex(BP[i][1])<<"  "
	 <<setw(12)<<setprecision(6)<<tmin<<"  "
	 <<setw(12)<<setprecision(4)<<dmin<<endl;
    }
  }

  if(Na>NS)
    cerr<<" "<<Na<<" pairs found. First "<<NS
	<<" written to file \"pairs.c++\".\n";
  else
    cerr<<" "<<Na<<" pairs found & written to file \"pairs.c++\".\n";

//   for(Bi=B,i=0; i<N; Bi++,i++)
//     if(i>=NS) Bi.flag_to_be_ignored();
//   falcON TREE2(B,N,0.);
//   TREE2.grow(porigin);
//   TREE2.dump_nodes("tre2.cells","tre2.leafs");

  delete BODIES;
}
