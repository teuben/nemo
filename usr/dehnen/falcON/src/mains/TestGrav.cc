// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// TestGrav.cc                                                                 |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <body.h>
#include <falcON.h>
#include <public/Pi.h>
#include <public/exit.h>

using std::ofstream;
using std::setfill;
using std::setw;
using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;

using namespace nbdy;

double ms_Force_gamma(const double gamma)
{
  register double x=5-3*gamma,f;
  if(x < 0.) return 0.;
  f =1/x;
  f-=4/(x+=1);
  f+=6/(x+=1);
  f-=4/(x+=1);
  f+=1/(x+=1);
  return (3-gamma)*f;
}

double ms_density_gamma(const double gamma)
{
  register double x=3-3*gamma,f;
  if(x < 0.) return 0.;
  f = 1/x;
  f-= 8/(x+=1);
  f+=28/(x+=1);
  f-=56/(x+=1);
  f+=70/(x+=1);
  f-=56/(x+=1);
  f+=28/(x+=1);
  f-= 8/(x+=1);
  f+= 1/(x+=1);
  return cube(3-gamma)*f/square(4*Pi);
}

int main(int argc, char* argv[])
{
#ifdef TWODIMENSIONAL
  if(argc < 5) {
    cerr<<
      " \"TestGrav MOD N S EPS [THE K Ns"
#ifdef ALLOW_INDI
      " Nsoft Nref"
#endif
      " Nc MAC s Ng Ncb0 Ncb1 Ncc Ncs DUMP]\" with \n"
      " MOD = 3/4           : Kuzmin disk / homogeneous disk\n"
#else
  if(argc < 6) {
    cerr<<
      " \"TestGrav MOD gam N S EPS [Ng THE K Ns"
#ifdef ALLOW_INDI
      " Nsoft Nref"
#endif
      " Nc MAC s Ncb0 Ncb1 Ncc Ncs DUMP]\" with \n"
      " MOD = 0/1/2/3/4     : hom. sphere / Plummer / gamma-model"
      " / Kuzmin disk / hom. disk\n"
      " gam                 : parameter of model (if any)\n"
#endif
      " N                   : number of positions\n"
      " S = long            : seed for RNG\n"
#ifdef ALLOW_INDI
      " EPS                 : fixed / maximum / individual softening length\n"
#else
      " EPS                 : softening length\n"
#endif
      " Ng   (default   0)  : # additional grow()s\n"
      " THE  (default "<<std::setw(3)<<Default::theta
	<<")  : accuracy parameter\n"
      " K    (default   "<<DEFAULT_KERNEL<<")  : P_n softening kernel\n"
      " Ns   (default   N)  : # sinks\n"
#ifdef ALLOW_INDI
      " Nsoft(default   0)  : use individual softening with Nsoft\n"
      " Nref (default Nsoft): use individual softening with Nref\n"
#endif
      " Nc   (default "<<std::setw(3)<<Default::Ncrit
	               <<")  : don't split cell with <= Nc bodies\n"
      " MAC  (default   "<<DEFAULT_MAC
	<<")  : 0/1/2/3: theta=const, f(M), f(M/rmax^2), f(M/rmax)\n"
      " s    (default   1)  : interweave interaction & evaluation phase?\n"
      " Ncb0 (default "<<std::setw(3)<<Default::direct[0]<<")  : N_cb^pre\n"
      " Ncb1 (default "<<std::setw(3)<<Default::direct[1]<<")  : N_cb^post\n"
      " Ncc  (default "<<std::setw(3)<<Default::direct[2]<<")  : N_cc^post\n"
      " Ncs  (default "<<std::setw(3)<<Default::direct[3]<<")  : N_cs\n"
      " DUMP (default   0)  : [0/1] dump nodes to tree.cells and tree.souls\n";
    nbdy::exit(1);
  }

  register clock_t  cpu0 = clock(), cpu1;
#ifdef ALLOW_INDI
  real              NSOFT=zero;
  bool              use_indiv = false;
  soft_type         Soft;
  unsigned          NREF;
#endif
  const real        rmax = 1.e3;
  int               p=1, MOD, T=0;
  MAC_type          MAC = Default::mac;
  kern_type         K   = Default::kernel;
  long              Seed;
  unsigned          N, Ns, split=1, Ng=0;
  int               Ncrit=Default::Ncrit,
                    DIR[4]={Default::direct[0],Default::direct[1],
			    Default::direct[2],Default::direct[3]};
  real              EPS, TH=Default::theta, MTOT=one;

  MOD  = atoi(argv[p++]);
#ifndef TWODIMENSIONAL
  real GAM;
  GAM  = atof(argv[p++]);
#endif
  N    = atoi(argv[p++]); Ns=N;
  Seed = atoi(argv[p++]);
  EPS  = atof(argv[p++]);
  if(argc>p) Ng     = atoi(argv[p++]);
  if(argc>p) TH     = atof(argv[p++]);
  if(argc>p) K      = kern_type(atoi(argv[p++]) % 10);
  if(argc>p) Ns     = atoi(argv[p++]);
#ifdef ALLOW_INDI
  if(EPS < 0.) use_indiv = true;
  if(argc>p) NSOFT  = atof(argv[p++]);
  NREF = NSOFT>0? (NSOFT>1? nbdy::uint(NSOFT+half) : 1) : 0;
  if(argc>p) NREF   = atoi(argv[p++]);
#endif
  if(argc>p) Ncrit  = atoi(argv[p++]);
  if(argc>p) MAC  = MAC_type (atoi(argv[p++]) % 10);
  if(argc>p) split  = atoi(argv[p++]);
  if(argc>p) DIR[0] = atoi(argv[p++]);
  if(argc>p) DIR[1] = atoi(argv[p++]);
  if(argc>p) DIR[2] = atoi(argv[p++]);
  if(argc>p) DIR[3] = atoi(argv[p++]);
#ifndef TWODIMENSIONAL
#endif
  if(argc>p) T     = atoi(argv[p++]);
  bodies *BODIES   = 
#ifdef ALLOW_INDI
    (NSOFT || use_indiv)? 
    new bodies(N,bodies::DEFBITS() | io::e) :
#endif
    new bodies(N);
#ifdef ALLOW_INDI
  Soft = (NSOFT||use_indiv)? nbdy::individual : nbdy::global;
#endif
  EPS  = nbdy::abs(EPS);
  vect  *A = new vect[N];
  real *RH = new real[N];
  srand48(Seed);

#ifndef TWODIMENSIONAL
  const    real g1=1.-GAM, g2=2.-GAM, ig2=1./g2, g3=3-GAM, ig3=1./g3;
#endif
  register real tm,r=zero,R,phi,fr,P,rh=zero;
#ifndef TWODIMENSIONAL
  register real cth;
#endif
  register double sph,cph;
  register double mi = MTOT/real(N);
  LoopBodies(bodies,BODIES,Bi) {
    do {
      switch(MOD) {
#ifndef TWODIMENSIONAL
      case 0:  // homogeneous
	r  = pow(drand48(),1./3.);
	P  = 0.5*(3*r*r-1);
	fr =-1.;
	rh = 0.75/Pi;
	break;
      case 1:  // Plummer
	tm = pow(drand48(),2./3.);
	tm/= 1.-tm;                    // r^2
	r  = sqrt(tm);
	P  =-1./sqrt(tm+1);
	fr = P*P*P;
	rh = 0.75/Pi*pow(-P,5u);
	break;
      case 2:  // gamma-model
	tm = pow(drand48(),double(ig3));
	r  = tm / (1.-tm);
	P  = (GAM==2.)?  log(tm) : -ig2*(1-pow(tm,g2));
	fr =-pow(tm,g1) * square(1-tm) / r;
	rh = 0.25*g3/Pi * pow(1.-tm,4u)/pow(tm,GAM);
	break;
#endif
      case 3:
	tm = drand48();
	r  = sqrt((2-tm)*tm)/(1.-tm);
	P  =-1./sqrt(r*r+1);
	fr = P*P*P;
	rh =-0.5/Pi*fr;
	break;
      case 4:
	tm = drand48();
	r  = sqrt(tm);
	P  = 2*(r*r/3.-2);
	fr =-4*r/3.;
	rh = 1.0/Pi;
	break;
       default:
	cerr<<" unknown model\n";
	nbdy::exit(1);
      }
    } while(r>rmax);
#ifdef TWODIMENSIONAL
    R = r;
#else
    cth = (MOD==3) ? zero : 2*drand48()-1;
    R   = r*sqrt(1.-cth*cth);
    Bi.pos(2) = r * cth;
#endif
    phi       = TPi*drand48();
    sph       = sin(phi);
    cph       = cos(phi);
    Bi.pos(0) = R * sph;
    Bi.pos(1) = R * cph;
    Bi.mass() = mi;
    if(index(Bi) < Ns) Bi.flag_as_sink();
    else               Bi.unflag_sink();
#ifdef ALLOW_INDI
    if(use_indiv) Bi.eps() = EPS;
#endif
    A [index(Bi)] = MTOT * fr * pos(Bi);   // exact force
    RH[index(Bi)] = rh;             // exact density
  }
//   // create file with masses & positions
//   ofstream out("positions.dat");
//   out<<N<<"\n";
//   for(Bi=B,i=0; i!=N; ++Bi,++i)
//     out<<setprecision(12)<<mi<<" "<<pos(Bi)<<"\n";
//   out.close();
//   //                                    
  cpu1 = clock();
  cout<<"\n"<<setprecision(6)
      <<" time needed for set up of X_i:           "
      <<(cpu1 - cpu0)/real(CLOCKS_PER_SEC)<<endl;
  cpu0 = cpu1;
#ifdef ALLOW_INDI
  falcON falcon(BODIES,EPS,TH,K,Soft,MAC);
#else
  falcON falcon(BODIES,EPS,TH,K,MAC);
#endif
  for(register unsigned ig=0; ig<=Ng; ++ig) {
    falcon.grow(Ncrit);
    cpu1 = clock();
    cout<<setprecision(6)
	<<" time needed for falcON::grow():          "
	<<(cpu1 - cpu0) / real(CLOCKS_PER_SEC)<<endl;
    cpu0 = cpu1;
  }
  if(T) falcon.dump_nodes("tree.cells","tree.souls");

  cpu0 = clock();
#ifdef ALLOW_INDI
  falcon.approximate_gravity(split,NSOFT,NREF,0.,DIR);
#else
  falcon.approximate_gravity(split,DIR);
#endif
  cpu1 = clock();
  cout<<setprecision(6)
      <<" time needed for falcON::approximate():   "
      <<(cpu1 - cpu0) / real(CLOCKS_PER_SEC)<<endl;
  falcon.stats(cout);
  cout<<"\n";

  register real ase=0.,damx=0.,dacc;
  register vect SF=zero;
  LoopBodies(bodies,BODIES,Bi) {
    dacc = norm(A[index(Bi)] - acc(Bi));
    ase += dacc;
    damx = max(damx,dacc);
    SF  += mass(Bi)*acc(Bi);
//     cerr<<setfill(' ')
// 	<<setw(13)<<setprecision(8)<<Bi.pos()[0]<<" "
// 	<<setw(13)<<setprecision(8)<<Bi.pos()[1]<<" "
// 	<<setw(13)<<setprecision(8)<<Bi.pos()[2]<<" "
// 	<<setw(13)<<setprecision(8)<<Bi.acc()[0]<<" "
// 	<<setw(13)<<setprecision(8)<<Bi.acc()[1]<<" "
// 	<<setw(13)<<setprecision(8)<<Bi.acc()[2]<<endl;
  }
  ase /= double(N);
  switch(MOD){
  case 1: ase *= 13.125; break;
  case 3: ase *= 17.5;   break;
  case 4: ase *= 1.125;  break;
  case 0: ase /= 0.6;    break;
#ifndef TWODIMENSIONAL
  case 2: if(GAM<1.66) ase /= (MTOT*MTOT*ms_Force_gamma(GAM)); break;    
#endif      
  default: break;
  }
#ifndef TWODIMENSIONAL
  if(MOD==2 && GAM>=1.66)
    cout<<setprecision(10)<<" ASE(F)           = "<< ase <<"\n";
  else
#endif      
    cout<<setprecision(10)<<" ASE(F)/<F^2>     = "<< ase <<"\n";
  cout<<" max (dF)^2       = "<< damx<<"\n"
      <<" Sum m_i acc_i    = "<< SF  <<endl;
//  cout<<endl;
//   for(;;) {
//     cerr<<"i = "; cin>>i;
//     cerr<<" "<<acc(B+i)<<" "<<neighbours(B+i)<<"\n";
//   }
#if(0)
#ifdef ALLOW_INDI
  if(NSOFT || use_indiv) {
    std::ofstream FILE("test.in");
    if(FILE.is_open()) {
      FILE<<setfill(' ')<<N<<' '<<NSOFT<<"\n";
      LoopBodies(bodies,BODIES,Bi)
	FILE
	  <<setw(13)<<setprecision(8)<<Bi.pos()[0]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.pos()[1]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.pos()[2]<<"  "
	  <<setw(13)<<setprecision(8)<<Bi.pot()   <<"  "
	  <<setw(13)<<setprecision(8)<<Bi.acc()[0]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.acc()[1]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.acc()[2]<<"  "
	  <<setw(13)<<setprecision(8)<<Bi.eps()<<"\n";
      FILE.close();
    }
  } else {
#endif
    std::ofstream FILE("test.gl");
    if(FILE.is_open()) {
      FILE<<setfill(' ')<<N<<' '<<EPS<<"\n";
      LoopBodies(bodies,BODIES,Bi)
	FILE
	  <<setw(13)<<setprecision(8)<<Bi.pos()[0]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.pos()[1]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.pos()[2]<<"  "
	  <<setw(13)<<setprecision(8)<<Bi.pot()   <<"  "
	  <<setw(13)<<setprecision(8)<<Bi.acc()[0]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.acc()[1]<<" "
	  <<setw(13)<<setprecision(8)<<Bi.acc()[2]<<"\n";
      FILE.close();
    }
#ifdef ALLOW_INDI
  }
#endif
#endif
  delete BODIES;

}
