// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// TestGrav.cc                                                                 |
//                                                                             |
// Copyright (C) 2000-2004, 2008 Walter Dehnen                                 |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
#define falcON_RepAction 1
#define falcON_NONEMO
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#define falcON_NONEMO

#include <main.h>
#include <body.h>
#include <forces.h>
#include <numerics.h>

using std::ofstream;
using std::setfill;
using std::setw;
using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;

using namespace falcON;
namespace {
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

  struct myRNG {
    static const double iRMAX;
    myRNG(unsigned seed) { std::srand(seed); }
    double operator() () const { return rand() * iRMAX; }
  };
  const double myRNG::iRMAX = 1./(1.+double(RAND_MAX));
}

void falcON::main(int argc, const char* argv[]) falcON_THROWING
{
  if(argc < 6) {
    cerr<<
      " \"TestGrav MOD gam N S EPS [Ngrow theta K"
#  ifdef falcON_INDI
#  ifdef falcON_ADAP
      " Nsoft Nref emin"
#  endif
      " epar"
#  endif
      " Ncrit"
#if (0)
      " Ncut"
#endif
      " Nact MAC s Ncb0 Ncb1 Ncc Ncs Rmax DUMP]\" with \n"
      " MOD = 0/1/2/3/4     : hom. sphere / Plummer / gamma-model"
      " / Kuzmin disk / hom. disk\n"
      " gam                 : parameter of model (if any)\n"
      " N                   : number of positions\n"
      " S = long            : seed for RNG\n"
#ifdef falcON_INDI
      " EPS                 : fixed / maximum individual softening length\n"
#else
      " EPS                 : softening length\n"
#endif
      " Ngrow(default   0)  : # additional grow()s\n"
      " theta(default "<<std::setw(4)<<Default::theta
	<<") : accuracy parameter\n"
      " K    (default   "<<falcON_KERNEL_TEXT<<")  : P_n softening kernel\n"
#ifdef falcON_INDI
#ifdef falcON_ADAP
      " Nsoft(default   0)  : use individual softening with Nsoft\n"
      " Nref (default Nsoft): use individual softening with Nref\n"
      " emin (default   0)  : lower limit for eps_i\n"
#endif
      " epar (default   0)  : if > 0: mean fixed individual softening length\n"
#endif
      " Ncrit(default "<<std::setw(3)<<Default::Ncrit
	               <<")  : don't split cell with <= Nc bodies\n"
#if (0)
      " Ncut (default   0)  : if != 0: N_cut in re-growing of tree\n"
#endif
      " Nact (default   N)  : # active bodies\n"
      " MAC  (default   "<<falcON_MAC_TEXT
	<<")  : 0/1/2/3: theta=const, f(M), f(M/rmax^2), f(M/rmax)\n"
      " s    (default   1)  : interweave interaction & evaluation phase?\n"
      " Ncb0 (default "<<std::setw(3)<<Default::direct[0]<<")  : N_cb^pre\n"
      " Ncb1 (default "<<std::setw(3)<<Default::direct[1]<<")  : N_cb^post\n"
      " Ncc  (default "<<std::setw(3)<<Default::direct[2]<<")  : N_cc^post\n"
      " Ncs  (default "<<std::setw(3)<<Default::direct[3]<<")  : N_cs\n"
      " Rmax (default 1e3)  : max radius\n"
      " DUMP (default   0)  : [0/1] dump nodes to tree.cells and tree.leafs\n";
    std::exit(1);
  }
  register clock_t  cpu0 = clock(), cpu1;
#ifdef falcON_INDI
#ifdef falcON_ADAP
  real              Nsoft=zero,emin=zero;
  unsigned          NREF;
#endif
  real              EPAR = zero;
  bool              indiv_soft = false;
#endif
  real              rmax = 1.e3;
  int               p=1, MOD, dump=0;
  MAC_type          MAC = Default::mac;
  kern_type         K   = Default::kernel;
  unsigned          Seed;
  unsigned          N, Nbod[BT_NUM]={0}, Nact;
  unsigned          split=1, Ngrow=0;
  int               Ncrit=Default::Ncrit, 
                    DIR[4]={Default::direct[0],Default::direct[1],
			    Default::direct[2],Default::direct[3]};
  real              EPS, theta=Default::theta, MTOT=one;

  MOD  = atoi(argv[p++]);
  real GAM;
  GAM  = atof(argv[p++]);
  N    = atoi(argv[p++]); Nbod[bodytype::std]=N;
  Nact = N;
  Seed = atoi(argv[p++]);
  EPS  = atof(argv[p++]);
  if(argc>p) Ngrow  = atoi(argv[p++]);
  if(argc>p) theta  = atof(argv[p++]);
  if(argc>p) K      = kern_type(atoi(argv[p++]) % 10);
#ifdef falcON_INDI
#ifdef falcON_ADAP
  if(argc>p) Nsoft  = atof(argv[p++]);
  if(Nsoft)  indiv_soft = true;
  NREF = Nsoft>0? (Nsoft>1? falcON::uint(Nsoft+half) : 1) : 0;
  if(argc>p) NREF   = atoi(argv[p++]);
  if(argc>p) emin   = atof(argv[p++]);
#endif
  if(argc>p) EPAR   = atof(argv[p++]);
  if(EPAR>0) indiv_soft = true;
#endif
  if(argc>p) Ncrit  = atoi(argv[p++]);
#if (0)
  if(argc>p) Ncut   = atoi(argv[p++]);
#endif
  if(argc>p) Nact   = atoi(argv[p++]); if(Nact==0 || Nact>N) Nact=N;
  if(argc>p) MAC    = MAC_type (atoi(argv[p++]) % 10);
  if(argc>p) split  = atoi(argv[p++]);
  if(argc>p) DIR[0] = atoi(argv[p++]);
  if(argc>p) DIR[1] = atoi(argv[p++]);
  if(argc>p) DIR[2] = atoi(argv[p++]);
  if(argc>p) DIR[3] = atoi(argv[p++]);
  if(argc>p) rmax   = atoi(argv[p++]);
  if(argc>p) dump   = atoi(argv[p++]);
  bodies __BODIES(Nbod), *BODIES=&__BODIES;
#ifdef falcON_INDI
  if(indiv_soft) BODIES->add_fields(fieldset::e);
#endif
  EPS  = falcON::abs(EPS);
  vect  *A = new vect[N];
  real *RH = new real[N];
  real *RA = new real[N];
  myRNG Rand(Seed);
  const real g1=1.-GAM, g2=2.-GAM, ig2=1./g2, g3=3-GAM, ig3=1./g3;
  real       tm,r=zero,R,phi,fr,P,rh=zero, cth;
  double     sph,cph,rho=0., mi = MTOT/real(N);
  LoopAllBodies(BODIES,Bi) {
    do {
      switch(MOD) {
      case 0:  // homogeneous
	r  = pow(Rand(),1./3.);
	P  = 0.5*(3*r*r-1);
	fr =-1.;
	rh = 0.75/Pi;
	break;
      case 1:  // Plummer
	tm = pow(Rand(),2./3.);
	tm/= 1.-tm;                    // r^2
	r  = sqrt(tm);
	P  =-1./sqrt(tm+1);
	fr = P*P*P;
	rh = 0.75/Pi*pow(-P,5u);
	break;
      case 2:  // gamma-model
	tm = pow(Rand(),double(ig3));
	r  = tm / (1.-tm);
	P  = (GAM==2.)?  log(tm) : -ig2*(1-pow(tm,g2));
	fr =-pow(tm,g1) * square(1-tm) / r;
	rh = 0.25*g3/Pi * pow(1.-tm,4u)/pow(tm,GAM);
	break;
      case 3:
	tm = Rand();
	r  = sqrt((2-tm)*tm)/(1.-tm);
	P  =-1./sqrt(r*r+1);
	fr = P*P*P;
	rh =-0.5/Pi*fr;
	break;
      case 4:
	tm = Rand();
	r  = sqrt(tm);
	P  = 2*(r*r/3.-2);
	fr =-4*r/3.;
	rh = 1.0/Pi;
	break;
      default:
	cerr<<" unknown model\n";
	std::exit(1);
      }
    } while(rmax>zero && r>rmax);
    cth = (MOD>2) ? zero : 2*Rand()-1;
    R   = r*sqrt(1.-cth*cth);
    Bi.pos()[2]   = r * cth;
    phi           = TPi*Rand();
    sph           = sin(phi);
    cph           = cos(phi);
    Bi.pos()[0]   = R * sph;
    Bi.pos()[1]   = R * cph;
    Bi.mass()     = mi;
    Bi.flag_as_active();
    A [bodyindex(Bi)] = MTOT * fr * pos(Bi);   // exact force
    RH[bodyindex(Bi)] = rh;                    // exact density
    RA[bodyindex(Bi)] = r;                     // radius
    rho          += rh;
  }
  rho = N / rho;                           // inverse mean density
#ifdef falcON_INDI
  if(indiv_soft) {
    std::cerr<<" indiv_soft = "<<indiv_soft<<'\n';
    std::cerr<<" epar       = "<<EPAR<<'\n';
    register double emean=0.;
    LoopAllBodies(BODIES,Bi) {
      Bi.eps() = min(EPS, real(EPAR*pow(RH[bodyindex(Bi)]*rho,-1./3.)));
      emean   += eps(Bi);
    }
    emean /= double(N);
    std::cerr<<" mean eps   = "<<emean<<'\n';
  }
#endif
  const bool all = Nact >= N;
  if(!all) {
    int *I = falcON_NEW(int,N);
    HeapIndex(RA,N,I);
    const real RqMax = square(RA[I[Nact-1]]);
    falcON_DEL_A(I);
    LoopAllBodies(BODIES,Bi)
      if( norm(pos(Bi)) > RqMax )
	Bi.unflag_active();
      else
	Bi.flag_as_active();
  }

//   // create file with masses & positions
//   ofstream out("positions.dat");
//   out<<N<<"\n";
//   LoopAllBodies(BODIES,Bi)
//     out<<setprecision(12)<<mi<<" "<<pos(Bi)<<"\n";
//   out.close();
//   //                                    
  cpu1 = clock();
  cout<<"\n"<<setprecision(6)
      <<" time needed for set up of X_i:                  "
      <<(cpu1 - cpu0)/real(CLOCKS_PER_SEC)<<endl;
  cpu0 = cpu1;
#ifdef falcON_INDI
  forces falcon(BODIES,EPS,theta,K,indiv_soft,one,MAC,one,DIR);
#else
  forces falcon(BODIES,EPS,theta,K,one,MAC,DIR);
#endif
  for(register unsigned ig=0; ig<=Ngrow; ++ig) {
#if (0)
    if(Ncut) falcon.re_grow(Ncut,Ncrit);
    else     
#endif
    falcon.grow(Ncrit);
    cpu1 = clock();
    cout<<setprecision(6)
	<<" time needed for forces::grow():                 "
	<<(cpu1 - cpu0) / real(CLOCKS_PER_SEC)<<endl;
    cpu0 = cpu1;
  }

  cpu0 = clock();
#ifdef falcON_ADAP
  falcon.approximate_gravity(split,all,Nsoft,NREF,emin,0.);
#else
  falcon.approximate_gravity(split,all);
#endif
  if(dump) falcon.dump_nodes("tree.cells","tree.leafs");

  cpu1 = clock();
  cout<<setprecision(6)
      <<" time needed for forces::approximate_gravity():  "
      <<(cpu1 - cpu0) / real(CLOCKS_PER_SEC)<<endl;
  falcon.stats(cout);
  cout<<"\n";

// #define WRITE_OUT_ACC 1
#ifdef WRITE_OUT_ACC
  std::ofstream accout("test.acc");
#endif
  register real ase=0.,damx=0.,dacc;
  register vect_d SF(0.);
  register bool pot_pos=0;
  LoopAllBodies(BODIES,Bi) if(is_active(Bi)) {
    vect&Acc(A[bodyindex(Bi)]);
    dacc = norm(Acc - acc(Bi));
    ase += dacc;
    damx = max(damx,dacc);
    SF  += mass(Bi)*acc(Bi);
    if(pot(Bi) > zero) pot_pos = true;
#ifdef WRITE_OUT_ACC
    accout<<setfill(' ')
	  <<setw(17)<<setprecision(12)<<Bi.pos()[0]<<" "
	  <<setw(17)<<setprecision(12)<<Bi.pos()[1]<<" "
	  <<setw(17)<<setprecision(12)<<Bi.pos()[2]<<"  "
	  <<setw(17)<<setprecision(12)<<Bi.pot()   <<"  "
	  <<setw(17)<<setprecision(12)<<Bi.acc()[0]<<" "
	  <<setw(17)<<setprecision(12)<<Bi.acc()[1]<<" "
	  <<setw(17)<<setprecision(12)<<Bi.acc()[2]<<"  "
	  <<setw(17)<<setprecision(12)<<Acc[0]<<" "
	  <<setw(17)<<setprecision(12)<<Acc[1]<<" "
	  <<setw(17)<<setprecision(12)<<Acc[2]<<'\n';
#endif
  }
#ifdef WRITE_OUT_ACC
  accout.close();
#endif
  ase /= double(Nact);
  switch(MOD){
  case 1: ase *= 13.125; break;
  case 3: ase *= 17.5;   break;
  case 4: ase *= 1.125;  break;
  case 0: ase /= 0.6;    break;
  case 2: if(GAM<1.66) ase /= (MTOT*MTOT*ms_Force_gamma(GAM)); break;    
  default: break;
  }
  if(MOD==2 && GAM>=1.66)
    cout<<setprecision(10)<<" ASE(F)           = "<< ase <<"\n";
  else
    cout<<setprecision(10)<<" ASE(F)/<F^2>     = "<< ase <<"\n";
  cout<<" max (dF)^2       = "<< damx<<"\n"
      <<" Sum m_i acc_i    = "<< SF  <<endl;
  if(pot_pos)
    cerr<<" WARNING: some potentials were > 0"<<endl;
//  cout<<endl;
//   for(;;) {
//     cerr<<"i = "; cin>>i;
//     cerr<<" "<<acc(B+i)<<" "<<neighbours(B+i)<<"\n";
//   }
#if(0)
#ifdef falcON_INDI
  if(Nsoft || use_indiv) {
    std::ofstream FILE("test.in");
    if(FILE.is_open()) {
      FILE<<setfill(' ')<<N<<' '<<Nsoft<<"\n";
      LoopAllBodies(BODIES,Bi)
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
      LoopAllBodies(BODIES,Bi)
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
#ifdef falcON_INDI
  }
#endif
#endif

}
