// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "cfalcon"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <ctime>
#include <body.h>                                  // the bodies                
#include <public/neighbours.h>                     // finding neighbours        
#include <public/bodyfunc.h>                       // body functions            
#include <externacc.h>                             // external potential        
#include <forces.h>                                // falcON       
#include <assert.h>
#include <math.h>
#include "cfalcon.h"
using namespace jclut;
using namespace falcON;
using falcON::real;
cfalcon::cfalcon()
{
}
extern "C" {
  // load protocol
  int falcon_gravity_(const int * nbody,
                       const float * pos, const float * mass, 
                       float * acc, float * phi,
                       const float * eps,
                       const float * G,
                       const float * theta,
                       const int * kernel_type,
                       const int * ncrit) {
    cfalcon::addGravity(*nbody,pos,mass,acc,phi,*eps,
                        *G,*theta,*kernel_type,*ncrit);
    return 1;
  }
  int falcon_gravity2_(const int * nbody,
                        const float * pos, const float * mass, 
                        const int * nbody_tp,const float * pos_tp,
                        float * acc, float * phi,
                        const bool *  self_p,
                        const float * eps,
                        const float * G,
                        const float * theta,
                        const int * kernel_type,
                        const int * ncrit) {
    std::cerr << "nbody="<< *nbody << " nbody_tp=" << *nbody_tp << " self_p=" << *self_p <<" eps="<<*eps<<"\n";
    cfalcon::addGravity2(*nbody,pos,mass,*nbody_tp,pos_tp,acc,phi,*self_p,*eps,
                        *G,*theta,*kernel_type,*ncrit);
    return 1;
  }
}
bool cfalcon::addGravity(const int nbody,
                         const float * pos, const float * mass, 
                         float * acc, float * phi,
                         const float eps,
                         const float G,
                         const float theta,
                         const int kernel_type,
                         const int ncrit
                         )
{
  bool status = true;
  
  unsigned int Nbod[bodytype::NUM]={0};
  Nbod[falcON::bodytype::std]=nbody;
  falcON::fieldset SRCE(fieldset::m | fieldset::x);
  vect X0(0.);
  vect * RC(&X0);
  bool SOFT = (eps<0);
  if(SOFT) SRCE |= fieldset::e;
  //snapshot *  SHOT = new falcON::snapshot(0.,Nbod,SRCE|fieldset::a|fieldset::p);
  
  int cpt=0;
  unsigned DIR[4]={Default::direct[0],Default::direct[1],
                   Default::direct[2],Default::direct[3]};
  bodies __BODIES(Nbod), *BODIES=&__BODIES;
  
   LoopAllBodies(BODIES,B) {
    B.pos()[0] =  pos[cpt*3+0];
    B.pos()[1] =  pos[cpt*3+1];
    B.pos()[2] =  pos[cpt*3+2];    
    B.mass()   =  mass[cpt];
    B.flag_as_active();
    //B.eps() = eps;
    cpt++;
  }

  MAC_type          MAC = Default::mac;
  forces FALCON(BODIES,eps,theta,kern_type(kernel_type),SOFT,G,MAC,eps,one,DIR);

  if (FALCON.NewtonsG() != zero) {
    FALCON.grow(ncrit, RC);
    FALCON.approximate_gravity(1);
    
    // get back acceleration and potential
    cpt=0;
    LoopAllBodies(BODIES,B) {
      acc[cpt*3+0] = B.acc()[0];
      acc[cpt*3+1] = B.acc()[1];
      acc[cpt*3+2] = B.acc()[2];
      phi[cpt]     = B.pot();
      cpt++;
    }
  }
  return status;
}
// -------------------------------------------------------------------------
bool cfalcon::addGravity2(const int nbody,
                         const float * pos, const float * mass, 
                         const int nbody_tp,
                         // test particles
                         const float * pos_tp,
                         float * acc, float * phi,
                         const bool self_p,
                         const float eps,
                         const float G,
                         const float theta,
                         const int kernel_type,
                         const int ncrit
                         )
{
  bool status = true;
  
  unsigned int Nbod[bodytype::NUM]={0};
  //Nbod[falcON::bodytype::std]=nbody;
  falcON::fieldset SRCE(fieldset::m | fieldset::x);
  vect X0(0.);
  vect * RC(&X0);
  bool SOFT = (eps<0);
  if(SOFT) SRCE |= fieldset::e;
  int N=0;
  unsigned DIR[4]={Default::direct[0],Default::direct[1],
                   Default::direct[2],Default::direct[3]};
  float sink_mass,min_mass;
  std::cerr << "nbody="<< nbody << " nbody_tp=" << nbody_tp << " self_p=" << self_p <<" eps="<<eps<<"\n";
  
  if (self_p) {  // we ASSUME that particles and Test Particles are the same
    if (nbody != nbody_tp ) {
      std::cerr <<
	"\nSELF Potential flag enable, so we suppose that SRC particles and" <<
	" SINK particles are the same,\nBUT nbody=["<< nbody << "] <> nbody_tp=[" <<nbody_tp <<
	"],program aborted...\n\n";
      std::exit(1);
    } 
    else {
      N = nbody;
    }
  } // if (self_p) ...
  else { // No self Potential
    N = (nbody)+(nbody_tp);
  }
  std::cerr <<"NBODY tot =" << N << "\n";
  Nbod[falcON::bodytype::std]=N;
  
  bodies __BODIES(Nbod), *BODIES=&__BODIES;
  
  // Initialize flags [0]=src [1]=sink
  min_mass=mass[0];
  int i=0;
  int j=0;
  
   LoopAllBodies(BODIES,B) {
     
     // proceed ALL/SRC particles
     if (i<nbody) {
       if ( !self_p)         // it's SRC part     
         B.unflag_active();  // SRC are not active
       else                   // it's SINK part    
         B.flag_as_active(); // SINK are active   
       if (mass[i] < min_mass) 
         min_mass = mass[i];
     
       B.mass()   = mass[i];
     
       B.pos()[0] =  pos[i*3+0];
       B.pos()[1] =  pos[i*3+1];
       B.pos()[2] =  pos[i*3+2];  
     }
     else { 
       // proceed SINK particles (without masses)
       if (! self_p) {
         assert(j<nbody_tp);
         //sink_mass = min_mass *  1.e-10/(*nbody_tp);
         sink_mass = min_mass / (nbody_tp);
         //sink_mass = min_mass / (*nbody);
         B.flag_as_active();
         B.mass()   = sink_mass;
         B.pos()[0] = pos_tp[j*3+0];
         B.pos()[1] = pos_tp[j*3+1];
         B.pos()[2] = pos_tp[j*3+2];
         j++; // next j index
       } 
     }
     i++; // next indexun
  }

  MAC_type          MAC = Default::mac;
  forces FALCON(BODIES,eps,theta,kern_type(kernel_type),SOFT,G,MAC,eps,one,DIR);

  if (FALCON.NewtonsG() != zero) {
    FALCON.grow(ncrit, RC);
    FALCON.approximate_gravity(1);
    
    // get back acceleration and potential
    i=0;
    LoopAllBodies(BODIES,B) {
      if (is_active(B)) {
        acc[i*3+0] = B.acc()[0];
        acc[i*3+1] = B.acc()[1];
        acc[i*3+2] = B.acc()[2];
        phi[i]     = B.pot();
        i++;
        assert(i<=nbody_tp);
      }
    }
  }
  return status;
}
