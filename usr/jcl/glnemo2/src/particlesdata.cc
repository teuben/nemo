// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "particlesdata.h"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include "globaloptions.h"
namespace glnemo {

// ============================================================================
// Constructor                                                                 
ParticlesData::ParticlesData(const ALLOC _model):cmodel(_model)
{
  pos      = NULL;
  vel      = NULL;
  vel_norm = NULL;
#if 1
  rneib    = NULL;
  rho      = NULL;
  temp     = NULL;
#endif
  timu     = NULL;
  nbody    = NULL;
  nemobits = NULL;
  valid_rho=false;
  allocVar();
  coo_max[0] = coo_max[1] = coo_max[2] = 0.0;
  i_max[0]   = i_max[1]   = i_max[2]   = 0;
}
// ============================================================================
// copy Constructor                                                            
const ParticlesData& ParticlesData::operator=(const ParticlesData& m)
{
  if (m.nbody) {
    if (cmodel)
      nbody = (int *) mallocate((char *) nbody, sizeof(int), true);
    else {
      if (nbody) delete nbody;
      nbody = new int;
    }
    *nbody = *m.nbody;
    // positions "pos"
    if (m.pos) {
      if (cmodel)
        pos = (float *) mallocate((char *) pos, sizeof (int) * 3 * (*nbody), true);
      else {
        if (pos) delete [] pos;
        pos = new float[3 * (*nbody)];
      }
      //memcpy((float *) pos, (float *) m.pos, sizeof(float)* (*nbody) * 3);
      for (int i=0; i < (*nbody); i++) {
        pos[i*3]    = m.pos[i*3];
        pos[i*3+1]  = m.pos[i*3+1];
        pos[i*3+2]  = m.pos[i*3+2];
      }
    }
    // velocities "vel"
    if (m.vel) {
      if (cmodel)
        vel = (float *) mallocate((char *) vel, sizeof (float) * 3 * (*nbody), true);
      else {
        if (vel) delete [] vel;
        vel = new float[3 * (*nbody)];
      }
      //memcpy((float *) vel, (float *) m.vel, sizeof(float)* (*nbody) * 3);
      for (int i=0; i < (*nbody); i++) {
        vel[i*3]    = m.vel[i*3];
        vel[i*3+1]  = m.vel[i*3+1];
        vel[i*3+2]  = m.vel[i*3+2];
      }
    }
    else {
      if (vel) {
        if (cmodel) { free ((float *) vel); }
        else        { delete [] vel;        }
      }
      vel = NULL;
    }
    // velocity norm "vel_norm"
    if (m.vel_norm) {
      max_vel_norm = m.max_vel_norm;
      if (cmodel)
        vel_norm = (float *) mallocate((char *) vel_norm, sizeof (float) * (*nbody), true);
      else {
        if (vel_norm) delete [] vel_norm;
        vel_norm = new float[*nbody];
      }
      //memcpy((float *) vel_norm, (float *) m.vel_norm, sizeof(float)* (*nbody));
      for (int i=0; i < (*nbody); i++) {
        vel_norm[i] = m.vel_norm[i];
      }
    }
    else {
      if (vel_norm) {
        if (cmodel) { free ((float *) vel_norm); }
        else        { delete [] vel_norm;        }
      }
      vel_norm = NULL;
    }
#if 1
    // Density "rho"
    if (m.rho) {
      max_rho = m.max_rho;
      min_rho = m.min_rho;
      if (cmodel)
	rho = (float *) mallocate((char *) rho, sizeof (float) * (*nbody), true);
      else {
	if (rho) delete [] rho;
	rho = new float[*nbody];
      }
      //memcpy((float *) rho, (float *) m.rho, sizeof(float)* (*nbody));
      for (int i=0; i < (*nbody); i++) {
	rho[i] = m.rho[i];
      }
      for (int i=0; i < 100; i++) {
        density_histo[i] =m.density_histo[i];
      }
    }
    else {
      if (rho) {
	if (cmodel) { free ((float *) rho); }
        else        { delete [] rho;        }
      }
      rho = NULL;
    }
    // neibourgh radius "rneib"
    if (m.rneib) {
      if (cmodel)
	rneib = (float *) mallocate((char *) rneib, sizeof (float) * (*nbody), true);
      else {
	if (rneib) delete [] rneib;
	rneib = new float[*nbody];
      }
      //memcpy((float *) rneib, (float *) m.rneib, sizeof(float)* (*nbody));
      for (int i=0; i < (*nbody); i++) {
	rneib[i] = m.rneib[i];
      }
    }
    else {
      if (rneib) {
	if (cmodel) { free ((float *) rneib); }
        else        { delete [] rneib;        }
      }
      rneib = NULL;
    }
    // Temperature "temp"
    if (m.temp) {
      max_temp = m.max_temp;
      min_temp = m.min_temp;
      if (cmodel)
        temp = (float *) mallocate((char *) temp, sizeof (float) * (*nbody), true);
      else {
        if (temp) delete [] temp;
        temp = new float[*nbody];
      }
      //memcpy((float *) temp, (float *) m.temp, sizeof(float)* (*nbody));
      for (int i=0; i < (*nbody); i++) {
        temp[i] = m.temp[i];
      }
      for (int i=0; i < 100; i++) {
        density_histo[i] =m.density_histo[i];
      }
    }
    else {
      if (temp) {
        if (cmodel) { free ((float *) temp); }
        else        { delete [] temp;        }
      }
      temp = NULL;
    }

#endif
    valid_temp = m.valid_temp;
    valid_rho = m.valid_rho;
    //memcpy((float *) coo_max, (float *) m.coo_max, sizeof(float)* 3);
    //memcpy((int *) i_max, (int *) m.i_max, sizeof(int)* 3);
    for (int i=0; i <3; i++) {
      coo_max[i] = m.coo_max[i];
      i_max[i]   = m.i_max[i];
    }

    if (m.nemobits) {
      if (cmodel)
        nemobits = (int *) mallocate((char *) nemobits, sizeof(int), true);
      else {
        if (nemobits) delete nemobits;
        nemobits = new int;
      }
      *nemobits = *m.nemobits;
    }

    if (m.timu) {
      if (cmodel)
        timu = (float *) mallocate((char *) timu, sizeof(float), true);
      else {
        if (timu) delete timu;
        timu = new float;
      }
      *timu = *m.timu;
    }
  }
  return *this;
}
// ============================================================================
// Destructor                                                                  
// we use the "free" C statement because io_nemo use malloc to alloc data and  
// you cannot "delete" a pointer allocated with malloc                         
ParticlesData::~ParticlesData()
{
  if (pos) {
    if (cmodel) { free ((float *) pos); }
    else        { delete [] pos;        }
  }
  if (vel) {
    if (cmodel) { free ((float *) vel); }
    else        { delete [] vel;        }
  }
  if (vel_norm) {
    if (cmodel) { free ((float *) vel_norm); }
    else        { delete [] vel_norm;        }
  }

  if (rho) {
    if (cmodel) { free ((float *) rho); }
    else        { delete [] rho;        }
  }
  if (rneib) {
    if (cmodel) { free ((float *) rneib); }
    else        { delete [] rneib;        }
  }
  if (temp) {
    if (cmodel) { free ((float *) temp); }
    else        { delete [] temp;        }
  }
  if (timu) {
    if (cmodel) { free ((float *) timu); }
    else        { delete timu;           }
  }
  if (nbody) {
    if (cmodel) { free ((int   *) nbody); }
    else        { delete nbody;           }
  }
  if (nemobits) {
    if (cmodel) { free ((int   *) nemobits); }
    else        { delete nemobits;           }
  }
  pos      = NULL;
  vel      = NULL;
  vel_norm = NULL;
  rho      = NULL;
  temp     = NULL;
  rneib    = NULL;
  timu     = NULL;
  nbody    = NULL;
  nemobits = NULL;
}

// ============================================================================
// ParticlesData::mallocate                                                    
char * ParticlesData::mallocate(char * p, int lg, bool force)
{
  char * ptr;

  if (p == NULL || force) {
    if (force && p ) {
      free ((char *) p);
    }
    ptr = (char *) malloc(sizeof(char)*lg);
    if (!ptr) {
      std::cerr << "[allocate_pointer], allocation memory error, aborted\n";
      exit(1);
    }
    return ptr;
  }
  else /* assume p is already  allocated */
    return p;
}
// ============================================================================
// ParticlesData::allocVar                                                     
int ParticlesData::allocVar()
{
  if (cmodel) {
    nbody    = (int   *) mallocate((char *) nbody   , sizeof(int));
    nemobits = (int   *) mallocate((char *) nemobits, sizeof(int));
    timu     = (float *) mallocate((char *) timu    , sizeof(float));
  } else {
    nbody    = new int;
    nemobits = new int;
    timu     = new float;
  }
  *nbody = *nemobits = -1;
  *timu  = -1.;
  return 1;
}
// ============================================================================
// ParticlesData::computeVelNorm()                                             
void ParticlesData::computeVelNorm()
{
  if (vel) {
    if (! vel_norm) {
      if (cmodel)
        vel_norm = (float *) mallocate((char *) vel_norm, sizeof(float)*(*nbody),true);
      else
        vel_norm = new float[*nbody];
    }
    bool first=true;
    for (int i=0; i<(*nbody); i++) {
      float vx=vel[i*3];
      float vy=vel[i*3 + 1];
      float vz=vel[i*3 + 2];
      vel_norm[i] = sqrt(vx*vx + vy*vy + vz*vz);
      if (first) {
        first=false;
        max_vel_norm = vel_norm[0];
      }
      else 
        if (vel_norm[i] > max_vel_norm) max_vel_norm = vel_norm[i];
    }
  }
}
// ============================================================================
// ParticlesData::computeMinMaxRho()                                           
int ParticlesData::computeMinMaxRho()
{
  if (rho) {
    max_rho = -1E9;
    min_rho = 1E9;
    // compute Min/Max
    for (int i=0; i<(*nbody); i++) {
      if (rho[i] != -1 ) {
	max_rho=std::max(max_rho,rho[i]);
	min_rho=std::min(min_rho,rho[i]);
      }
    }
    if (max_rho == -1E9 && min_rho == 1E9) {
      valid_rho = false;
    } 
    else {
      GlobalOptions::rho_exist = true;
      valid_rho = true;
    }
    std::cerr << "min rho="<<min_rho<<" max rho="<<max_rho<<"\n";
    // compute density histogram                                      
    // in a array of 100 bins, going from log(min rho) to log(max rho)
    // we compute the index in that range for the particles's density,
    // and we increment the index for that bin                        
    for (int i=0; i<100; i++) {
      density_histo[i] = 0;
    }
    if (valid_rho) {
      for (int i=0; i<(*nbody); i++) {
        if (rho[i] != -1 && rho[i] != 0) {
          int index=(log(rho[i])-log(min_rho))*99./(log(max_rho)-log(min_rho));
          //int index=((rho[i])-(min_rho))*99./((max_rho)-(min_rho));
          assert(index<100);
          density_histo[index]++;
        }
      }
    }
  }
  return 1;
}
// ============================================================================
// ParticlesData::computeMinMaxTemp()                                           
int ParticlesData::computeMinMaxTemp()
{
  if (temp) {
    max_temp = -1E9;
    min_temp = 1E9;
    // compute Min/Max
    for (int i=0; i<(*nbody); i++) {
      if (temp[i] != -1 ) {
        max_temp=std::max(max_temp,temp[i]);
        min_temp=std::min(min_temp,temp[i]);
      }
    }
    if (max_temp == min_temp) {
      valid_temp = false;
    } 
    else {
      valid_temp = true;
    }
    std::cerr << "min temp="<<min_temp<<" max temp="<<max_temp<<"\n";
    // compute density histogram                                      
    // in a array of 100 bins, going from log(min temp) to log(max temp)
    // we compute the index in that range for the particles's density,
    // and we increment the index for that bin                        
    for (int i=0; i<100; i++) {
      temp_histo[i] = 0;
    }
    if (valid_temp) {
      for (int i=0; i<(*nbody); i++) {
        if (temp[i] != -1 && temp[i] != 0) {
          int index=(log(temp[i])-log(min_temp))*99./(log(max_temp)-log(min_temp));
          //int index=((temp[i])-(min_temp))*99./((max_temp)-(min_temp));
          assert(index<100);
          temp_histo[index]++;
        }
      }
    }
  }
  return 1;
}
// ============================================================================

}
