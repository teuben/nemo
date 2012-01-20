// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "particlesdata.h"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <limits>
#include <iomanip>
#include "globaloptions.h"
#define MAX(A,B) ((A)>(B)?(A):(B))
#define PRINT_D if (1)
namespace glnemo {

// ============================================================================
// Constructor                                                                 
ParticlesData::ParticlesData(const ALLOC _model)
{
  pos      = NULL;
  vel      = NULL;
  vel_norm = NULL;

  rneib    = NULL;
  rho      = NULL;
  temp     = NULL;
  pressure = NULL;
  timu     = NULL;
  nbody    = NULL;
  nemobits = NULL;
  id.clear();
  cmodel   = _model;
  allocVar();
  coo_max[0] = coo_max[1] = coo_max[2] = 0.0;
  max_size=0.;
  i_max[0]   = i_max[1]   = i_max[2]   = 0;
  ipvs = 1; // density by default
}
// ============================================================================
// copy Constructor                                                            
const ParticlesData& ParticlesData::operator=(const ParticlesData& m)
{
  if (m.nbody) {
    cmodel = m.cmodel;
    ipvs   = m.ipvs;
    if (cmodel)
      nbody = (int *) mallocate((char *) nbody, sizeof(int), true);
    else {
      if (nbody) delete nbody;
      nbody = new int;
    }
    *nbody = *m.nbody;
    // Ids
    id = m.id;
    
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
    // Density "rho"
    if (m.rho) {
      if (rho) delete rho;	
      rho = new PhysicalData(PhysicalData::rho,*nbody);
      *rho = *m.rho;
    }
    else {
      if (rho) {
	if (cmodel) { free ((float *) rho); }
        else        { delete rho;        }
      }
      rho = NULL;
    }
    // neibourgh radius "rneib"
    if (m.rneib) {
      if (rneib) delete rneib;
      rneib = new PhysicalData(PhysicalData::neib,*nbody);
      *rneib = *m.rneib;
    }
    else {
      if (rneib) {
	if (cmodel) { free ((float *) rneib); }
        else        { delete rneib;        }
      }
      rneib = NULL;
    }
    // Temperature "temp"
    if (m.temp) {
      if (temp) delete temp;
      temp = new PhysicalData(PhysicalData::temperature,*nbody);
      *temp = *m.temp;
    }
    else {
      if (temp) {
        if (cmodel) { free ((float *) temp); }
        else        { delete temp;           }
      }
      temp = NULL;
    }
    // Pressure "pressure"
    if (m.pressure) {
      if (pressure) delete pressure;
      pressure = new PhysicalData(PhysicalData::pressure,*nbody);
      *pressure = *m.pressure;
    }
    else {
      if (pressure) {
        if (cmodel) { free ((float *) pressure); }
        else        { delete  pressure;        }
      }
      pressure = NULL;
    }

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
    else        { delete rho;        }
  }
  if (rneib) {
    if (cmodel) { free ((float *) rneib); }
    else        { delete rneib;        }
  }
  if (temp) {
    if (cmodel) { free ((float *) temp); }
    else        { delete temp;        }
  }
  if (pressure) {
    if (cmodel) { free ((float *) pressure); }
    else        { delete pressure;        }
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
  pressure = NULL;
  rneib    = NULL;
  timu     = NULL;
  nbody    = NULL;
  nemobits = NULL;
  id.clear();
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
// ParticlesData::getPhysData()                                                
// Return data at the index of the Physical value selected ipvs                
// index = 0 -> rho, 1 -> temp, 2 -> pressure
PhysicalData * ParticlesData::getPhysData(int index) const
{
  if (index==-1) {
    index=ipvs;
  }
  PhysicalData * ptr;
  switch (index) {
  case 0:
    //assert(rneib!=NULL);
    ptr = rneib;
    break;
  case 1:
    //assert(rho!=NULL);
    ptr = rho;
    break;
  case 2:
    //assert(temp!=NULL);
    ptr = temp;
    break;
  case 3:
    //assert(pressure!=NULL);
    ptr = pressure;
    break;
  case 4: // temperature sorted by density
    //assert(temp!=NULL);
    ptr = temp;
    break;
    
    default:
    assert(0);
  }
  return ptr;
}
// ============================================================================
// ParticlesData::computeMaxSize()                                           
void ParticlesData::computeMaxSize()
{
  if (pos && max_size==0.) {
    coo_max[0] = pos[0];
    coo_max[1] = pos[1];
    coo_max[2] = pos[2];
    coo_min[0] = pos[0];
    coo_min[1] = pos[1];
    coo_min[2] = pos[2];
    for (int i=1; i<(*nbody); i++) {
      coo_max[0] = std::max(coo_max[0],pos[i*3+0]); 
      coo_max[1] = std::max(coo_max[1],pos[i*3+1]); 
      coo_max[2] = std::max(coo_max[2],pos[i*3+2]);
      coo_min[0] = std::min(coo_min[0],pos[i*3+0]); 
      coo_min[1] = std::min(coo_min[1],pos[i*3+1]); 
      coo_min[2] = std::min(coo_min[2],pos[i*3+2]);
    }
    //float max=std::max(std::max(coo_max[0],coo_max[1]),coo_max[2]);
    //float min=std::min(std::min(coo_min[0],coo_min[1]),coo_min[2]);
    max_size=sqrt(pow(coo_max[0]-coo_min[0],2)+pow(coo_max[1]-coo_min[1],2)+pow(coo_max[2]-coo_min[2],2));
    //max_size = max-min;
    if (1) {
      PRINT_D std::cerr << "Max coordinates \n";
      PRINT_D std::cerr << coo_max[0] << " " << coo_max[1] << " " << coo_max[2] << "\n";
      PRINT_D std::cerr << coo_min[0] << " " << coo_min[1] << " " << coo_min[2] << "\n";
      PRINT_D std::cerr << "max_size = " << max_size << "\n";
    }
  }
}
// ============================================================================
// PhysicalData Class implementation                                           
// ============================================================================

// ============================================================================
// Constructor                                                                 
PhysicalData::PhysicalData(const PHYS _type,const int _nbody,const ALLOC _model)
{
  data = NULL;
  type = _type;
  nbody = _nbody;
  cmodel = _model;
  if (nbody>0) {
    if (cmodel)
      data = (float *) ParticlesData::mallocate((char *) data, sizeof (float) * (nbody), true);
    else {
      data = new float[nbody];
    }
  }
  valid = false;
}
// ============================================================================
// copy Constructor                                                            
const PhysicalData& PhysicalData::operator=(const PhysicalData& m)
{
  nbody = m.nbody;
  valid = m.valid;
  cmodel= m.cmodel;
  type  = m.type;
  if (nbody) {
    max = m.max;
    min = m.min;
    if (cmodel) {
      if (data) free ((float *) data);
      data = (float *) ParticlesData::mallocate((char *) data, sizeof (float) * (nbody), true);
    }
    else {
      if (data) delete [] data;
      data = new float[nbody];
    }
    for (int i=0; i < (nbody); i++) {
      data[i] = m.data[i];
    }
    for (int i=0; i < 100; i++) {
      data_histo[i] =m.data_histo[i];
    }
  }
  return *this;
}
// ============================================================================
// Destructor                                                                  
PhysicalData::~PhysicalData()
{
  if (data) {
    if (cmodel) { free ((float *) data); }
    else        { delete [] data;        }
  }
}
// ============================================================================
// PhysicalData::computeMinMax()                                           
int PhysicalData::computeMinMax()
{
  if (data) {
    max = -1E9;
    min = 1E9;
    // compute Min/Max
    for (int i=0; i<(nbody); i++) {
      if (data[i] != -1 ) {
	max=std::max(max,(double) data[i]);
	min=std::min(min,(double) data[i]);
      }
    }       
    std::cerr << "-------------------------------\n";
    valid = true;
    if ((max == -1E9 && min == 1E9 )) {
      valid = false;
    }
    if (min>0 && max>0) {
      if ((max == -1E9 && min == 1E9 )|| (max == min) ||
          max >  std::numeric_limits<double>::max() ||
          max <  std::numeric_limits<double>::min() ||
          min >  std::numeric_limits<double>::max() ||
          min <  std::numeric_limits<double>::min() 
        ) {
        valid = false;
      } 
    }
    switch (type) {
    case PhysicalData::neib :
      std::cerr << "Hsml        range :\n";
      break;
    case PhysicalData::rho : 
      if (valid) GlobalOptions::rho_exist         = true;
      else       GlobalOptions::rho_exist         = false;
      std::cerr << "Density     range :\n";
      break;
    case PhysicalData::temperature : 
      if (valid) GlobalOptions::temperature_exist = true;
      else       GlobalOptions::temperature_exist = false;
      std::cerr << "Temperature range :\n";
      break;  
    case PhysicalData::temperaturesd : 
      if (valid) GlobalOptions::temperature_exist = true;
      else       GlobalOptions::temperature_exist = false;
      std::cerr << "Temperature range :\n";
      break;     
    case PhysicalData::pressure :
      if (valid) GlobalOptions::pressure_exist    = true;
      else       GlobalOptions::pressure_exist    = false;
      std::cerr << "Pressure    range :\n";
      break;                        
    } 
    std::cerr << "VALID =["<<valid<<"]\n";
    std::cerr <<" min = "<< std::scientific << std::setw(10) << min
              <<"\n max = "<< std::setw(10)<<max<<"\n";
    if (max >  std::numeric_limits<double>::max()) {
      std::cerr << "max is inf....\n"; 
    }
    // compute density histogram                                      
    // in a array of 100 bins, going from log(min data) to log(max data)
    // we compute the index in that range for the particles's density,
    // and we increment the index for that bin                        
    for (int i=0; i<100; i++) {
      data_histo[i] = 0;
    }
    if (min <=0) {
      std::cerr << "Min is negative, rescaling data...\n";
      double offset;
      if (min==0) {
        min=1E-9;
        offset=1E-9;
      }
      else offset=min*-1.0001;//(1+1.E-12);
      double mmin=0.0;
      bool first=true;
      for (int i=0; i<(nbody); i++) {
        if (data[i] != -1 ) {
          float dd=data[i];
          data[i] += offset;
          if (first) {
            mmin = data[i];
            if (data[i] != 0.0)
              first=false;
          } else {
            if (data[i]!=0.0)
              mmin = std::min(mmin,(double)data[i]);
            else {
              std::cerr << "i="<<i<<" data="<<dd<<"\n";
            }
          }
        }
      }
      min=1E-12;
      min=mmin;
      max=max+offset;
      std::cerr <<" new min = "<< std::scientific << std::setw(10) << min
                <<"\n new max = "<< std::setw(10)<<max<<"\n";
    }
    if (valid) {
      for (int i=0; i<(nbody); i++) {
        if (data[i] != -1 && data[i] != 0) {
          int index=(log(data[i])-log(min))*99./(log(max)-log(min));
          //int index=((data[i])-(min))*99./((max)-(min));
          assert(index<100);
          data_histo[index]++; // data_histo stores the number of particles foreach percentage
                               // of physical value from 0 to 99 %
        }
      }
    }
  }
  return 1;
}
} // namespace glnemo
