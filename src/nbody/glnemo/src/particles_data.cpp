// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "particles_data.h"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <math.h>
// ============================================================================
// Constructor                                                                 
ParticlesData::ParticlesData()
{
  pos      = NULL;
  vel      = NULL;
  vel_norm = NULL;
  timu     = NULL;
  nbody    = NULL;
  nemobits = NULL;
  tree_depth = NULL;	
  coo_max[0] = coo_max[1] = coo_max[2] = 0.0;
  i_max[0]   = i_max[1]   = i_max[2]   = 0;
  tree_size_max = 1000000;
}
// ============================================================================
// copy Constructor                                                            
const ParticlesData& ParticlesData::operator=(const ParticlesData& m)
{
  if (m.nbody) {
    nbody = (int *) mallocate((char *) nbody, sizeof(int), true);
    *nbody = *m.nbody;

   tree_depth = (int *) mallocate((char *) tree_depth, sizeof(int)* (*nbody), true);

   if (m.pos) {
      pos = (float *) mallocate((char *) pos, sizeof (int) * 3 * (*nbody), true);
      //memcpy((float *) pos, (float *) m.pos, sizeof(float)* (*nbody) * 3);
      for (int i=0; i < (*nbody); i++) {
        pos[i*3]    = m.pos[i*3];
        pos[i*3+1]  = m.pos[i*3+1];
        pos[i*3+2]  = m.pos[i*3+2];
      }
    }

    if (m.vel) {
      vel = (float *) mallocate((char *) vel, sizeof (float) * 3 * (*nbody), true);
      //memcpy((float *) vel, (float *) m.vel, sizeof(float)* (*nbody) * 3);
      for (int i=0; i < (*nbody); i++) {
        vel[i*3]    = m.vel[i*3];
        vel[i*3+1]  = m.vel[i*3+1];
        vel[i*3+2]  = m.vel[i*3+2];
      }
    } 
    else {
        if (vel) 
          free ((float *) vel);     
        vel = NULL;
    }

    if (m.vel_norm) {
      vel_norm = (float *) mallocate((char *) vel_norm, sizeof (float) * (*nbody), true);
      //memcpy((float *) vel_norm, (float *) m.vel_norm, sizeof(float)* (*nbody));
      for (int i=0; i < (*nbody); i++) {
        vel_norm[i] = m.vel_norm[i];
      }
    } 
    else {
      if (vel_norm) 
	free ((float *) vel_norm);
      vel_norm = NULL;
    }
 
    //memcpy((float *) coo_max, (float *) m.coo_max, sizeof(float)* 3);
    //memcpy((int *) i_max, (int *) m.i_max, sizeof(int)* 3);
    for (int i=0; i <3; i++) {
      coo_max[i] = m.coo_max[i];
      i_max[i]   = m.i_max[i];
    }

    if (m.nemobits) {
      nemobits = (int *) mallocate((char *) nemobits, sizeof(int), true);
      *nemobits = *m.nemobits;
    }

    if (m.timu) {
      timu = (float *) mallocate((char *) timu, sizeof(float), true);
      *timu = *m.timu;
    }
}
 return *this;
}
// ============================================================================
// Destructor                                                                  
// we use the "free" C statement because io_nemo use malloc to alloc data and  
// you cannot "delete" a pointer allocate with malloc                          
ParticlesData::~ParticlesData()
{
  if (pos) 
 	free ((float *) pos);
  if (vel) 
	free ((float *) vel);
  if (vel_norm) 
	free ((float *) vel_norm);
  if (timu) 
	free ((float *) timu);
  if (nbody) 
	free ((int   *) nbody);
  if (nemobits) 
	free ((int   *) nemobits);
  if (tree_depth) 
	free ((int   *) tree_depth);

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
     std::exit(1);
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
  nbody    = (int   *) mallocate((char *) nbody   , sizeof(int));
  //tree_depth = (int   *) mallocate((char *) tree_depth, sizeof(int)* (*nbody));
  //tree_depth = (int   *) mallocate((char *) tree_depth, sizeof(int));
  nemobits = (int   *) mallocate((char *) nemobits, sizeof(int));
  timu     = (float *) mallocate((char *) timu    , sizeof(float));
  *nbody = *nemobits = -1;
  *timu  = -1.;
  return 1;
}
// ============================================================================
// ParticlesData::allocTree                                                     
int ParticlesData::allocTree()
{
  tree_depth = (int   *) mallocate((char *) tree_depth, sizeof(int)* (*nbody));
  for (int i=0; i<*nbody; i++) {
    tree_depth[i] = 1;
  }
  return 1;
}
// ============================================================================
// ParticlesData::computeVelNorm()                                             
void ParticlesData::computeVelNorm()
{
  if (vel) {
    if (! vel_norm) {
      vel_norm = (float *) mallocate((char *) vel_norm, sizeof(float)*(*nbody),true);
    }
    for (int i=0; i<(*nbody); i++) {
      float vx=vel[i*3];
      float vy=vel[i*3 + 1];
      float vz=vel[i*3 + 2];
      vel_norm[i] = sqrt(vx*vx + vy*vy + vz*vz);
    }
  }
}
// ============================================================================
