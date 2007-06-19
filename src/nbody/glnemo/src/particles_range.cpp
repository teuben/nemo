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
//                                                                             
// ParticlesRange class definition                                             
//                                                                             
// Parse selected string                                                       
// ============================================================================

#include <string.h>
// Nemo stuffs
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include <assert.h>

#include "particles_range.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"

//int ParticlesRange::nb_select=0;

// ============================================================================
// constructor                                                                 
ParticlesRange::ParticlesRange():VirtualParticlesSelect()
{
  v_type = 1; // record VirtualParticlesSetect type
  index_tab = NULL;
  ni_index  = 0;
}
// ============================================================================
// destructor                                                                  
ParticlesRange::~ParticlesRange()
{
  //nb_select--;
  if (index_tab) {
    delete [] index_tab;
  }
}

// ============================================================================
// ParticlesRange::defaultIndexTab()                                           
// fill the index tab array with the default number of particles               
inline int ParticlesRange::defaultIndexTab()
{
  if (index_tab) {
    delete [] index_tab;
  }
  index_tab = new int[npart];
  ni_index = 0;
  for (int i=0; i<npart; i+=step_part) {
    index_tab[ni_index++] = first_part+i;
  }
  return 1;
}
// ============================================================================
// ParticlesRange::getIndex()                                                  
// return index of the particle                                                
int ParticlesRange::getIndex(int index)
{
  return (first_part+index);
}

// ============================================================================
// ParticlesRange::parseSelectedString                                         
// Use nemoinpi engine to find out how many particles in the selected string,  
// THANKS to Peter Teuben, A LOT !!!                                           
int ParticlesRange::parseSelectedString(char * select_string, const int nbody, 
				      ParticlesSelectVector * psv)
{
  int * int_array = new int[nbody];
  
  PRINT_D std::cerr << "In parseSemicolon2...["<< select_string 
               <<  "and nbody = " << nbody << "]\n";
  if (  strcmp(select_string,"all")) {
    npart = nemoinpi(select_string, int_array, nbody);
    if (npart <=0 ) {
      std::cerr << "nemoinpi = [" << select_string << "] npart = "<<npart
           <<" with nbody=["<<nbody<<"]\n";
      exit(1);
    }
#if 0
    // Correct npart value if selected_range is out of nbody
    int nbody_out = 0;
    for (int i=0; i<npart; i++) {
      int p_index=int_array[i];
      if (p_index < nbody) {
        nbody_out++;
      }
    }
    PRINT_D std::cerr << "NBODY out = " << nbody_out 
                 << " and nbody =" << nbody <<"\n";
    npart = nbody_out;
#endif       
  }
  else {  // select all the particles
    npart=nbody;
  }
  delete [] int_array;  // useless anymore
  PRINT_D std::cerr << "In parseSemicolon2 npart =["<< npart << "]\n";
  // rescale particle range
  //if (nb_select == 1 ) { // first object
  if (psv->size() == 0 ) { // first object
    first_part = 0;
    last_part  = npart-1;
  } else {
    //first_part = (*psv)[nb_select-2].vps->last_part + 1;
    first_part = (*psv)[psv->size()-1].vps->last_part + 1;
    last_part  = first_part + npart - 1;
  }
  step_part = 1;   
  // allocate memory for index_tree

  defaultIndexTab();
  return npart;
}
// ============================================================================
