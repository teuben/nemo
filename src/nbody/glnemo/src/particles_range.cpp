// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
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

int ParticlesRange::nb_select=0;

// ============================================================================
//
ParticlesRange::ParticlesRange()
{
  //QColor modulo_col[5] = { Qt::white, Qt::green, Qt::yellow, Qt::red, Qt::blue };
  int modulo_col[5][3] = { 
                            { 255, 75, 39 },
                            { 214,214, 52 },
                            { 114,214, 32 },
                            {  58, 61,214 },
                            { 214, 47,197 }
                          };
  nb_select++;
  PRINT_D cerr << "ParticlesRange Constructor : " << nb_select << "\n";
  npart=-1;
  first_part=-1;
  last_part=-1;
  step_part=-1;
  //col = modulo_col[nb_select%5];
  col = QColor(modulo_col[(nb_select-1)%5][0],
               modulo_col[(nb_select-1)%5][1],
               modulo_col[(nb_select-1)%5][2]);
  is_visible=TRUE;
}
// ============================================================================
//
ParticlesRange::~ParticlesRange()
{
  //nb_select--;
}
// ============================================================================
//
void ParticlesRange::printRange()
{
  PRINT_D cerr << "Npart       = " << npart      << "\n";
  PRINT_D cerr << "First_part  = " << first_part << "\n";
  PRINT_D cerr << "Last_part   = " << last_part  << "\n";
  PRINT_D cerr << "Step_part   = " << step_part  << "\n";
  //  cerr << "Color       = " << QString(col) << "\n";
  PRINT_D cerr << "Is visible? = " << is_visible << "\n";
}
// ============================================================================
// parse 'select_string' according to the 'nemoinpi' rules.                    
// return the string at the position after the next 'coma' otherwise NULL
//
char * ParticlesRange::parseString(const char * select_string, const int nbody, 
				 ParticlesRangeVector * prv)
{
  char * status;
#if 0
  cerr << "before PRV size :" << prv->size() << "\n";
  prv->push_back(*this);
  cerr << "after PRV size  :" << prv->size() << "\n";
#endif
  PRINT_D cerr << "In parseString...["<< select_string << "]\n";

  char * c = strchr(select_string,',');
  int sup;
  if ( c) {
    status = c+1;
    sup = c-select_string;
  } else {
    status = NULL;
    sup = strlen(select_string)+1;
  }
  char tmp[100];
  strncpy(tmp,select_string,sup);
  tmp[sup] = '\0';
  parseSelectedString(tmp,nbody,prv);

  return status;
}

// ============================================================================
// Use nemoinpi engine to find out how many particles in the
// selected string, THANKS Peter, A LOT !!!
//
int ParticlesRange::parseSelectedString(char * select_string, const int nbody, 
				      ParticlesRangeVector * prv)
{
  int * int_array = new int[nbody];
  
  PRINT_D cerr << "In parseSemicolon2...["<< select_string <<
  "and nbody = " << nbody << "]\n";
  if (  strcmp(select_string,"all")) {
    npart = nemoinpi(select_string, int_array, nbody);
    if (npart <=0 ) {
      cerr << "nemoinpi = [" << select_string << "] npart = "<<npart
      <<" with nbody=["<<nbody<<"]\n";
      exit(1);
    }
    // Correct npart value if selected_range is out of nbody
    int nbody_out = 0;
    for (int i=0; i<npart; i++) {
      int p_index=int_array[i];
      if (p_index < nbody) {
        nbody_out++;
      }
    }
    PRINT_D cerr << "NBODY out = " << nbody_out << " and nbody =" <<
    nbody <<"\n";
    npart = nbody_out;
    delete int_array;  // useless anymore
  }
  else {  // select all the particles
    npart=nbody;
  }
  PRINT_D cerr << "In parseSemicolon2 npart =["<< npart << "]\n";
  // rescale particle range
  // 
  if (nb_select == 1 ) { // first object
    first_part = 0;
    last_part  = npart-1;
  } else {
    //cerr << "(*prv)["<<nb_select-2<<"].last_part=" << (*prv)[nb_select-2].last_part <<"\n";
    first_part = (*prv)[nb_select-2].last_part + 1;
    last_part  = first_part + npart - 1;
  }
  step_part = 1;   
  return npart;
}
// ============================================================================
