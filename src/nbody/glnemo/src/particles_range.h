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
#ifndef PARITCLES_RANGE_H
#define PARITCLES_RANGE_H
#include <qcolor.h>
#include <qnamespace.h>
#include <iostream>
#include <vector>

using namespace std;

class ParticlesRange;

typedef vector <ParticlesRange> ParticlesRangeVector;

class ParticlesRange {
 public:

  ParticlesRange();

  ~ParticlesRange();

  // method
  char * parseString(const char * select_string, const int  
                     nbody,ParticlesRangeVector * prv );
  void printRange();

  // particles index variable
  int npart;            // #particles
  int first_part;       // index of the first particle
  int last_part;        // index of the last particle
  int step_part;        // incremental step between first and last
  QColor col;           // particles color
  bool   is_visible;    // TRUE if particles are visible

  static int nb_select; 

 private:
  // method
  int  parseSelectedString(char * select_sting, const int nbody,
                         ParticlesRangeVector * prv);


};
#endif
// ============================================================================
