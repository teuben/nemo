// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014
//           Centre de donneeS Astrophysiques de Marseille (CeSAM)              
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Aix Marseille Universite, CNRS, LAM 
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================

/* 
	@author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
 */

#ifndef CPART_H
#define CPART_H

#include <string>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "cfortio.h"
#include "snapshotramses.h"

namespace uns {
class CParticles;
}

namespace ramses {

class CPart {
public:
  CPart(const std::string,const bool _v=true);
  ~CPart();
  void setBoundary(float x[6]) {
    xmin=x[0];
    xmax=x[1];
    ymin=x[2];
    ymax=x[3];
    zmin=x[4];
    zmax=x[5];
    
  }
  bool isValid();
  int loadData(uns::CParticles * particles,
               const unsigned int req_bits, const unsigned int comp_bits);
  int getNbody(int * dm, int * stars)    { 
    *dm    = ndm_box;
    *stars = nstar_box; 
    return nselect;
  }
    
private:
  bool verbose,valid;
  std::string infile,indir;
  int nselect;
  int npart,nstar,ncpu,ndim,nbody,ndm;
  int ndm_box, nstar_box;
  std::string s_run_index;

  float xmin,xmax,ymin,ymax,zmin,zmax;

  CFortIO part;
  int readHeader();
};
} // namespace ramses
#endif
