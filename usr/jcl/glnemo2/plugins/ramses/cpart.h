// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/* 
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */

#ifndef CPART_H
#define CPART_H

#include <string>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "cfortio.h"
#include <QObject>

namespace ramses {
class CPart: public QObject {
  Q_OBJECT
public:
  CPart(const std::string,const int select,const bool _v=true);
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
  int loadData(bool , bool, float * pos=NULL, float * vel=NULL,const int *index=NULL,
               const int nsel=0,   const bool load_vel=false, const int namr_box=0);
  int getNbody(int * dm, int * stars)    { 
    *dm    = ndm_box;
    *stars = nstar_box; 
    return nselect;
  }
signals:
    void stringStatus(const QString);
    void intStatus(const int);    
    
private:
  bool verbose,valid;
  std::string infile,indir;
  int select,nselect;
  int npart,nstar,ncpu,ndim,nbody,ndm;
  int ndm_box, nstar_box;
  std::string s_run_index;

  float xmin,xmax,ymin,ymax,zmin,zmax;

  CFortIO part;
  int readHeader();
};
} // namespace ramses
#endif
