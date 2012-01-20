// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================

/* 
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef CAMR_H
#define CAMR_H

#include <string>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "cfortio.h"
#include <QObject>

namespace ramses {
class CAmr : public QObject {
  Q_OBJECT
public:
    CAmr(const std::string,const bool _v=true);
    
    ~CAmr();
    void setBoundary(float x[8]) {
      xmin=x[0];
      xmax=x[1];
      ymin=x[2];
      ymax=x[3];
      zmin=x[4];
      zmax=x[5];
      
      if (x[7]==0.) {
        lmax= nlevelmax;
      } else {
        lmax = (int) x[7];
      } 
      lmin = std::min((int) x[6],lmax-1);                   
    }
    bool isValid();
    int loadData(float * pos=NULL, float * vel=NULL, float * rho=NULL,
                 float * rneib=NULL, float * temp=NULL,const int *index=NULL,
                 const int nsel=0,   const bool load_vel=false);
    int getNbody()    { return nbody;}
    
signals:
    void stringStatus(const QString);
    void intStatus(const int);    
    
private:
    // some variables
    
    bool verbose,valid;
    std::string infile,indir;
    int select,nselect;
    int nbody;
    std::string s_run_index,ordering;
  
    float xmin,xmax,ymin,ymax,zmin,zmax;
    int lmin,lmax;
    CFortIO  amr, hydro;
    int readHeader();
    // amr header variables
    static const double XH, mH, kB;
    int ncpu, ndim, nx, ny ,nz, nlevelmax, ngridmax, nboundary, ngrid_current;
    int twotondim;
    double xbound[3];
    // hydro
    int nvarh;
    double scale_nH;
};
} // end of namespace
#endif // CAMR_H
