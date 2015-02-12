// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                       
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================

/* 
	@author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
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
#include <cmath>
#include <map>
#include <climits>
#include <cstdlib>

namespace ramses {

typedef struct  {
  double time;
  double boxlen, omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini;
  double aexp,hexp,aexp_old,epot_tot_int,epot_tot_old;
} Header;

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
    Header * getHeader() { return &header; }
    double getMapInfo(std::string _s) {
      return mapinfo[_s];
    }
    int getNdim() { return ndim; }
signals:
    void stringStatus(const QString);
    void intStatus(const int);    
    
private:
    // some variables
    // return randomizely a negative or positive sign
    inline float getSign() {
      if (qrand() < RAND_MAX/2) return -1.;
      else return +1;
    }
    // return a random number between [0:val]
    inline double getRandom(const double val) {
      return val*(double) qrand()/(double) RAND_MAX;
    }
    // return a random number between [0:val]
    inline double getRandomGaussian(const double val) {

      double distance=val*(double) qrand()/(double) RAND_MAX;
      double g=2.;
      double pi=atan(1.0)*4.;
      double ig=1./g;
      double isqrtpi=1./sqrt(2.*pi);
      double r=(exp(-(distance)*(distance)*ig*ig*0.5)*isqrtpi*ig);
      return (r);
    }
    bool verbose,valid;
    std::string infile,testhydrofile,indir;
    int select,nselect;
    int nbody;
    std::string s_run_index,ordering;
  
    float xmin,xmax,ymin,ymax,zmin,zmax;
    int lmin,lmax;
    CFortIO  amr, hydro;
    // header
    int readHeader();
    Header header;
    // info file .txt
    std::string info_file;
    std::map<std::string, double > mapinfo;
    bool readInfoFile();
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
