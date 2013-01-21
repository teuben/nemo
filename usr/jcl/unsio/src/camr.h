// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2013                                       
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
#include <cstdlib>
#include "cfortio.h"
#include "snapshotramses.h"

namespace uns {
class CParticles;
}
namespace ramses {

typedef struct  {
  double time;
  double boxlen, omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini;
  double aexp,hexp,aexp_old,epot_tot_int,epot_tot_old;
} Header;

class CAmr  {
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
      std::cerr << "min = "<< (int) x[6] << " lmax="<<lmax<<" lmin="<<lmin<<"\n";
      //exit(1);
    }
    bool isValid();
    int loadData(uns::CParticles * particles,
                 const unsigned int req_bits);
    int getNbody()    { return nbody;}
    Header * getHeader() { return &header; }


private:
    // some variables
    
    bool verbose,valid;
    std::string infile,testhydrofile, indir;
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
    Header header;
    // hydro
    int nvarh;
    double scale_nH;

};
} // end of namespace
#endif // CAMR_H
