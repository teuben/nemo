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
//
//   Stars and Gas
//
//     class - classification. For example,
//                         class < -100,000   Stellar halo particle
//             -100,000 <= class < 0          Stellar disk particle
//                         class = 0          Central BH
//                     0 < class              SPH particle
//
//     state - current state of particle
//             1 - star formation feedback
//             2 - star formation wind
//
//     mass  - mass
//     pos   - position vector
//     vel   - velocity vector
//     eps   - gravitational softening length
//     phi   - gravitational potential (per unit mass)
//     f     - force vector (per unit mass, i.e. acceleration)
//     rho   - local mass density
//     sfmt  - star formation epoch
//     sfmz  - metalicity
//
//   Gas
//
//     u      - internal energy (per unit mass)
//     udot1  - du/dt adiabatic
//     udot2  - du/dt viscous
//     udot3  - du/dt non-adiabatic
//     udot4  - du/dt star formation
//     d      - gas density
//     h      - gas smoothing length
//     divv   - divergence of velocity field
//     wgtmol - mean molecular weight
//     frcneu - fraction of neutral mass
// ============================================================================
#ifndef FTMIO_H
#define FTMIO_H

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include "componentrange.h"
#include "particlesdata.h"
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
namespace ftm {

class FtmComponent {
  public:
  FtmComponent(int *, int *, int, int , int);
  int ihole, igas, idisk;
  //CompRange halo,disk, hole, gas;
  glnemo::ComponentRange halo,disk, hole, gas;
};

class FtmIO{
public:
    FtmIO(const std::string);

    ~FtmIO();

    int open();
    int close();
    int read();
    int read(glnemo::ParticlesData * pdata, const int *index, const int nsel,
             const bool); 
    float * getMass()   const { return mass; }
    float * getPos()    const { return pos; }
    float * getVel()    const { return vel; }
    float   getTime()   const { return float(tframe);}
    int     getNtotal() const { return n0;}
    const FtmComponent * getComp() const { return comp;}
    
private:
  std::string filename;
  std::ifstream in;
  //header
  int dmpindx, ndim, eqnindx, lsfm,
      n0, n1, n2;             // no (total),n1 (DM + stars),n2 (gas)
  double version, gamma, poly, tframe;
  char dmp_date[30];
  //data
  float * mass, * pos, * vel;
  int * classe, * index1;
  float frame;
  //IO
  int readFRecord();
  int readHeader();
  int smartRead(char *, unsigned int);
  int skipBytes(unsigned int);
  void swapArrayIndex3D(float *, float *, int);
  void swapArrayIndex3D(float *, float *, int, int);
  void swapArrayIndex3D(float *, float *, int, int, const int*, const int);
  // skip Block
  inline void skipBlock() {
    int len1 = readFRecord();
    in.seekg(len1,std::ios::cur);
    int len2 = readFRecord();
    assert(in.good() && len1==len2);
  }
  // organize
  void classComponents();
  int ihole, igas, idisk;
  FtmComponent * comp;
  //control
  bool is_open, is_read;
  bool status;
  //fortran offset record length
  static int frecord_offset;
};
}
#endif
