// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2013
// e-mail:   Jean-Charles.Lambert@oamp.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           Pole de l'Etoile, site de Chateau-Gombert
//           38, rue Frederic Joliot-Curie
//           13388 Marseille cedex 13 France
//           CNRS U.M.R 6110
// ============================================================================

/**
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef SNAPSHOTRAMSES_H
#define SNAPSHOTRAMSES_H

#include "snapshotinterface.h"
//#include "camr.h"
//#include "cpart.h"

namespace ramses {
class CAmr;
class CPart;
}
namespace uns {

class CParticles {
public:
  CParticles() {
    ntot=ngas=ndm,nstars=0;
    load_bits=0;
  }
  std::vector <float> pos,vel,mass,hsml,rho,temp,age,metal;
  std::vector <int> indexes,id;
  int ntot, ngas, ndm, nstars;
  unsigned int load_bits;
};

class CSnapshotRamsesIn: public CSnapshotInterfaceIn {

public:
  CSnapshotRamsesIn(const std::string, const std::string, const std::string, const bool verb=false);

  ~CSnapshotRamsesIn();
  // pure virtual function implemented
  ComponentRangeVector * getSnapshotRange();
  int nextFrame(uns::UserSelection &);
  bool getData(const std::string,int *n,float **);
  bool getData(const std::string,       float * );
  bool getData(const std::string,int *n,int   **);
  bool getData(const std::string,       int   * );
  bool getData(const std::string, const std::string ,int *,float **);
  bool getData(const std::string, const std::string ,int *,int   **);
 int close();

private:
 ramses::CAmr * amr;
 ramses::CPart * part;
 CParticles * particles;
 bool first_loc;

 int reorderParticles(uns::UserSelection & );
};
}
#endif // SNAPSHOTRAMSES_H

