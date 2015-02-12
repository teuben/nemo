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

#ifndef TIPSYIO_H
#define TIPSYIO_H
#include <iostream>
#include <fstream>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <string>
#include <vector>
#include <assert.h>
#include "componentrange.h"

namespace tipsy {

#define MAXDIM 3

typedef float Real;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} ;

//static struct gas_particle *gas_particles;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} ;

//static struct dark_particle *dark_particles;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} ;

//static struct star_particle *star_particles;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

class TipsyIO {
public:
  TipsyIO(const std::string&);
  ~TipsyIO();

  int open(const std::string);
  int close();

  int loadData(std::vector <int> * id,float * pos, float * vel, float * rho, float * rneib, float * temp,const int *index,
               const int nsel,   const bool load_vel, const bool take_gas, const bool take_halo, const bool take_stars);
  float * getMass()   const { return mass; }
  float * getPos()    const { return pos; }
  float * getVel()    const { return vel; }
  float   getTime()   const { return float(tframe);}
  int     getNtotal() const { return npartTotal;}
  const glnemo::ComponentRangeVector getCRV() const { return crv;}

private:
  // variables
  XDR xdr;
  std::string filename,file0;
  FILE * in;
  int multiplefiles;
  bool lonely_file;
  //data
  float * mass, * pos, * vel, * intenerg,  * intenergp,* tempp, * rhop;
  float tframe;
  struct dump header;
  int npartTotal, npart;
  glnemo::ComponentRangeVector  crv;


  //control
  bool is_open, is_read;
  bool status;

  // methods
  int readHeader();
  int xdr_header(XDR *xdrs, struct dump *header);
  int xdr_gas(XDR *xdrs,struct gas_particle *p);
  int xdr_dark(XDR *xdrs,struct dark_particle *p);
  int xdr_star(XDR *xdrs,struct star_particle *p);
  int storeComponents();

};
}

#endif // TIPSYIO_H
