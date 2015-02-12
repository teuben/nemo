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


#include "tipsyio.h"
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <string.h>
#include <cmath>

namespace tipsy {

// ============================================================================
// Constructor
TipsyIO::TipsyIO(const std::string & _f)
{
  filename = _f;
  status=false;
  is_open=false;
  is_read=false;
  mass   = NULL;
  pos    = NULL;
  vel    = NULL;
  intenerg=NULL;
}

// ============================================================================
// Destuctor
TipsyIO::~TipsyIO()
{
  if (mass)     delete [] mass;
  if (pos)      delete [] pos;
  if (vel)      delete [] vel;
  if (intenerg) delete [] intenerg;
}
// ============================================================================
// open() :
// open file and return :
// 0 : success
// 1 : unable to open
// 2 : not a TIPSY file
int TipsyIO::open(const std::string myfile)
{
  int fail=0;
  std::cerr << "In TipsyIO open file : " << myfile << "\n";
  in = fopen(myfile.c_str(),"r");
  if ( ! in) {
   //fclose(in);
   //in.clear(); // mandatory under win32 !
   fail = 1;               // unable to open
   std::cerr << "In TipsyIO, failed to open...\n";
  }
  else {
    is_open=true;
    fail = readHeader();    // try to read header
    if (!fail) status=true; // valid header
    else close();           // not valid header
  }
  return fail;
}
// ============================================================================
// readHeader
int TipsyIO::readHeader()
{
  int status=0;
  xdrstdio_create(&xdr,in,XDR_DECODE);

  /* read header */
  xdr_header(&xdr,&header);

  if (header.ndim==3) {
    status=0; // it's a valid tipsy file
    tframe = header.time;
    storeComponents();
  } else {
    status=1; // not a tipsy file
  }
  return status;
}
// ============================================================================
// close() file
int TipsyIO::close()
{
  if (is_open) fclose(in);
  is_open = false;
  return 1;
}
// ============================================================================
//
int TipsyIO::loadData(std::vector <int> * id, float * pos, float * vel, float * rho, float * rneib,
                      float * temp,const int *index, const int nsel,
                      const bool load_vel, const bool take_gas, const bool take_halo, const bool take_stars)
{
  struct star_particle sp;
  struct dark_particle dp;
  struct gas_particle gp;

  int i=0,cpt=0; // particles counter

  // read gas
  while (cpt<header.nsph) {
    xdr_gas(&xdr,&gp);
    if (take_gas) {
      //mass[i]  =gp.mass;
      pos[i*3+0]=gp.pos[0];
      pos[i*3+1]=gp.pos[1];
      pos[i*3+2]=gp.pos[2];
      if (load_vel) {
        vel[i*3+0]=gp.vel[0];
        vel[i*3+1]=gp.vel[1];
        vel[i*3+2]=gp.vel[2];
      }
      rho[i]   =gp.rho;
      temp[i]=gp.temp;
      rneib[i]=gp.hsmooth;
      i++; // one more particle to store
    }
    cpt++; // one more particle to read
  }

  // read dark matter halo
  cpt=0;
  while (cpt<header.ndark) {
    xdr_dark(&xdr,&dp);
    if (take_halo) {
      //mass[i]  =gp.mass;
      pos[i*3+0]=dp.pos[0];
      pos[i*3+1]=dp.pos[1];
      pos[i*3+2]=dp.pos[2];
      if (load_vel) {
        vel[i*3+0]=dp.vel[0];
        vel[i*3+1]=dp.vel[1];
        vel[i*3+2]=dp.vel[2];
      }
      i++; // one more particle to store
    }
    cpt++; // one more particle to read
  }

  // read stars
  cpt=0;
  while (cpt<header.nstar) {
    xdr_star(&xdr,&sp);
    if (take_stars) {
      //mass[i]  =gp.mass;
      pos[i*3+0]=sp.pos[0];
      pos[i*3+1]=sp.pos[1];
      pos[i*3+2]=sp.pos[2];
      if (load_vel) {
        vel[i*3+0]=sp.vel[0];
        vel[i*3+1]=sp.vel[1];
        vel[i*3+2]=sp.vel[2];
      }
      i++; // one more particle to store
    }
    cpt++; // one more particle to read
  }
  assert(i==nsel);
  return 1;
}

// ============================================================================
//
int TipsyIO::storeComponents()
{
  glnemo::ComponentRange cr;
  // all
  npartTotal = header.nbodies;
  cr.setData(0,npartTotal-1);
  cr.setType("all");
  crv.clear();
  crv.push_back(cr);
  // components
  int start=0;
  //gas
  if (header.nsph) {
    cr.setData(start,start+header.nsph-1,"gas");
    crv.push_back(cr);
    start+=header.nsph;
  }
  //dm
  if (header.ndark) {
    cr.setData(start,start+header.ndark-1,"halo");
    crv.push_back(cr);
    start+=header.ndark;
  }
  //stars
  if (header.nstar) {
    cr.setData(start,start+header.nstar-1,"stars");
    crv.push_back(cr);
    start+=header.nstar;
  }
  return 1;
}

// ============================================================================
//
int TipsyIO::xdr_header(XDR *xdrs, struct dump *header)
{
  int pad=0;

  if (xdr_double(xdrs,&header->time) != TRUE) return 0;
  if (xdr_int(xdrs,&header->nbodies) != TRUE) return 0;
  if (xdr_int(xdrs,&header->ndim) != TRUE) return 0;
  if (xdr_int(xdrs,&header->nsph) != TRUE) return 0;
  if (xdr_int(xdrs,&header->ndark) != TRUE) return 0;
  if (xdr_int(xdrs,&header->nstar) != TRUE) return 0;
  if (xdr_int(xdrs,&pad) != TRUE) return 0;
  return 1;
}

// ============================================================================
//
int TipsyIO::xdr_gas(XDR *xdrs,struct gas_particle *p)
{
  if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->rho) != TRUE) return 0;
  if (xdr_float(xdrs,&p->temp) != TRUE) return 0;
  if (xdr_float(xdrs,&p->hsmooth) != TRUE) return 0;
  if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
  if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
  return 1;
}

// ============================================================================
//
int TipsyIO::xdr_dark(XDR *xdrs,struct dark_particle *p)
{
  if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
  if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
  return 1;
}

// ============================================================================
//
int TipsyIO::xdr_star(XDR *xdrs,struct star_particle *p)
{
  if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
  if (xdr_float(xdrs,&p->tform) != TRUE) return 0;
  if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
  if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
  return 1;
}
} // end of namespace tipsy

