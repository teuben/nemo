// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008 / 2009                                    
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// ----------------------------------------------------------------------------
// gadget2nemo.cc                                                              
// 15-Oct-08 : V 3.1 new version, gadget2 support, happy g++ 4.3.2  (JCL)    
// 07-Oct-09 : V 3.2 new gadget2 reader
// ----------------------------------------------------------------------------
#include <iostream>                                   // C++ I/O
#include <fstream>                                    // C++ file I/O
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <snapshot/snapshot.h>
#include <assert.h>

#include "gadgetio.h"
#include "userselection.h"

extern "C" {
#include <nemo.h>                                     // NEMO basics
  int io_nemo(char *, char *,...);
}

using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
// NEMO parameters                                                              
const char * defv[] = {
  "in=???\n             GADGET (1 or 2 little/big endian) input       ",
  "out=???\n            NEMO output                                   ",
  "select=disk\n        requested data:                             \n"
  "                   disk, gas, halo, bulge, stars, bndry, all     \n"
  "                   you can request several data separated with   \n"
  "                   \",\" like select=disk,stars,gas                ",
  "comp=mxv\n           component requested to be saved             \n"
  "                   m->mass, x->positions, v->velocities          ",
  "verb=f\n             verbose mode                                  ",
  "VERSION=3.2\n         compiled on <"__DATE__"> JCL                 ",
  NULL
};
const char * usage="Convert GADGET (1 or 2 little/big endian) snapshot to NEMO snapshot";

int    *iobits=NULL;
float * mass=NULL, * pos=NULL, * vel=NULL;
::string comp;

//------------------------------------------------------------------------------
// outNemo                                                                      
int outNemo(gadget::GadgetIO * gadget_io, 
            glnemo::UserSelection * user_select,
            ::string out)
{
  iobits = new int;
  *iobits = TimeBit;
  int * N = new int;
  // allocate memory
  if (strchr(comp,'x')) {
    (*iobits) |= PosBit;
    pos  = gadget_io->getPos();
  }
  if (strchr(comp,'v')) {
    (*iobits) |= VelBit;
    vel  = gadget_io->getVel();
  }
  
  if (strchr(comp,'m')) {
    (*iobits) |= MassBit;
    mass =  gadget_io->getMass();
  }

  // save nemo spnashot
  *N = user_select->getNSel();
  const float * tps = gadget_io->getTime();
  if ( ! io_nemo(out,(char *)"save,float,n,t,m,x,v,b,info",
		 &N,&tps,&mass,&pos,&vel,&iobits)) {
    cerr << "Unable to save snapshot [" << out << "]\n";
    cerr << "Aborted\n";
  } else {
  }
  
  return 1;
}

//------------------------------------------------------------------------------
// main                                                                         
int main(int argc, char ** argv )
{
  ::string in,out;
  std::string interface_type="Gadget";
  bool verb;
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));

  // Get parametersg
  in          = getparam((char *) "in");
  out         = getparam((char *) "out");
  comp        = getparam((char *) "comp");
  std::string select=getparam((char *) "select");
  verb        = getbparam((char *) "verb");
  

  // Read GADGET snapshot
  gadget::GadgetIO * gadget_io = new gadget::GadgetIO(in,verb);
  int fail = gadget_io->open(in);
  if (!fail) {
    std::ostringstream stm;
    stm << gadget_io->getVersion();
    interface_type += " " + stm.str();
    glnemo::ComponentRangeVector crv = gadget_io->getCRV();    // get component range       
    if (verb)  {                                                // verbose ?                 
      std::cerr << "Version : " << interface_type << "\n";
      glnemo::ComponentRange::list(&crv);                      // display list of components
    }
    glnemo::UserSelection * user_select = 
      new glnemo::UserSelection();                             // new user select obj       
    user_select->setSelection(select,&crv);                    // select according to crv   
    gadget_io->read(user_select->getIndexesTab(),user_select->getNSel()); // read Gadget    
    outNemo(gadget_io,user_select,out);                        // save to NEMO format
  } 
  else {
    std::cerr << "File["<<in<<"] is not a Gadget file, aborting...\n";
    exit(1);
  }

  //   finish NEMO
  finiparam();
}
// ----------- End Of [gadget2nemo.cc] ------------------------------------
