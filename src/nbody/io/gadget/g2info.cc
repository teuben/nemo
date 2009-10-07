// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

// -----------------------------------------------------------------------------
// g2info.cc                                                                    
// 14-oct-08 : 1.01 : happy g++ 4.3.2                                           
// -----------------------------------------------------------------------------
#include <iostream>                                   // C++ I/O
#include <fstream>                                    // C++ file I/O
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <snapshot/snapshot.h>
#include "gadgetio.h"
#include <nemo.h>

using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
// NEMO parameters
const char * defv[] = { 
  "in=???\n             GADGET input                                 ",
  "select=disk\n        info requested                               ",
  "verb=f\n             verbose info                                 ",
  "VERSION=1.01\n       compiled on <"__DATE__"> JCL                ",
  NULL
};
const char * usage="Display Gadget2 file information";


//------------------------------------------------------------------------------
//                            GADGET Functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  char * in;
  bool verb;

  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));

  // Get parameters
  in                 = getparam((char *) ("in"));
  std::string select = getparam((char *) ("select"));
  verb               = getbparam((char *) ("verb"));

  // Read GADGET snapshot
  gadget::GadgetIO * gadget_io = new gadget::GadgetIO(in);
  int fail = gadget_io->open(in);
  if (!fail) {
    std::ostringstream stm;
    stm << gadget_io->getVersion();
    std::cerr <<"Gadget Version : "<<stm.str()<<"\n";
    glnemo::ComponentRangeVector crv = gadget_io->getCRV();

    if (verb) glnemo::ComponentRange::list(&crv);
    int index=glnemo::ComponentRange::print(&crv,select);
    if (index >= 0 ) {
      std::cout << *(gadget_io->getTime()) << " " << crv[index].range << "\n";
    }
  } 
  else {
    std::cerr << "File["<<in<<"] is not a Gadget file, aborting...\n";
    exit(1);
  }

  //   finish NEMO
  finiparam();
}
// ----------- End Of [g2info.cc] ------------------------------------
