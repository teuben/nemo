// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#include <iostream>                                   // C++ I/O     
#include <fstream>                                    // C++ file I/O
#include <sstream>
#include <cstdio>                    // changed from stdio.h  WD, Sep 2008
#include <cstdlib>                   // changed from stdlib.h WD, Sep 2008
#include <assert.h>

#include "uns.h"
#define DEBUG 0
#include "unsdebug.h"

#define _vectmath_h // put this statement to avoid conflict with C++ vector class
extern "C" {
#include <nemo.h>                                     // NEMO basics
  int io_nemo(const char *, const char *,...);
}

using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
// NEMO parameters
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           input file (gadget|nemo)          ",
  "time=all\n         selected time",
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL  ",
  NULL
};
const char * usage="print out rho and temperature for the gas component";

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  if (argc) {;} // remove compiler warning :)
  // Get parameters
  char * simname   = getparam((char *) "in"    );
  char * select_t  = getparam((char *) "time"  );
  // instantiate a new uns object
  //s::Cuns * uns = new uns::Cuns(simname,select_c,select_t);
  uns::CunsIn * uns = new uns::CunsIn(simname,"gas",select_t);
  if (uns->isValid()) {
    while(uns->snapshot->nextFrame()) {
        bool ok;
	int n,nbody;
        float time;
        // get the input number of bodies according to the selection
        ok=uns->snapshot->getData("nsel",&nbody);
        // get the simulation time
        ok=uns->snapshot->getData("time",&time);
	std::cerr << "nbody=" << nbody << " time="<<time <<"\n";
        // rho
        float * rho=NULL;
        ok=uns->snapshot->getData("rho",&n,&rho);
        std::cerr << "n=" << n <<"\n";
	assert(n== nbody);
        // temp
	float * temp=NULL;
        ok=uns->snapshot->getData("temp",&n,&rho);
        std::cerr << "n=" << n <<"\n";
	assert(n== nbody);
	for (int i=0; i<n; i++) {
	  std::cout << rho[i] << " " << temp[i] << "\n";
	}

    }
  } else {
    std::cerr << "Unknown UNS file format["<<simname<<"]\n";
  }
  //   finish NEMO
  finiparam();
}
// ----------- End Of [stress_io_nemo.cc] ------------------------------------
