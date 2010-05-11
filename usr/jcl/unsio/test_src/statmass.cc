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
    "select=all\n       component selected (disk,stars,halo,gas,range)",
    "time=all\n         selected time",
    "VERSION=1.0\n       compiled on <"__DATE__"> JCL  ",
    NULL
};
const char * usage="Compute total mass";

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
    //   start  NEMO
    initparam(const_cast<char**>(argv),const_cast<char**>(defv));
    if (argc) {;} // remove compiler warning :)
    // Get parameters
    char * simname   = getparam((char *) "in"    );
    char * select_c  = getparam((char *) "select");
    char * select_t  = getparam((char *) "time"  );

    // instantiate a new uns object
    //s::Cuns * uns = new uns::Cuns(simname,select_c,select_t);
    uns::CunsIn * uns = new uns::CunsIn(simname,select_c,select_t);
    if (uns->isValid()) {
        while(uns->snapshot->nextFrame()) {
          bool ok;
          int nbody;
          float time;
          // get the input number of bodies according to the selection
          ok=uns->snapshot->getData("nsel",&nbody);
          // get the simulation time
          ok=uns->snapshot->getData("time",&time);
          std::cerr << "nbody=" << nbody << " time="<<time <<"\n";
          
          // get mass
          float * m;
          ok=uns->snapshot->getData("mass",&nbody,&m);
          std::cerr << "nbody=" << nbody << " time="<< time <<"\n";
          double masstot=0;
          for (int i=0; i<nbody; i++) {
            masstot += m[i];
          }
          std::cout << time << " " << masstot << "\n";
        }
    }
    else {
        std::cerr << "Unknown UNS file format["<<simname<<"]\n";
    }
    //   finish NEMO
    finiparam();
}
// ----------- End Of [stress_io_nemo.cc] ------------------------------------
