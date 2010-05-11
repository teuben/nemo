// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010                                       
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
#include <cstdio>                   
#include <cstdlib>                  
#include <assert.h>
#include <cstdio>
#include "uns.h"

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
  "in=???\n           input file (gadget|nemo)        ",
  "out=???\n          output file                     ",
  "select=???\n       component selected (disk,stars,halo,gas,range)",
  "type=gadget2\n     type of the output file         ",
  "time=all\n         selected time                   ",
  "verbose=f\n        verbose mode                    "
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL   ",
  NULL
};
const char * usage="test unsio library";

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  if (argc) {;} // remove compiler warning :)
  
  // Get input parameters
  std::string simname (getparam ((char *) "in"      ));
  std::string outname (getparam ((char *) "out"     ));
  std::string type    (getparam ((char *) "type"    ));
  std::string select_c(getparam ((char *) "select"  ));
  std::string select_t(getparam ((char *) "time"    ));
  bool        verbose =getbparam((char *) "verbose" );
  
  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns = new uns::CunsIn(simname,select_c,select_t,verbose);
  
  if (uns->isValid()) { // input file is known by UNS lib        
    int cpt=0;
    while(uns->snapshot->nextFrame()) { // there is a new frame
      std::cerr << "Input file is of type :"<<uns->snapshot->getInterfaceType()<<"\n";
      bool ok;
      int cnbody,cfirst,clast,nbody;      
      float * pos, * vel, * mass, time;
      // get the input number of bodies according to the selection
      ok=uns->snapshot->getData("nsel",&nbody);
      // get the simulation time
      ok=uns->snapshot->getData("time",&time);
      // get POS from input snapshot
      ok=uns->snapshot->getData("pos" ,&cnbody,&pos);
      // get VEL from input snapshot
      ok=uns->snapshot->getData("vel" ,&cnbody,&vel);
      // get MASS from input snapshot
      ok=uns->snapshot->getData("mass",&cnbody,&mass);
      
      std::cerr << "nbody=" << nbody << " time="<<time <<"\n";

      // OUTPUT operations
      // create an output filename : basename +  integer
      // example : myoutput.0 myoutput.1 ...... etc
      stringstream number;
      number << cpt++;
      std::string out_name = std::string(outname)+"."+number.str();;
      
      // -----------------------------------------------
      // Instantiate a UNS output snapshot (for writing)
      uns::CunsOut * unsout = new uns::CunsOut(out_name,type,verbose);
      
      // save time
      unsout->snapshot->setData("time",time);

      // according to user's input request ("select" parameter)
      // check if halo component exist from input snapshot
      if (uns->snapshot->getRangeSelect("halo",&cnbody,&cfirst,&clast)) {
        fprintf(stderr,"Halo :nbody %d first = %d last = %d\n",cnbody,cfirst,clast);
        unsout->snapshot->setData("halo",cnbody,mass+cfirst,pos+cfirst*3,vel+cfirst*3,false);
      }
      // according to user's input request ("select" parameter)
      // check if disk component exist from input snapshot
      if (uns->snapshot->getRangeSelect("disk",&cnbody,&cfirst,&clast)) {
        fprintf(stderr,"disk :nbody %d first = %d last = %d\n",cnbody,cfirst,clast);
        for (int i=0; i<cnbody; i++) {
          //fprintf(stderr,"disk mass = %f\n",mass[cfirst+i]);
        }
        unsout->snapshot->setData("disk",cnbody,&mass[cfirst],&pos[cfirst*3],&vel[cfirst*3],false);
      }
      // according to user's input request ("select" parameter)
      // check if bulge component exist from input snapshot
      if (uns->snapshot->getRangeSelect("bulge",&cnbody,&cfirst,&clast)) {
        fprintf(stderr,"bulge :nbody %d first = %d last = %d\n",cnbody,cfirst,clast);

        unsout->snapshot->setData("bulge",cnbody,mass+cfirst,pos+cfirst*3,vel+cfirst*3,false);
      }
      // according to user's input request ("select" parameter)
      // check if stars component exist from input snapshot      
      if (uns->snapshot->getRangeSelect("stars",&cnbody,&cfirst,&clast)) {
        fprintf(stderr,"stars :nbody %d first = %d last = %d\n",cnbody,cfirst,clast);
        unsout->snapshot->setData("stars",cnbody,mass+cfirst,pos+cfirst*3,vel+cfirst*3,false);
      }

      // according to user's input request ("select" parameter)
      // check if gas component exist from input snapshot      
      if ((uns->snapshot->getRangeSelect("gas",&cnbody,&cfirst,&clast) && // gas exist
          type     == "gadget2") ||                                       // AND outype is gadget2
          (select_c == "gas" && type=="nemo" &&                           // OR only gas/nemo selected
           uns->snapshot->getRangeSelect("gas",&cnbody,&cfirst,&clast))) {// AND gas exist
        fprintf(stderr,"Gas :nbody %d first = %d last = %d\n",cnbody,cfirst,clast);
        unsout->snapshot->setData("gas",cnbody,mass+cfirst,pos+cfirst*3,vel+cfirst*3,false);

        // Density
        float * rho;
        ok=uns->snapshot->getData("rho",&cnbody,&rho);
        if (ok)  unsout->snapshot->setData("rho",cnbody,rho,false);

        // Hydro Smooth Length
        float * hsml;
        ok=uns->snapshot->getData("hsml",&cnbody,&hsml);
        if (ok)   unsout->snapshot->setData("hsml",cnbody,hsml,false);
        
        // internal energy
        float * u;
        ok=uns->snapshot->getData("u",&cnbody,&u);
        if (ok)  unsout->snapshot->setData("u",cnbody,u,false);
      }
      // for NEMO out only
      unsout->snapshot->setData("all",nbody,mass,pos,vel,false);
      
      // save snapshot
      unsout->snapshot->save();
      delete unsout; // remove object      
    }
  } else {
    std::cerr << "Unknown UNS file format["<<simname<<"]\n";
  }
  delete uns;
  //   finish NEMO
  finiparam();
}
// ----------- End Of [unsio_demo.cc] ------------------------------------

