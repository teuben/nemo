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
#include <iomanip>
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
  "type=gadget2\n     type of the output file  (gadget2|nemo)",
  "time=all\n         selected time                   ",
  "verbose=f\n        verbose mode                    "
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL   ",
  NULL
};
const char * usage="test unsio library";

//------------------------------------------------------------------------------
// processComponent
// read pos,vel,mass of the components
// if component exist AND it has been selected, then respecting comp's data are
// prepared to be saved
void processComponent(std::string comp, uns::CunsIn * uns,uns::CunsOut * unsout)
{
  float * pos, * vel, * mass;
  int * id;
  int n1,n2,n3,n4;
  bool ok1,ok2,ok3,ok4;


  ok1 = uns->snapshot->getData(comp,"pos" ,&n1,&pos );
  ok2 = uns->snapshot->getData(comp,"vel" ,&n2,&vel );
  ok3 = uns->snapshot->getData(comp,"mass",&n3,&mass);
  
  if (ok1 && ok2 && ok3) {
    assert(n1==n2);
    assert(n1==n3);
    std::cerr << "--> "<< std::left << std::setfill('.')<<
        std::setw(8) << comp << ":" << std::setfill(' ')<<
        std::right   << std::setw(10) << n1 <<"\n";
    unsout->snapshot->setData(comp,n1,mass,pos,vel,false);
    
  }
  ok4 = uns->snapshot->getData(comp,"id"  ,&n4,&id  );
  if (ok4) unsout->snapshot->setData(comp,"id",n4,id,false);
  float * rho;
  ok4 = uns->snapshot->getData(comp,"rho"  ,&n4,&rho  );
  if (ok4) unsout->snapshot->setData(comp,"rho",n4,rho,false);
  float * hsml;
  ok4 = uns->snapshot->getData(comp,"hsml"  ,&n4,&hsml  );
  if (ok4) unsout->snapshot->setData(comp,"hsml",n4,hsml,false);
  
}
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
  uns::CunsIn * unsin = new uns::CunsIn(simname,select_c,select_t,verbose);
  
  if (unsin->isValid()) { // input file is known by UNS lib        
    int cpt=0;
    while(unsin->snapshot->nextFrame()) { // there is a new frame
      std::string itype = unsin->snapshot->getInterfaceType();
      std::cerr << "Input file is of type :"<<itype<<"\n";
      bool ok;
      int cnbody,cfirst,clast,nbody;      
      float time;
      // get the input number of bodies according to the selection
      ok=unsin->snapshot->getData("nsel",&nbody);
      // get the simulation time
      ok=unsin->snapshot->getData("time",&time);
      //
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

      // processing
      processComponent("halo" ,unsin,unsout);
      processComponent("disk" ,unsin,unsout);
      processComponent("bulge",unsin,unsout);
      processComponent("stars",unsin,unsout);
      processComponent("bndry",unsin,unsout);
      processComponent("gas"  ,unsin,unsout);
      processComponent("all"  ,unsin,unsout); // only all particles selected
      
      // according to user's input request ("select" parameter)
      // check if gas component exist from input snapshot     
      bool is_gas=unsin->snapshot->getRangeSelect("gas",&cnbody,&cfirst,&clast);
      if (( is_gas                    &&         // gas exist
          type      == "gadget2")     ||         // AND outype is gadget2
          (select_c == "gas" && type=="nemo" &&  // OR only gas and nemoout selected
           is_gas)) {                            // AND gas exist
        // Density
        float * rho;
        ok=unsin->snapshot->getData("rho",&cnbody,&rho);
        if (ok)  unsout->snapshot->setData("rho",cnbody,rho,false);

        // Hydro Smooth Length
        float * hsml;
        ok=unsin->snapshot->getData("hsml",&cnbody,&hsml);
        if (ok)   unsout->snapshot->setData("hsml",cnbody,hsml,false);
        
        // internal energy
        float * u;
        ok=unsin->snapshot->getData("u",&cnbody,&u);
        if (ok)  unsout->snapshot->setData("u",cnbody,u,false);
      }
      // for NEMO out only
      if (itype=="Nemo") {
        //processComponent("all" ,unsin,unsout);
      }      
      // save snapshot
      unsout->snapshot->save();
      delete unsout; // remove object      
    }
  } else {
    std::cerr << "Unknown UNS file format["<<simname<<"]\n";
  }
  delete unsin;
  //   finish NEMO
  finiparam();
}
// ----------- End Of [unsio_demo.cc] ------------------------------------


