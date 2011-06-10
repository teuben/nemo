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
#include <nemo.h>       

using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
// NEMO parameters
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           input file (gadget|nemo)        ",
  "out=???\n          output file                     ",
  "select=???\n       component selected (disk,stars,halo,gas,range,all)",
  "times=all\n         selected time                   ",
  "first=f\n           add a trailing numbering to the first output file",
  "offset=0.01\n      +/- time offset",
  "verbose=f\n        verbose mode                    "
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL   ",
  NULL
};
const char * usage="Convert an UNS file to Nemo file format";

//------------------------------------------------------------------------------
// processComponent
// read pos,vel,mass of the components
// if component exist AND it has been selected, then respecting comp's data are
// prepared to be saved
void processComponent(std::string select,std::string comp, uns::CunsIn * uns,uns::CunsOut * unsout)
{
  float * pos, * vel, * mass;
  int n1,n2,n3;
  bool ok1,ok2,ok3;
  
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
 
  float * rho, * hsml;
  bool ok;
  int nn;
  // Try to get Rho
  ok = uns->snapshot->getData("rho" ,&nn,&rho );
  if (ok && n1 == nn) {
    unsout->snapshot->setData("rho",n1,rho,false);
  }  
  // Try to get Hsml
  ok = uns->snapshot->getData("hsml" ,&nn,&hsml );
  if (ok && n1 == nn) 
    unsout->snapshot->setData("hsml",n1,hsml,false);  
  // Try to get Ids
  int * id;
  ok = uns->snapshot->getData("id" ,&nn,&id );
  if (ok && n1 == nn) 
    unsout->snapshot->setData("id",n1,id,false);  
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
  std::string select_c(getparam ((char *) "select"  ));
  std::string select_t(getparam ((char *) "times"    ));
  bool first=getbparam((char *) "first"     );
  float       offset= getdparam((char *) "offset"   );
  bool        verbose =getbparam((char *) "verbose" );
    
  bool one_file=false;
  bool stop=false;
  bool special_nemo=false;
  if (outname=="-" || outname==".") special_nemo=true;
  // in case of an input simulation from the database
  // and with just one time requested,
  // we create a range of time to speedup the searching
  if (select_t!="all" && select_t.find(":",0)==std::string::npos) {
    float match_time;
    stringstream ss;
    ss << select_t;
    ss >> match_time; // convert string time to float
    ss.str(std::string()); // empty stringstream
    ss.clear();            // empty stringstream (mandatory after >>)
    ss << match_time-offset<<":"<<match_time+offset;
    select_t = ss.str();
    one_file=true;
    std::cerr << "Modified selected time =["<<select_t<<"]\n";
  }

  uns::CunsOut * unsout=NULL; // out object
  bool first_out=true;
  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * unsin = new uns::CunsIn(simname,select_c,select_t,verbose);
  
  if (unsin->isValid()) { // input file is known by UNS lib        
    int cpt=0;
    while(unsin->snapshot->nextFrame("mxvRHI")&&!stop) { // there is a new frame
      std::string itype = unsin->snapshot->getInterfaceType();
      std::cerr << "Input file is of type :"<<itype<<"\n";
      bool ok;
      int nbody;      
      float time;
      // get the input number of bodies according to the selection
      ok =unsin->snapshot->getData("nsel",&nbody);
      // get the simulation time
      ok=unsin->snapshot->getData("time",&time);
      //      
      std::cerr << "nbody=" << nbody << " time="<<time <<"\n";
      if (nbody>0) { // there are particles
        // OUTPUT operations
        // create an output filename : basename +  integer
        // example : myoutput.0 myoutput.1 ...... etc
        stringstream number;
        number << cpt++;
        std::string out_name=std::string(outname);;
        if (! special_nemo) { // ! standard output && ! "."
          if (one_file || (cpt==1 && !first)) {
            out_name=std::string(outname);
            if (one_file) stop = true; // do not continue
          } else {
            stringstream ss;
            ss << std::string(outname) << "." << setw(5) << setfill('0') << number.str();
            //out_name=std::string(outname)+"."+number.str();
            out_name=ss.str();
          }
          // create a new UNS out object
          unsout = new uns::CunsOut(out_name,"nemo",verbose);      
        } else {
          if (first_out) {
            first_out = false;
            // instantiate only once unsout, because outname="-"
            unsout = new uns::CunsOut(out_name,"nemo",verbose);
          }
        }
        std::cerr << "output filename=["<<out_name<<"]\n";
        
        // save time
        unsout->snapshot->setData("time",time);
        // processing
        processComponent(select_c,"all"  ,unsin,unsout); // only all particles selected
        // save snapshot
        unsout->snapshot->save();
        
        if (!special_nemo) {
          delete unsout; // remove object      
        }
      }
    }
  } else {
    std::cerr << "Unknown UNS file format["<<simname<<"]\n";
  }
  delete unsin;
  //   finish NEMO
  finiparam();
}
// ----------- End Of [unsio_demo.cc] ------------------------------------



