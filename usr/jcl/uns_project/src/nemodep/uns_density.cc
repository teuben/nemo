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
//#include <cstdio>                   
//#include <cstdlib>                  
#include <assert.h>
#include "uns.h"
#include <nemo.h>
#include "cfalcon.h"

using namespace std; // prevent writing statment like 'std::cerr'
using namespace jclut;

//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
// NEMO parameters
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           UNS input file (gadget|nemo)                       ",
  "out=???\n          output file                                        ",
  "select=???\n       component selected (disk,stars,halo,gas,range)     ",
  "m=0\n              0: Ferrers  method, 1:hackdens method              ",
  "K=32\n             number of neighbours                               ",
  "N=1\n              order of Ferrers kernel                            ",
  "ncrit=\n           ncrit value, if empty=max(1,K/4)                   ",
  "time=all\n         selected time                                      ",
  "first=f\n           add a trailing numbering to the first output file ",
  "offset=0.01\n      +/- time offset                                    ",
  "verbose=f\n        verbose mode                                       ",
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL                      ",
  NULL
};
const char * usage="estimate mass density based on distance to Kth neighbour, add distance to the kth neighbours";

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  if (argc) {;} // remove compiler warning :)
  
  // Get input parameters
  std::string simname  (getparam ((char *) "in"      ));
  std::string outname  (getparam ((char *) "out"     ));
  std::string select_c (getparam ((char *) "select"  ));
  const unsigned K     (getiparam((char *) "K"       ));
  const int      N     (getiparam((char *) "N"       ));
  int ncrit=std::max(1,(int)K/4);
  if (hasvalue((char *) "ncrit" )) {
    ncrit = getiparam((char *) "ncrit");
  }
  const unsigned method(getiparam((char *) "m"       ));
  std::string select_t (getparam ((char *) "time"    ));
  bool first=           getbparam((char *) "first"    );
  float       offset=   getdparam((char *) "offset"   );
  bool        verbose = getbparam((char *) "verbose"  );
  
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
  
  // 
  uns::CunsOut * unsout=NULL; // UNS out object
  bool first_out=true;
  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns = new uns::CunsIn(simname,select_c,select_t,verbose);
  
  if (uns->isValid()) { // input file is known by UNS lib        
    int cpt=0;
    while(uns->snapshot->nextFrame("maxvpI")&&!stop) { // there is a new frame
      std::cerr << "Input file is of type :"<<uns->snapshot->getInterfaceType()<<"\n";
      bool ok;
      int cnbody,nbody;      
      float * pos, * vel, * mass, * acc, time;
      // get the input number of bodies according to the selection
      ok=uns->snapshot->getData("nsel",&nbody);
      assert(ok==true);
      // get the simulation time
      ok=uns->snapshot->getData("time",&time);
      // get POS from input snapshot
      ok=uns->snapshot->getData("pos" ,&cnbody,&pos);
      assert(ok==true);
      // get VEL from input snapshot
      ok=uns->snapshot->getData("vel" ,&cnbody,&vel);
      // get MASS from input snapshot
      ok=uns->snapshot->getData("mass",&cnbody,&mass);
      assert(ok==true);
      // get ACC from input snapshot
      ok=uns->snapshot->getData("acc",&cnbody,&acc);
      std::cerr << "nbody=" << nbody << " time="<<time <<"\n";
      // OUTPUT operations
      // create an output filename : basename +  integer
      // example : myoutput.0 myoutput.1 ...... etc
      if (nbody >0) {
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
        // Instantiate a density object  
        CDensity * density = new CDensity(nbody,pos,mass);     
        density->compute(method,K,N,ncrit); // estimate density
        
        // get back data
        float * rho  = density->getRho();
        float * hsml = density->getHsml();
        // feed up UNS out arrays
        unsout->snapshot->setData("time",time);
        //unsout->snapshot->setData("all" ,nbody,mass,pos,vel,false);
        if (pos) {
          unsout->snapshot->setData("pos",nbody,pos,false);
        }
        if (vel) {
          unsout->snapshot->setData("vel",nbody,vel,false);
        }
        if (mass) {
          unsout->snapshot->setData("mass",nbody,mass,false);
        }
        if (acc) {
          unsout->snapshot->setData("acc",nbody,acc,false);
        }
        unsout->snapshot->setData("rho" ,nbody,rho ,false);
        unsout->snapshot->setData("hsml",nbody,hsml,false);
        // add pot (for manu)
        float * pot=NULL;
        ok=uns->snapshot->getData("pot" ,&cnbody,&pot);
        if (ok) {
          unsout->snapshot->setData("pot",cnbody,pot,false);
        }
        // Try to get Ids
        int * id,nn;
        ok = uns->snapshot->getData("id" ,&nn,&id );
        if (ok && nbody == nn) 
          unsout->snapshot->setData("id",nn,id,false);  
        
        // save snapshot
        unsout->snapshot->save();
        if (!special_nemo) {
          delete unsout; // remove object      
        }
        delete density;
      }
    }
  }
}
