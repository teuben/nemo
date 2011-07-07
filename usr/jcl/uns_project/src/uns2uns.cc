// ============================================================================
// Copyright Jean-Charles LAMBERT - 2011                                       
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
#include <nemo.h>       
#include <vector>
#include <map>

// ------------------------------------------------------------
// Include file
#include "uns.h"
#include "csnaptools.h"
using namespace std; 
using namespace jclut;

// NEMO parameters
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           uns input file ",
  "out=???\n          uns output file (gagdet|nemo)",
  "select=???\n       select input particles to be saved (range or component)",
  "type=???\n         type of the outfile (gadget1|gadget2|nemo)",
  "bits=\n            physicals quantities that you want to save\n",
  "times=all\n         selected time                   ",
  "first=f\n           add a trailing numbering to the first output file",
  "offset=0.01\n      +/- time offset",
  "verbose=f\n        verbose mode                    "
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL   ",
  NULL
};
const char * usage="A universal converting tool for UNS snapshots";
//------------------------------------------------------------------------------
// saveComponent
// save components according to their indexes stored in a vector
template <class T> void saveComponent(T * array, std::string typearray,
                                      const int n, 
                                      std::vector<int> &vec, 
                                      const int dim, 
                                      std::string comp, 
                                      uns::CunsOut * unsout)
{
  float * array_o = new float[vec.size()*dim]; // particles to be saved
  // loop on vector index to select particles to be save
  for (unsigned int i=0;i<vec.size();i++) {
    assert(vec[i]<n);
    for (int j=0;j<dim;j++) {
      array_o[i*dim+j] = array[vec[i]*dim+j];
    }
  }
  unsout->snapshot->setData(comp,typearray,vec.size(),array_o,false);
  delete [] array_o;
}
//------------------------------------------------------------------------------
// processFromRangeComp
// we parse range of particles from component parameters
// this method should be used with RANGE like input file, and COMPONENT output
// [select=] input string must be built like this
// ex : select=gas@0:10000+halo@200000:3000000
//  
// ALL particles MUST have been selected during CUNSIN instantiation
// read pos,vel,mass of the components
// if component exist AND it has been selected, then respecting comp's data are
// prepared to be saved
void processFromRangeComp(std::string select, uns::CunsIn * unsin,uns::CunsOut * unsout)
{
  float * pos, * vel, * mass;
  int n,nbody;
  bool ok;
  // parse select string which MUST have '@' separator
  if (select.find("@",0) == std::string::npos) { // "@" ! found
    std::cerr 
        << "\n\nprocessFromRangeComp : algorithm ERROR !!\n"
        << "you want to convert a RANGE to COMPONENT file format, but\n"
        << "your parameter [select="<<select<<"] is misformed and should be like:\n"
        << "example: select=gas@0:100+disk@200:300\n"
        << "aborting...\n\n";
    std::exit(1);
  }
  // get nbody in
  ok =unsin->snapshot->getData("nsel",&nbody);  
  // read range of particles for the component
  // and put then in a vector
  std::map <std::string, std::vector<int> > myMapOfVec;
  // build a map<component,vector of indexes>
  myMapOfVec = CSnaptools::mapStringVectorIndexes(select,nbody);
  // loop on each component of the map
  for( std::map <std::string, std::vector<int> >::iterator iter = myMapOfVec.begin(); iter != myMapOfVec.end(); ++iter ) {    
    std::string      comp = (*iter).first;   // component
    std::vector<int> vec  = (*iter).second;  // vector of indexes
    std::cerr << "["<<comp<<"] vec size = " << vec.size() << "\n";
    for (unsigned int i=0; i<vec.size();i++) {
      //std::cerr << "i="<<i<<"        -> " << vec[i] << "\n";
    }    
    // read position
    ok = unsin->snapshot->getData("all","pos" ,&n,&pos );
    if (ok) {
      saveComponent(pos,"pos",n,vec,3,comp,unsout);      
    }
    // read velocities
    ok = unsin->snapshot->getData("all","vel" ,&n,&vel );
    if (ok) {
      saveComponent(vel,"vel",n,vec,3,comp,unsout);      
    }
    // read masses
    ok = unsin->snapshot->getData("all","mass" ,&n,&mass );
    if (ok) {
      saveComponent(mass,"mass",n,vec,1,comp,unsout);      
    }
  }
}
//------------------------------------------------------------------------------
// processFromComp
// try to read all data of the components
// if component exist AND it has been selected, then respecting comp's data are
// prepared to be saved
void processFromComp(std::string comp, uns::CunsIn * unsin,uns::CunsOut * unsout)
{  
  float * pos, * vel, * mass;
  int nbody,n1,n2,n3;
  bool ok;
  // read position
  ok = unsin->snapshot->getData(comp,"pos" ,&n1,&pos );
  if (ok) {
    unsout->snapshot->setData(comp,"pos",n1,pos,false);
  }
  // read velocities
  ok = unsin->snapshot->getData(comp,"vel" ,&n2,&vel );
  if (ok) {
    unsout->snapshot->setData(comp,"vel",n2,vel,false);
  }
  // read masses
  ok = unsin->snapshot->getData(comp,"mass" ,&n3,&mass );
  if (ok) {
    unsout->snapshot->setData(comp,"mass",n3,mass,false);
  }
  // read Ids
  int * id;
  ok = unsin->snapshot->getData(comp,"id" ,&nbody,&id );
  if (ok) { 
    unsout->snapshot->setData(comp,"id",nbody,id,false);
  }
  if (comp=="all") {
    float * rho, * hsml;
    int nn;
    nbody=std::max(n1,std::max(n2,n3));
    // Try to get Rho
    ok = unsin->snapshot->getData("rho" ,&nn,&rho );
    if (ok && nbody == nn) {
      unsout->snapshot->setData(comp,"rho",nbody,rho,false);
    }  
    // Try to get Hsml
    ok = unsin->snapshot->getData("hsml" ,&nn,&hsml );
    if (ok && nbody == nn)  {
      unsout->snapshot->setData(comp,"hsml",nbody,hsml,false);  
    }
  }
  if (comp=="gas") {
    float * rho, * hsml, * u , * temp, * metal;
    // Try to get Rho
    ok = unsin->snapshot->getData(comp,"rho" ,&nbody,&rho );
    if (ok) {
      unsout->snapshot->setData(comp,"rho",nbody,rho,false);
    }  
    // Try to get Hsml
    ok = unsin->snapshot->getData(comp,"hsml" ,&nbody,&hsml );
    if (ok) 
      unsout->snapshot->setData(comp,"hsml",nbody,hsml,false);  
    // Try to get u
    ok = unsin->snapshot->getData(comp,"u" ,&nbody,&u );
    if (ok) 
      unsout->snapshot->setData(comp,"u",nbody,u,false);  
    // Try to get temperature
    ok = unsin->snapshot->getData(comp,"temp" ,&nbody,&temp );
    if (ok) 
      unsout->snapshot->setData(comp,"temp",nbody,temp,false);  
    // Try to get gas metalicity
    ok = unsin->snapshot->getData(comp,"gas_metal" ,&nbody,&metal );
    if (ok) 
      unsout->snapshot->setData(comp,"gas_metal",nbody,metal,false);  
  }
  if (comp == "stars") {
    float * age, * metal;
    // Try to get stars metalicity
    ok = unsin->snapshot->getData(comp,"stars_metal" ,&nbody,&metal );
    if (ok) 
      unsout->snapshot->setData(comp,"stars_metal",nbody,metal,false);
    // Try to get stars age
    ok = unsin->snapshot->getData(comp,"age" ,&nbody,&age );
    if (ok) 
      unsout->snapshot->setData(comp,"age",nbody,age,false);
  }
  
}
//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
int main(int argc, char ** argv )
{
  // global bariable
  std::string typein,file_structin, file_structout;
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv)); 
  if (argc) {;} // remove compiler warning :)
  
  // Get input parameters
  std::string simname   = (getparam ((char *) "in"      ));
  std::string outname   = (getparam ((char *) "out"     ));
  std::string typeout   = (getparam ((char *) "type"    ));
  std::string select    = (getparam ((char *) "select"  ));
  std::string bits      = (getparam ((char *) "bits"    ));
  std::string select_t  = (getparam ((char *) "times"   ));
  bool        first     = (getbparam((char *) "first"   ));
  float       offset    = (getdparam((char *) "offset"  ));
  bool        verbose   = (getbparam((char *) "verbose" ));
  
  bool one_file=false;
  bool stop=false;
  bool special_nemo=false;
  if (outname=="-" || outname==".") {
    if (typeout!="nemo") {
      std::cerr << "You cannot request out=["<<outname<<"] with type=["<<typeout<<"\n";
      std::exit(1);
    }
    special_nemo=true;
  }
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
  
  std::string select_c = select; // default component selection
  if (select.find("@",0) != std::string::npos) { // "@" found
    std::cerr << "Found [@] character in select string\n";
    // we assume that structure file is in=range and out=component
    select_c = "all";
  }
  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * unsin = new uns::CunsIn(simname,select_c,select_t,verbose);
  if (unsin->isValid()) { // input file is known by UNS lib        
    int cpt=0;
    while(unsin->snapshot->nextFrame(bits)&&!stop) { // there is a new frame
      typein = unsin->snapshot->getInterfaceType();
      file_structin = unsin->snapshot->getFileStructure();
      std::cerr << "Input file type      :"<<typein<<"\n";
      std::cerr << "Input file structure :"<<file_structin<<"\n";
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
          unsout = new uns::CunsOut(out_name,typeout,verbose);      
          file_structout = unsout->snapshot->getFileStructure();
        } else {
          if (first_out) {
            first_out = false;
            // instantiate only once unsout, because outname="-"
            unsout = new uns::CunsOut(out_name,typeout,verbose);
            file_structout = unsout->snapshot->getFileStructure();            
          }
        }
        std::cerr << "output filename=["<<out_name<<"]\n";
        std::cerr << "Input  file structure  ="<< file_structin <<"\n";
        std::cerr << "Output file structure  ="<< file_structout <<"\n";
        // save time
        unsout->snapshot->setData("time",time);
        // processing
        if (file_structout=="range") {
          // Whatever input file structure but output file structure is RANGE
          // range= field must be up
          processFromComp("all",unsin,unsout); // only all particles selected
        } else {
          if (file_structin == "range"      && 
              file_structout== "component") {
            std::cerr << "processFromRangeComp\n";
            processFromRangeComp(select,unsin,unsout);
          } else {
            if (file_structin == "component" && 
                file_structout== "component") {
              std::cerr << "processFromComp\n";
              processFromComp("gas"  ,unsin,unsout);
              processFromComp("disk" ,unsin,unsout);
              processFromComp("stars",unsin,unsout);
              processFromComp("bulge",unsin,unsout);
              processFromComp("halo" ,unsin,unsout);
              processFromComp("bndry",unsin,unsout);
            }
          }
        }
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
// ----------- End Of [uns2uns.cc] ------------------------------------
