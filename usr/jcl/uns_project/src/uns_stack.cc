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
#include <cstdio>  
#include <cstdlib> 
#include <assert.h>
#include <cmath>
#include <nemo.h>
#include <iomanip>
#include <vector>

// ------------------------------------------------------------
// Include file
#include "uns.h"
#include "csnaptools.h"
// ------------------------------------------------------------
// Nemo variable
const char * defv[] = {
  "in1=???\n		      UNS input snapshot",  
  "in2=???\n		      UNS input snapshot",
  "out=???\n	              output snapshot",
  "deltar=0.0,0.0,0.0\n	      position of in1 w.r.t. in2",
  "deltav=0.0,0.0,0.0\n	      velocity of in1 w.r.t. in2",
  "shift=1\n                  Shift over 1st or 2nd one?",
  "zerocm=f\n                 Centering On Mass after stacking?",
  "verbose=f\n                verbose on/off",
  "VERSION=1.O\n              compiled on <"__DATE__"> JCL  ",
  NULL,
};
const char * usage="stack two UNS systems on top of each other";
using namespace std;
using namespace uns;
using namespace jclut;
// ------------------------------------------------------------
// addArray
bool addArray(std::string comp, std::string name, int dim, CunsIn * uns1,CunsIn * uns2, CunsOut * unsout, bool verbose)
{
  bool status=false;
  float * d=NULL, * d1=NULL, * d2=NULL;
  int n,n1=0,n2=0;
  bool ok1 = uns1->snapshot->getData(comp,name ,&n1,&d1);
  bool ok2 = uns2->snapshot->getData(comp,name ,&n2,&d2);
  n=0;
  if (ok1 || ok2) {
    std::cerr << "n1="<<n1<<" n2="<<n2;
    n = n1+n2;
    assert(n>0);
    d = new float[n*dim];
    if (verbose)
      std::cerr << "* INFO * addArray getData["<<comp<<","<<name<<"] n1="<<n1<<" n2="<<n2<<"\n";
    memcpy(d         ,d1,sizeof(float)*n1*dim); // copy first array
    memcpy(d+(n1*dim),d2,sizeof(float)*n2*dim); // copy second array
    // prepare data to be saved
    unsout->snapshot->setData(comp,name,n,d,false);
    delete [] d;
    status=true;
  } else {
    if (ok1 | ok2) {
      if (!ok1)
        std::cerr << "* ERROR * addArray on file<"<<uns1->snapshot->getFileName()<<">\n"
            "getData["<<comp<<","<<name<<"] failed, aborting....\n";
      if (!ok2)
        std::cerr << "* ERROR * addArray on file<"<<uns2->snapshot->getFileName()<<">\n"
            "getData["<<comp<<","<<name<<"] failed, aborting....\n";
      std::exit(1);
    } else {
      if (verbose)
        std::cerr << "* WARNING * addArray getData["<<comp<<","<<name<<"] failed, skipping...\n";
    }
  }
  return status;
}
// ------------------------------------------------------------
// addComponent
void addComponent(std::string comp, CunsIn * uns1,CunsIn * uns2, 
                  CunsOut * unsout, bool verbose)
{
  addArray(comp,"pos" ,3,uns1,uns2,unsout,verbose);
  addArray(comp,"vel" ,3,uns1,uns2,unsout,verbose);
  addArray(comp,"mass",1,uns1,uns2,unsout,verbose);
  if (comp == "gas") {
    addArray(comp,"u"    ,1,uns1,uns2,unsout,verbose);
    addArray(comp,"hsml" ,1,uns1,uns2,unsout,verbose);
    addArray(comp,"rho"  ,1,uns1,uns2,unsout,verbose);
  }
}
// ------------------------------------------------------------
// process
void process(CunsIn * uns1,CunsIn * uns2, char * out, char * dr, char * dv, bool com, int shift, bool verbose)
{
  // convert string to vector
  std::vector<float> deltar=CSnaptools::stringToVector<float>(dr,3,0.0);
  std::vector<float> deltav=CSnaptools::stringToVector<float>(dv,3,0.0);

  if (shift==1) {
    uns1->snapshot->shift("pos",deltar[0],deltar[1],deltar[2]);
    uns1->snapshot->shift("vel",deltav[0],deltav[1],deltav[2]);
  } else {
    uns2->snapshot->shift("pos",deltar[0],deltar[1],deltar[2]);
    uns2->snapshot->shift("vel",deltav[0],deltav[1],deltav[2]);
  }
  // -----------------------------------------------
  // Instantiate a UNS output snapshot (for writing)
  uns::CunsOut * unsout = new uns::CunsOut(out,uns1->snapshot->getInterfaceType(),verbose);   
  
  if (uns1->snapshot->getInterfaceType() == "Nemo") {
    addComponent("all",uns1,uns2,unsout,verbose);
  }
  else {
    addComponent("gas"  ,uns1,uns2,unsout,verbose);
    addComponent("halo" ,uns1,uns2,unsout,verbose);
    addComponent("disk" ,uns1,uns2,unsout,verbose);
    addComponent("bulge",uns1,uns2,unsout,verbose);
    addComponent("stars",uns1,uns2,unsout,verbose);
    addComponent("bndry",uns1,uns2,unsout,verbose);
  }
  if (com)
    std::vector<double> com=unsout->snapshot->moveToCom(); // shift to COM
  unsout->snapshot->save();  // save file
}

// ------------------------------------------------------------
// main program
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  // get input parameters
  char * in1     = getparam ((char *) "in1"      );
  char * in2     = getparam ((char *) "in2"      );
  char * out     = getparam ((char *) "out"      );
  char * dr      = getparam ((char *) "deltar"   );
  char * dv      = getparam ((char *) "deltav"   );
  bool   com     = getbparam((char *) "zerocm"   );
  int    shift   = getiparam((char *) "shift"    );  
  bool   verbose = getbparam((char *) "verbose"  );  
  
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns1 = new uns::CunsIn(in1,"all","all",verbose);
  
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns2 = new uns::CunsIn(in2,"all","all",verbose);
  
  // some checking
  if (!uns1->isValid()) {
    std::cerr << "File ["<<in1<<"] is not an UNS known file format\n";
    std::exit(1);
  }
  if (!uns2->isValid()) {
    std::cerr << "File ["<<in2<<"] is not an UNS known file format\n";
    std::exit(1);
  }
  // read first time
  bool ok1=uns1->snapshot->nextFrame();
  bool ok2=uns2->snapshot->nextFrame();
  
  if (ok1&&ok2) {        
    // snapshot type must be identical
    if (uns1->snapshot->getInterfaceType() != uns2->snapshot->getInterfaceType()) {
      std::cerr << "UNS files types are not identical, aborting....\n";
      std::exit(1);
    }
    process(uns1,uns2,out,dr,dv,com,shift,verbose);
  } else {
    std::cerr << "Can't read nextFrame ....\n";
  }
  
  //   finish NEMO
  finiparam();
}
