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
#include <cmath>
#include <nemo.h>
#include <iomanip>

// ------------------------------------------------------------
// Include file
#include "uns.h"

// ------------------------------------------------------------
// Nemo variable
const char * defv[] = {
  "in=???\n		      UNS input snapshot",  
  "select=???\n               select particles (range, or component name)\n"
  "                   component: gas,halo,disk,bulge,stars,bndry",
  "display=t\n                display array's content (t|f)",
  "maxline=2\n                max lines per components",  
  "times=all\n		      selected time",
  "verbose=f\n                verbose on/off",
  "VERSION=1.O\n              compiled on <"__DATE__"> JCL  ",
  NULL,
};
const char * usage="Print information about an UNS file";
using namespace std;
void displayInfo(bool display,int maxlines, std::string comp, uns::CunsIn * uns);
template <class T> void displayFormat(int maxlines,std::string text, T * array, int dim, int size, int np);
// ------------------------------------------------------------
//  displayInfo
void displayInfo(bool display,int maxlines, std::string comp, uns::CunsIn * uns)
{
  float * pos, * vel, * mass;
  int * id;
  int n1,n2,n3,n4;
  bool ok1,ok2,ok3,ok4;

  ok1 = uns->snapshot->getData(comp,"pos" ,&n1,&pos );
  ok2 = uns->snapshot->getData(comp,"vel" ,&n2,&vel );
  ok3 = uns->snapshot->getData(comp,"mass",&n3,&mass);
  ok4 = uns->snapshot->getData(comp,"id"  ,&n4,&id);
  if (ok1 || ok2 || ok3) {
    int nbody=max(max(n1,n2),n3);
    std::cerr << setw(50) << setfill('=') << ""<<"\n";
    std::cerr<< setfill(' ');
    std::cerr << left<< setw(8) << comp << ":" << setw(9) << right << nbody << "\n";
  }
  if (ok3 && display) {
    displayFormat(maxlines,"mass[1] = ",mass,1,n3, 3);
  }  
  if (ok1 && display) {
    displayFormat(maxlines,"pos [3] = ",pos ,3,n1, 1);
  }
  if (ok2 && display) {
    displayFormat(maxlines,"vel [3] = ",vel ,3,n2, 1);
  }
  if (ok4 && display) {
    displayFormat(maxlines,"id  [1] = ",id  ,1,n4, 3);
  }  
  if (comp == "gas" && ok1 && ok2 && ok3) {
    float * rho, * u, * hsml;
    ok1 = uns->snapshot->getData(comp,"rho" ,&n1,&rho );
    ok2 = uns->snapshot->getData(comp,"u"   ,&n2,&u );
    ok3 = uns->snapshot->getData(comp,"hsml",&n3,&hsml);
    if (ok1 && display) {
      displayFormat(maxlines,"rho [1] = ",rho ,1,n1, 3);
    }
    if (ok2 && display) {
      displayFormat(maxlines,"u   [1] = ",u   ,1,n2, 3);
    }
    if (ok3 && display) {
      displayFormat(maxlines,"hsml[1] = ",hsml,1,n3, 3);
    }
    
  }
}
// ------------------------------------------------------------
// displayFormat
template <class T>  void displayFormat(int maxlines,std::string text, T * array, int dim, int size, int np)
{
  std::cerr << scientific << left << setw(11) << text;
  // First line
  for (int k=0;k<std::min(size,np);k++) {
    for (int j=0;j<dim;j++) {
      std::cerr << array[k*dim+j] << " ";
    }
  }
  std::cerr << "\n";
  // other lines
  for (int i=1; i<std::min(maxlines,size); i+=np) {
    std::cerr << left << setw(11) << "";
    for (int k=0;k<std::min(size,np);k++) {
      for (int j=0;j<dim;j++) {
        std::cerr << array[(k+(i*np))*dim+j] << " ";
      }
    }
    std::cerr << "\n";
  }  
  std::cerr << left << setw(11) << "" << ". . .\n";
}

// ------------------------------------------------------------
// main program
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  // get input parameters
  char * simname     = getparam ((char *) "in"      );
  char * select_c    = getparam ((char *) "select"  );
  bool   display     = getbparam((char *) "display" );
  int    maxlines    = getiparam((char *) "maxline" );
  char * select_t    = getparam ((char *) "times"   );
  bool   verbose     = getbparam((char *) "verbose" );  
  
  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns = new uns::CunsIn(simname,select_c,select_t,verbose);
  
  if (uns->isValid()) { // input file is known by UNS lib        
    while(uns->snapshot->nextFrame()) { // there is a new frame
      std::string stype = uns->snapshot->getInterfaceType();
      std::cerr << setw(50) << setfill('*') << ""<<"\n";
      std::cerr << "File name : "<<uns->snapshot->getFileName()<<"\n";
      std::cerr << "File type : "<<stype<<"\n";
      int nbody; float time; bool ok;
      // get the input number of bodies according to the selection
      ok=uns->snapshot->getData("nsel",&nbody);
      // get the simulation time
      ok=uns->snapshot->getData("time",&time);

      std::cerr << "Nbody selected = " << nbody << "\nTime="<<time <<"\n";
      if (stype=="Nemo") {
        displayInfo(display,maxlines,"all",uns);
      } else {
        displayInfo(display,maxlines,"gas"  ,uns);
        displayInfo(display,maxlines,"halo" ,uns);
        displayInfo(display,maxlines,"disk" ,uns);
        displayInfo(display,maxlines,"bulge",uns);
        displayInfo(display,maxlines,"stars",uns);
        displayInfo(display,maxlines,"bndry",uns);
      }
    }
  }
  delete uns;
  //   finish NEMO
  finiparam();
}
// ------------------------------------------------------------
