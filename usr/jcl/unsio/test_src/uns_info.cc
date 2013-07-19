// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010-2013
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
  "select=all\n               select particles (range, or component name)\n"
  "                   component: gas,halo,disk,bulge,stars,bndry",
  "display=t\n                display array's content (t|f)",
  "bits=\n         (default: all , but you cat set \"mxvpXRIUMAHT\" physicals quantities that you want to display\n",
  "maxline=2\n                max lines per components",  
  "times=all\n		      selected time",
  "verbose=f\n                verbose on/off",
  "VERSION=2.1\n              compiled on <"__DATE__"> JCL  ",
  NULL,
};
const char * usage="Print information about an UNS file";
using namespace std;
void displayInfo(bool display,int maxlines, std::string comp, uns::CunsIn * uns);
template <class T> void displayFormat(int maxlines,std::string text, T * array, int dim, int size, int np);
bool component_exist=false;
// ------------------------------------------------------------
//  displayInfo
void displayInfo(bool display,int maxlines, std::string comp, uns::CunsIn * uns)
{
  float * pos, * vel, * mass, * pot , *acc, *eps;
  int * id;
  int nbody=0;
  bool ok=false;

  float * nullp;
  ok = uns->snapshot->getData(comp,"nbody" ,&nbody,&nullp);
  if (ok) component_exist=true;
  if (ok) {
    std::cout << setw(50) << setfill('=') << ""<<"\n";
    std::cout<< setfill(' ');
    std::cout<< left<< setw(8) << comp << ":" << setw(9) << right << nbody << "\n";
  }
  ok = uns->snapshot->getData(comp,"mass",&nbody,&mass);
  if (ok && display) {
    displayFormat(maxlines,"mass[1] = ",mass,1,nbody, 3);
  }  
  ok = uns->snapshot->getData(comp,"pos" ,&nbody,&pos );
  if (ok && display) {
    displayFormat(maxlines,"pos [3] = ",pos ,3,nbody, 1);
  }
  ok = uns->snapshot->getData(comp,"vel" ,&nbody,&vel );
  if (ok && display) {
    displayFormat(maxlines,"vel [3] = ",vel ,3,nbody, 1);
  }
  ok = uns->snapshot->getData(comp,"pot" ,&nbody,&pot );
  if (ok && display) {
    displayFormat(maxlines,"pot [1] = ",pot ,1,nbody, 3);
  }
  ok = uns->snapshot->getData(comp,"eps" ,&nbody,&eps );
  if (ok && display) {
    displayFormat(maxlines,"eps [1] = ",eps ,1,nbody, 3);
  }
  ok = uns->snapshot->getData(comp,"acc" ,&nbody,&acc );
  if (ok && display) {
    displayFormat(maxlines,"acc [3] = ",acc ,3,nbody, 1);
  }
  ok = uns->snapshot->getData(comp,"id"  ,&nbody,&id);
  if (ok && display) {
    displayFormat(maxlines,"id  [1] = ",id  ,1,nbody, 3);
  }  
  //if (comp == "gas") {
  float * rho, * u, * hsml, * temp, * metal;
  ok = uns->snapshot->getData(comp,"rho" ,&nbody,&rho );
  
  if (ok && display) {
    displayFormat(maxlines,"rho [1] = ",rho ,1,nbody, 3);
  }
  ok = uns->snapshot->getData(comp,"u"   ,&nbody,&u );
  if (ok && display) {
    displayFormat(maxlines,"u   [1] = ",u   ,1,nbody, 3);
  }
  ok = uns->snapshot->getData(comp,"hsml",&nbody,&hsml);
  if (ok && display) {
    displayFormat(maxlines,"hsml[1] = ",hsml,1,nbody, 3);
  }
  ok = uns->snapshot->getData(comp,"temp",&nbody,&temp);
  if (ok && display) {
    displayFormat(maxlines,"temp[1] = ",temp,1,nbody, 3);
  }
  ok = uns->snapshot->getData(comp,"metal",&nbody,&metal);
  if (ok && display) {
    displayFormat(maxlines,"metal[1] = ",metal,1,nbody, 3);
  }    
  //}
  //if (comp == "stars") {
  float * age;//, * metal;
  ok = uns->snapshot->getData(comp,"age" ,&nbody,&age );
  if (ok && display) {
    displayFormat(maxlines,"age [1] = ",age,1,nbody, 3);
  } 
  float * im;
  ok = uns->snapshot->getData(comp,"im" ,&nbody,&im );
  if (ok && display) {
    displayFormat(maxlines,"im [1] = ",im,1,nbody, 3);
  }
  float * ssl;
  ok = uns->snapshot->getData(comp,"ssl" ,&nbody,&ssl );
  if (ok && display) {
    displayFormat(maxlines,"ssl [1] = ",ssl,1,nbody, 3);
  }
  float * zs,* zsmt;
  int czs, czsmt, nzs, nzsmt;
  ok = uns->snapshot->getData(comp,"zs" ,&nzs,&zs );
  if (ok && display) {
    ok = uns->snapshot->getData("czs"   ,&czs );
    std::cerr << "nzs="<< nzs<<" czs ="<<czs<<"\n";
    displayFormat(maxlines,"zs [1] = ",zs,1,nzs, 3);
  }
  ok = uns->snapshot->getData(comp,"zsmt" ,&nzsmt,&zsmt );
  if (ok && display) {
    ok = uns->snapshot->getData("czsmt"   ,&czsmt );
    std::cerr << "nzsmt="<< nzsmt<<"  czsmt ="<<czsmt<<"\n";
    displayFormat(maxlines,"zsmt[1] = ",zsmt,1,nzsmt, 3);
  }
  float * cm; int ncm;
  ok = uns->snapshot->getData(comp,"cm" ,&ncm,&cm );
  if (ok && display) {
    displayFormat(maxlines,"cm [1] = ",cm,1,ncm, 3);
  }
//  ok = uns->snapshot->getData(comp,"metal" ,&nbody,&metal );
//  if (ok && display) {
//    displayFormat(maxlines,"metal[1] = ",metal,1,nbody, 3);
//  } 
  //}
}
// ------------------------------------------------------------
// displayFormat
template <class T>  void displayFormat(int maxlines,std::string text, T * array, int dim, int size, int np)
{
  std::cout << scientific << left << setw(11) << text;
  // First line
  for (int k=0;k<std::min(size,np);k++) {
    for (int j=0;j<dim;j++) {
      std::cout << array[k*dim+j] << " ";
    }
  }
  std::cout << "\n";
  // other lines
  for (int i=1; i<std::min(maxlines,size/min(size,np)); i+=min(size,np)) {
    std::cout << left << setw(11) << "";
    for (int k=0;k<std::min(size,np);k++) {
      for (int j=0;j<dim;j++) {
        std::cout << array[(k+(i*np))*dim+j] << " ";
      }
    }
    std::cout << "\n";
  }  
  std::cout << left << setw(11) << "" << ". . .\n";
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
  std::string bits   =(getparam ((char *) "bits"    ));
  int    maxlines    = getiparam((char *) "maxline" );
  char * select_t    = getparam ((char *) "times"   );
  bool   verbose     = getbparam((char *) "verbose" );  
  
  // -----------------------------------------------
  // instantiate a new UNS input object (for reading)
  uns::CunsIn * uns = new uns::CunsIn(simname,select_c,select_t,verbose);
  
  if (!display) bits="none"; // we don't read anything
  if (uns->isValid()) { // input file is known by UNS lib        
    while(uns->snapshot->nextFrame(bits)) { // there is a new frame
      std::string stype = uns->snapshot->getInterfaceType();
      std::string file_structure=uns->snapshot->getFileStructure();
      std::cout << setw(50) << setfill('*') << ""<<"\n";
      std::cout << "File name : "<<uns->snapshot->getFileName()<<"\n";
      std::cout << "File type : "<<stype<<"\n";
      int nbody; float time;
      // get the input number of bodies according to the selection
      uns->snapshot->getData("nsel",&nbody);
      // get the simulation time
      uns->snapshot->getData("time",&time);

      std::cout << "Nbody selected = " << nbody << "\nTime="<<time <<"\n";
      if (0 && file_structure=="range") {
        displayInfo(display,maxlines,"all",uns);
      } else {
        component_exist=false;
        displayInfo(display,maxlines,"gas"  ,uns);
        displayInfo(display,maxlines,"halo" ,uns);
        displayInfo(display,maxlines,"disk" ,uns);
        displayInfo(display,maxlines,"bulge",uns);
        displayInfo(display,maxlines,"stars",uns);
        displayInfo(display,maxlines,"bndry",uns);
        if (!component_exist) { // no comp, diplay all
          displayInfo(display,maxlines,"all",uns);
        }
      }
    }
  }
  delete uns;
  //   finish NEMO
  finiparam();
}
// ------------------------------------------------------------
