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
#include <vector>

#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>                                     // NEMO basics
#include <io_nemo.h>                                     // NEMO basics


using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
// NEMO parameters
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           INGO input file          ",
  "out=???\n          Nemo output file                        ",
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL  ",
  NULL
};
const char * usage="Save in NEMO format a glnemo indexes list input file";

std::vector <float> vr;    // vector to store positions 
std::vector <float> vv;    // vector to store velocities
std::vector <float> vm;    // vector to store masses    
std::vector <int>   vi;    // vector to store indexes    

float timex;

//------------------------------------------------------------------------------
// readIndexList
void readIndexList(std::string listname)
{
  std::ifstream         // File Handler
    fd;               // manipfile file desc
  
  // open file of velocities table
  fd.open(listname.c_str(),std::ios::in);
  if ( ! fd.is_open()) {
    std::cerr <<
      "Unable to open ["<<listname<<"] for input, aborting..\n\n";
    std::exit(1);
  }
  std::string line;
  // Read Header
  getline(fd,line); // 1st line
  getline(fd,line); // 2nd line
  getline(fd,line); // 3rd line
  std::istringstream sst(line);
  sst >> timex;

  // Read 
  while (! fd.eof()) {           // while ! eof
    std::string line;
    getline(fd,line);
    if ( ! fd.eof()) {
      int index;
      float value;
      std::istringstream ss(line);
      ss >> index; // read index
      vi.push_back(index);
      ss >> value;
      vm.push_back(value);
      for (int i=0;i<3;i++){
	ss >> value;
	vr.push_back(value);
      }
      for (int i=0;i<3;i++){
	ss >> value;
	vv.push_back(value);
      }
    }
  }
}

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));

  // Get parameters
  std::string listname (getparam((char *) "in"    ));
  std::string outnemo (getparam((char *) "out"   ));

  readIndexList(listname);
  int nbody = vm.size();
  int * n   = &nbody;
  float * t = &timex;

  float * pos = &(vr[0]);// new float[3*vi.size()];
  float * vel = &(vv[0]);//new float[3*vi.size()];
  float * mass= &(vm[0]); //new float[  vi.size()];

  io_nemo((char *) (outnemo.c_str()),(char *)"info,float,save,n,t,x,v,m",&n,&t,
	  &pos,&vel,&mass);

  io_nemo(const_cast<char *>(outnemo.c_str()),(char *)"close");

  //   finish NEMO
  finiparam();
}
// ----------- End Of [cell2nemo.cc] --------------------------------------------
