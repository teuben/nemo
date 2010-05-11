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
#include "uns.h"
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>                                     // NEMO basics
#include <io_nemo.h>                                     // NEMO basics


using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
// NEMO parameters
const char * defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n           UNS input file          ",
  "out=???\n          Nemo output file                        ",
  "index=???\n        glnemo2 indexes input file     ",
  "select=???\n       component selected (disk,stars,halo,gas,range)",
  "time=all\n         selected time",
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL  ",
  NULL
};
const char * usage="Save in NEMO format a glnemo indexes list input file";

std::vector <float> vr;    // vector to store positions 
std::vector <float> vv;    // vector to store velocities
std::vector <float> vm;    // vector to store masses    
std::vector <int>   vi;    // vector to store indexes    

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
  getline(fd,line);
  if (line != "#glnemo_index_list") {
    std::cerr <<"Input file ["<<listname<<" is not a know glnemo"
	      <<"index list file....aborting\n";
    std::exit(1);
  }
  // Read 
  while (! fd.eof()) {           // while ! eof
    std::string line;
    getline(fd,line);
    if ( ! fd.eof()) {
      int index;
      std::istringstream ss(line);
      ss >> index; // read index
      vi.push_back(index);
    }
  }
}

//------------------------------------------------------------------------------
// nemoOut
void nemoOut(uns::CunsIn * uns, int nbody,std::string outnemo,float * pos, float *vel, 
	     float * mass)
{

  int nsel;
  float timex, *x,*v,*m;
  uns->snapshot->getData("time",&timex);
  uns->snapshot->getData("pos" ,&nsel,&x);
  uns->snapshot->getData("vel" ,&nsel,&v);
  uns->snapshot->getData("mass",&nsel,&m);
  
  float * t = &timex;

  std::cerr << "nbody=" << nbody << " time="<<timex <<"\n";
  int ii=0;
  for (std::vector<int>::iterator i=vi.begin(); i<vi.end(); i++) {
    // mass
    mass[ii] = m[(*i)];
    // pos
    pos[ii*3+0] = x[(*i)*3+0];
    pos[ii*3+1] = x[(*i)*3+1];
    pos[ii*3+2] = x[(*i)*3+2];
    // vel
    vel[ii*3+0] = v[(*i)*3+0];
    vel[ii*3+1] = v[(*i)*3+1];
    vel[ii*3+2] = v[(*i)*3+2];

    ii++;
  }
  int    nn = vi.size(); // number of indexes
  assert(ii==nn);
  int   * n = &nn;
  io_nemo((char *) (outnemo.c_str()),(char *)"info,float,save,n,t,x,v,m",&n,&t,
	  &pos,&vel,&mass);

}
//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  //   start  NEMO
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  if (argc) {;} // remove compiler warning :)
  // Get parameters
  std::string simname (getparam((char *) "in"    ));
  std::string outnemo (getparam((char *) "out"   ));
  std::string listname(getparam((char *) "index" ));
  char * select_c  = getparam((char *) "select");
  char * select_t  = getparam((char *) "time"  );

  readIndexList(listname);
  float * pos = new float[3*vi.size()];
  float * vel = new float[3*vi.size()];
  float * mass= new float[  vi.size()];

  int ok=1;
  // instantiate a new uns object
  //s::Cuns * uns = new uns::Cuns(simname,select_c,select_t);
  uns::CunsIn * uns = new uns::CunsIn(simname.c_str(),select_c,select_t);
  if (uns->isValid()) {
    while(uns->snapshot->isNewFrame()&&ok>0) {
      uns::UserSelection user_select;
      // get CRV
      uns::ComponentRangeVector * crv=uns->snapshot->getSnapshotRange();
      // user select request
      user_select.setSelection(uns->snapshot->getSelectPart(),crv);
      // load next frame
      ok=uns->snapshot->nextFrame(user_select.getIndexesTab(),user_select.getNSel());
      if (ok>0) {
	int   nbody = user_select.getNSel();
	nemoOut(uns,nbody,outnemo,pos,vel,mass);
      }
    }
    io_nemo(const_cast<char *>(outnemo.c_str()),(char *)"close");
  } else {
    std::cerr << "Unknown UNS file format["<<simname<<"]\n";
  }

  //   finish NEMO
  finiparam();
}
// ----------- End Of [cell2nemo.cc] --------------------------------------------
