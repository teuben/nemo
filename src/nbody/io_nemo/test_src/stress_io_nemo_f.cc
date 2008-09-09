// =============================================================================
// Copyright Jean-Charles LAMBERT - 2005                                        
// e-mail:   Jean-Charles.Lambert@oamp.fr                                       
// address:  Dynamique des galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                           
//           2, place Le Verrier                                                
//           13248 Marseille Cedex 4, France                                    
//           CNRS U.M.R 6110                                                    
// =============================================================================
// stress_io_nemo_f.cc                                                          
//                                                                              
// program to stress io_nemo_f_ function, usually used from Fortram program     
// -----------------------------------------------------------------------------
#include <iostream>                                   // C++ I/O     
#include <fstream>                                    // C++ file I/O
#include <sstream>
#include <cstdio>                    // changed from stdio.h  WD, Sep 2008
#include <cstdlib>                   // changed from stdlib.h WD, Sep 2008
#include <assert.h>
#include <vector>
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
extern "C" {
#include <nemo.h>                                     // NEMO basics
  int io_nemo_f_(const char *, const int&, const int&, const char *,...);
  int close_io_nemo_f_(const char *, const int&);
}

using namespace std; // prevent writing statment like 'std::cerr'

class t_select {
 public:
  t_select() {};
  t_select(const t_select& m) {
    inemo_name = m.inemo_name;
    onemo_name = m.onemo_name;
    sel_part   = m.sel_part;
    nbody      = m.nbody;
  };
  const t_select& operator=(const t_select& m) {
    inemo_name = m.inemo_name;
    onemo_name = m.onemo_name;
    sel_part   = m.sel_part;
    nbody      = m.nbody;
    return *this;
  };

 public:
  std::string inemo_name;
  std::string onemo_name;
  std::string sel_part;
  int         nbody;
};
#define MAXFILES 100

class VClass {
 public:
  VClass() {};
  virtual ~VClass() {};
  virtual void proceed() {};
};
template <class T> class StressIO_Nemo: public VClass {
 public:
  //
  StressIO_Nemo(::string,std::string);
  ~StressIO_Nemo();
  void proceed();
  void proceedRW();
  void printListFile();
 private:
  int readListFile();
  bool read();
  bool write();
  bool readWrite();
 public:
  T * pos[MAXFILES], *vel[MAXFILES], *mass[MAXFILES], *tps[MAXFILES],
    * acc[MAXFILES], * pot[MAXFILES];
  int   *nbody[MAXFILES],*bits[MAXFILES];

    
 private:
  std::string in,ftype;
  int * in_status, * out_status;
  int nfiles;
  std::ostringstream inemo_name[MAXFILES];
  std::ostringstream onemo_name[MAXFILES];
  vector<t_select> sim_param;

};
//------------------------------------------------------------------------------
template <class T> StressIO_Nemo<T>::StressIO_Nemo(::string _in,std::string _ftype)
{
  in=_in;

  nfiles=readListFile();     // read input file
  printListFile();    // print out input parameters
  in_status  = new int[nfiles];
  out_status = new int[nfiles];

  // init
  for (int i=0; i< nfiles; i++) {
    //inemo_name[i] <<  in  << "_" << i;
    //onemo_name[i] <<  out << "_" << i;
    //std::cerr << "file: " << inemo_name[i].str() << "\n";
    pos  [i]   = new T[sim_param[i].nbody*3];
    vel  [i]   = new T[sim_param[i].nbody*3];
    mass [i]   = new T[sim_param[i].nbody];
    pot  [i]   = new T[sim_param[i].nbody];
    acc  [i]   = new T[sim_param[i].nbody*3];
    tps  [i]   = new T;
    nbody[i]   = new int;
    bits [i]   = new int[sim_param[i].nbody*3];
  }
  ftype = _ftype;
}
//------------------------------------------------------------------------------
template <class T> StressIO_Nemo<T>::~StressIO_Nemo()
{
  if (nfiles) {
    delete [] in_status;
    delete [] out_status;
  }
  for (int i=0; i< nfiles; i++) {
    if ( pos[i]) {
      delete [] pos[i];
    }
    if ( vel[i]) {
      delete [] vel[i];
    }
    if ( mass[i]) {
      delete [] mass[i];
    }
    if ( pot[i]) {
      delete [] pot[i];
    }
    if ( acc[i]) {
      delete [] acc[i];
    }
    if ( tps[i]) {
      delete tps[i];
    }
    if ( nbody[i]) {
      delete nbody[i];
    }
    if ( bits[i]) {
      delete [] bits[i];
    }
  }
}
//------------------------------------------------------------------------------
template <class T> int StressIO_Nemo<T>::readListFile()
{
  t_select tx;
  std::ifstream fi(in.c_str());
  
  if (! fi.is_open()) {
    std::cerr << "failed to open ["<<in<<"]\n";
    std::exit(1);
  } 
  int i=0;
  while (! fi.eof() && i <=MAXFILES) {
    i++;
    std::string x;
    fi >> x;            // read input filename
    if ( ! fi.eof()) {
      tx.inemo_name = x;
      fi >> tx.nbody;   // read nbody
      fi >> x;          // read  selected particle range
      tx.sel_part = x+"#";
      tx.onemo_name = tx.inemo_name + "_out";
      sim_param.push_back(tx); // add simulation parameters
    }
  }
  return sim_param.size(); 
}
//------------------------------------------------------------------------------
template <class T> void StressIO_Nemo<T>::printListFile()
{
  for (unsigned int i=0; i<sim_param.size(); i++) {
    std::cerr << "[" << i << "]\n";
    std::cerr << sim_param[i].inemo_name << "\n" 
	      << sim_param[i].onemo_name << "\n" 
	      << sim_param[i].sel_part  << "\n";
  }
}
//------------------------------------------------------------------------------
template <class T> void StressIO_Nemo<T>::proceed()
{
  bool run=true;
  while (run) {
    run=read();
    if (run) {
      write();
    }
  }
}
//------------------------------------------------------------------------------
template <class T> void StressIO_Nemo<T>::proceedRW()
{
  bool run=true;
  while (run) {
    run=readWrite();
  }
}
//------------------------------------------------------------------------------
template <class T> bool StressIO_Nemo<T>::readWrite()
{
  bool status=false;
  std::string read_param=ftype+",read,n,t,x,v,m,acc,pot,sp,info";
  std::string save_param=ftype+",save,n,t,x,v,m,acc,pot,sp,info";

  for (int i=0; i<nfiles; i++) {
    in_status[i] = io_nemo_f_(sim_param[i].inemo_name.c_str(), sim_param[i].inemo_name.length(),sim_param[i].nbody,
			   read_param.c_str(),
			   nbody[i],tps[i],pos[i],vel[i],mass[i],
			   acc[i],pot[i],sim_param[i].sel_part.c_str());
    if (in_status[i]) {
      std::cerr << "READ ["<<i<<"] file <" << sim_param[i].inemo_name << "> nbody="<< *nbody[i] << "\n";
      std::cerr << "Pos address =" << pos[i] << "\n";
      out_status[i] = io_nemo_f_(sim_param[i].onemo_name.c_str(),sim_param[i].inemo_name.length(),sim_param[i].nbody,
				 save_param.c_str(),
				 nbody[i],tps[i],pos[i],vel[i],mass[i],acc[i],pot[i]);
      if (out_status[i]) {
	std::cerr << "WRITE ["<<i<<"] file <" << sim_param[i].onemo_name << "> nbody="<< *nbody[i] << "\n";
	std::cerr << "Pos address =" << pos[i] << "\n";
	status=true;
      }
      else {
	close_io_nemo_f_(sim_param[i].onemo_name.c_str(),sim_param[i].inemo_name.length());
      }
    }
    else {
      close_io_nemo_f_(sim_param[i].inemo_name.c_str(), sim_param[i].inemo_name.length());
    }
  }
  return status;
}

//------------------------------------------------------------------------------
template <class T> bool StressIO_Nemo<T>::read()
{
  bool status=false;
  std::string read_param=ftype+",read,n,t,x,v,m,acc,pot,sp,info";

  for (int i=0; i<nfiles; i++) {
    in_status[i] = io_nemo_f_(sim_param[i].inemo_name.c_str(),sim_param[i].inemo_name.length(),
			      sim_param[i].nbody,
			      read_param.c_str(),
			      nbody[i],tps[i],pos[i],vel[i],mass[i],acc[i],pot[i],sim_param[i].sel_part.c_str());
    if (in_status[i]) {
      std::cerr << "READ ["<<i<<"] file <" << sim_param[i].inemo_name << "> nbody="<< *nbody[i] << "\n";
      std::cerr << "Pos address =" << pos[i] << "\n";
      status=true;
    }
    else {
      close_io_nemo_f_(sim_param[i].inemo_name.c_str(), sim_param[i].inemo_name.length());
    }
  }
  return status;
}
//------------------------------------------------------------------------------
template <class T> bool StressIO_Nemo<T>::write()
{
  bool status=false;
  std::string save_param=ftype+",save,n,t,x,v,m,acc,pot,sp,info";

  for (int i=0; i<nfiles; i++) {
    if (in_status[i]) {   // save only if reading successed
      out_status[i] = io_nemo_f_(sim_param[i].onemo_name.c_str(),sim_param[i].inemo_name.length(),sim_param[i].nbody,
				 save_param.c_str(),
				 nbody[i],tps[i],pos[i],vel[i],mass[i],acc[i],pot[i]);
      if (out_status[i]) {
	std::cerr << "WRITE ["<<i<<"] file <" << sim_param[i].onemo_name << "> nbody="<< *nbody[i] << "\n";
	std::cerr << "Pos address =" << pos[i] << "\n";
	status=true;
      }
      else {
	close_io_nemo_f_(sim_param[i].onemo_name.c_str(),sim_param[i].inemo_name.length());
      }
    }
  }
  return status;
}

//------------------------------------------------------------------------------
//                             M   A   I   N                                    
//------------------------------------------------------------------------------
// NEMO parameters
::string defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n             txt input file          ",
  "precision=double\n       double (default) | float",  
  "VERSION=1.0\n       compiled on <"__DATE__"> JCL  ",
  NULL
};
::string usage="stress test on io_nemo_f function";

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  ::string in;
  ::string ftype;

  //   start  NEMO
  initparam(argv,defv);

  // Get parameters
  in     = getparam("in");
  ftype  = getparam("precision");

  std::string sftype=ftype;
  VClass * ltest;
  if (sftype == "double" ) {
    ltest = new  StressIO_Nemo<double>(in,sftype);
  }
  else {
    ltest = new  StressIO_Nemo<float>(in,sftype);
  }
  ltest->proceed();
  delete ltest;
  //   finish NEMO
  finiparam();
}
// ----------- End Of [stress_io_nemo_f.cc] ------------------------------------
