// -*- C++ -*-																												|
//-----------------------------------------------------------------------------------------------------------------------------+
//																																|
// glnemo2.cc																											|
//																																|
// C++ code																												|
//																																|
// Copyright Jean-Charles LAMBERT - 2005-2010															|
// e-mail:   Jean-Charles.Lambert@oamp.fr																	|
// address:  Dynamique des galaxies																			|
//           Laboratoire d'Astrophysique de Marseille															|
//           2, place Le Verrier																							|
//           13248 Marseille Cedex 4, France																	|
//           CNRS U.M.R 7326																						|
//																																|
//-----------------------------------------------------------------------------------------------------------------------------+
//																																|
// glnemo2 manipulator transform gyrfalcON simulation program in a glnemo2 server	|
// allowing real time 3D rendering of the running gyrfalcON simulation via						|
// the glnemo2 program. (see $NEMOSRC/nbody/glnemo2).											|
//																																|
// The manipulator use socket to communicate with glnemo2 program. The network		|
// server part of the manipulator is multithreaded, that means that gyrfalcON				|
// simulation keep running (computing) while it is communicating with glnemo2				|
// client program. The communication between server (manipulator) and client				|
// (glnemo2) is "safe" that means that you can interrupt whenever you want				|
// glnemo2 program while it is "talking" to gyrfalcON program. All the nework				|
// ressources will be deallocated properly.																	|
//																																|
//-----------------------------------------------------------------------------------------------------------------------------+

#include <public/defman.h>
#include <nemo.h>
#include "master_server_thread.h"

namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class glnemo                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class glnemo2 : public manipulator {
  private:
	const   int      port;     // communication port
	const   int      max_port; // max communication port
	mutable  MasterServerThread * glnemo2_server;
	mutable  int current_step;
	//--------------------------------------------------------------------------
  public:
	const char* name    () const { return "glnemo2"; }
	const char* describe() const {
	  return message("Start glnemo2 server");
	}
	//--------------------------------------------------------------------------
	fieldset          need    () const { return fieldset::empty; }
	fieldset          change  () const { return fieldset::vectors; }
	fieldset          provide () const { return fieldset::empty; }
	bool manipulate(const snapshot*) const;
	//--------------------------------------------------------------------------
	// constructor:
	glnemo2(const double *pars,
		  int          npar,
		  const char   *file) :
	  port(npar>0?int (pars[0]):4000),    // default port = 4444
      max_port(npar>1?int (pars[1]):2)    // 8 ports
	{
          if (file) {;} // get rid of compiler warning
	  if (nemo_debug(1) || nemo_debug(2)) {
	std::cerr<<
	  " Manipulator \"glnemo2\":\n"
	  " Start glnemo2 server on dedicated port"
	  " par[0] = port, int  [ default 4444 - you must give a value]\n";
	  }
	  current_step=0;
	}
	//--------------------------------------------------------------------------
	~glnemo2() {
	  delete glnemo2_server;
	}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool glnemo2::manipulate(const snapshot*S) const {
	if (current_step==0) {  // this is firs time
	  std::string sim_name;
	  if (hasvalue((char *)"out")) sim_name = getparam((char *) "out");
	  else                 sim_name = "noname";
	  // instantiate Master GLnemo2 Server Thread

	  glnemo2_server = new MasterServerThread(sim_name,port,max_port,S);
	}
	else {
	  glnemo2_server->updateData();
	}
	current_step++;
	return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

#ifndef ALREADY_DEF_MAN
__DEF__MAN(glnemo2)
#endif
