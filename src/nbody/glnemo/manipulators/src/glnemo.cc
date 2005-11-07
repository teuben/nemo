// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// glnemo.cc                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Jean-Charles LAMBERT - 2005                                       |
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      |
// address:  Dynamique des galaxies                                            |
//           Laboratoire d'Astrophysique de Marseille                          |
//           2, place Le Verrier                                               |
//           13248 Marseille Cedex 4, France                                   |
//           CNRS U.M.R 6110                                                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// glnemo manipulator transform gyrfalcON simulation program in a glnemo server|
// allowing real time 3D rendering of the running gyrfalcON simulation via     |
// the glnemo program. (see $NEMOSRC/nbody/glnemo).                            |
//                                                                             |
// The manipulator use socket to communicate with glnemo program. The network  |
// server part of the manipulator is multithreaded, that means that gyrfalcON  |
// simulation keep running (computing) while it is communicating with glnemo   |
// client program. The communication between server (manipulator) and client   |
// (glnemo) is "safe" that means that you can interrupt whenever you want      |
// glnemo program while it is "talking" to gyrfalcON program. All the nework   |
// ressources will be deallocated properly.                                    |
//                                                                             |
//-----------------------------------------------------------------------------+

#include <public/defman.h>
#include <nemo.h>
#include "master_server_thread.h"

namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class glnemo                                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class glnemo : public manipulator {
  private:
    const   int      port;   // communication port  
    mutable  MasterServerThread * glnemo_server;
    mutable  int current_step;
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "glnemo"; }
    const char* describe() const {
      return message("Start glnemo server");
    }
    //--------------------------------------------------------------------------
    fieldset          need    () const { return fieldset::o; }
    fieldset          change  () const { return fieldset::vectors; }
    fieldset          provide () const { return fieldset::o; }
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    // constructor:                                                             
    glnemo(const double *pars,
		  int          npar,
		  const char   *file) :
      port(npar>0?int (pars[0]):4444)    // default port = 4444
    {
      if (nemo_debug(1) || nemo_debug(2)) {
	std::cerr<<
	  " Manipulator \"glnemo\":\n"
	  " Start glnemo server on dedicated port"
	  " par[0] = port, int  [ default 4444 - you must give a value]\n";
      }
      current_step=0;

    }
    //--------------------------------------------------------------------------
    ~glnemo() {
      std::cerr << "In glnemo destructor...\n";
      delete glnemo_server;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  bool glnemo::manipulate(const snapshot*S) const {
    if (current_step==0) {  // this is firs time
      std::string sim_name;
      if (hasvalue("out")) sim_name = getparam("out");
      else                 sim_name = "noname";
      // instantiate Master GLnemo Server Thread
      glnemo_server = new MasterServerThread(sim_name,port,S);
    }
    else {
      glnemo_server->updateData(S);
    }
    current_step++;
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

#ifndef ALREADY_DEF_MAN
__DEF__MAN(glnemo)
#endif  
