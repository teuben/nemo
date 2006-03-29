// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// lagrange.cc                                                                 |
//                                                                             |
// Copyright (C) 2004-2006 Walter Dehnen                                       |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/defman.h>
#include <public/basic.h>
#include <public/tools.h>
#include <public/io.h>
#include <ctime>
#include <cmath>

namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class lagrange                                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class lagrange : public manipulator {
  private:
    int              N;
    double          *M;
    mutable double  *R;
    mutable output   OUT;
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "lagrange"; }
    const char* describe() const {
      return message("computes lagrange radii w.r.t. origin "
		     "and writes them to file %s",OUT.file());
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::m | fieldset::x; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    lagrange(const double*pars,
	     int          npar,
	     const char  *file) :
      N   ( npar ),
      M   ( npar>0? falcON_NEW(double,npar) : 0 ),
      R   ( npar>0? falcON_NEW(double,npar) : 0 ),
      OUT ( file? file : "." )
    {
      for(int i=0; i!=N; ++i) M[i] = pars[i];
      if(file == 0 || npar == 0) {
	std::cerr<<
	  " Manipulator \"lagrange\":\n"
	  " computes lagrange radii (w.r.t. origin) for given mass fractions\n"
	  " (paramaters) and writes them to data given data file\n";
      }
      if(file == 0)
	error("Manipulator \"lagrange\": no output file given\n");
      if(!OUT.is_open())
	error("Manipulator \"lagrange\": couldn't open output\n");
      if(npar == 0)
	error("Manipulator \"lagrange\": no mass fractions given\n");
      OUT  << "#\n"
	   << "# output from Manipulator \"lagrange\"\n";
      if(RunInfo::cmd_known())
	OUT<< "#\n# command: \""<<RunInfo::cmd() <<"\"\n";
      OUT  << "# run at "<<RunInfo::time()<<'\n';
      if(RunInfo::user_known())
	OUT<< "#     by \""<<RunInfo::user()<<"\"\n";
      if(RunInfo::host_known())
	OUT<< "#     on \""<<RunInfo::host()<<"\"\n";
      if(RunInfo::pid_known())
	OUT<< "#     pid "<<RunInfo::pid()<<'\n';
      OUT  << "#\n# time      ";
      for(int i=0; i!=N; ++i)
	OUT << " r[" << std::setw(4) << 100*M[i] << "%]";
      OUT  << "\n#----------";
      for(int i=0; i!=N; ++i)
	OUT << "---------";
      OUT << std::endl;
      OUT.stream().setf(std::ios::left, std::ios::adjustfield);
    }
    //--------------------------------------------------------------------------
    ~lagrange() {
      if(M) falcON_DEL_A(M); 
      if(R) falcON_DEL_A(R);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  bool lagrange::manipulate(const snapshot*S) const {
    if(S->have_not(need())) {
      warning("Manipulator \"lagrange\" insufficient data\n");
      return false;
    }
    find_lagrange_rad(S->begin_all_bodies(),N,M,R);
    OUT << std::setw(12) << S->time();
    for(int i=0; i!=N; ++i)
      OUT  << ' ' << std::setw(8) << R[i];
    OUT << std::endl;
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

__DEF__MAN(lagrange)
