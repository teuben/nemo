//-----------------------------------------------------------------------------+
//                                                                             |
// GrowAccs.cc                                                                 |
//                                                                             |
// Copyright (C) 2004 Walter Dehnen                                            |
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
//                                                                             |
// Versions                                                                    |
// 0.0    17/08/2004  WD created                                               |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
#include <acceleration.h>
#include <timer.h>
#define __NO_AUX_DEFACC
#include <defacc.h>
//=============================================================================#
// define C++ implementation of growing accelerations                          |
//=============================================================================#
namespace {
  //----------------------------------------------------------------------------
  inline void SwallowRestofLine(std::istream& from)
    // swallow rest of line from istream
  {
    char c;
    do from.get(c); while(from.good() && c !='\n');
  }
  //----------------------------------------------------------------------------
  inline void __read_line(std::istream& from, char*line, int const&n)
    // read line until character `#'; swallow rest of line
  {
    // 1. find first non-space character whereby:
    //    - return empty line if EOF, EOL, '#' are encountered
    do {
      from.get(*line);
      if(from.eof() || (*line)=='\n') {
	*line=0;
	return;
      }
      if((*line)=='#') {
	SwallowRestofLine(from);
	*line=0;
	return;
      }
    } while(isspace(*line));
    // 2. read line until first white-space character
    //    swallow rest til EOL, if necessary 
    register char*l=line+1;
    const    char*L=line+n;
    while(l != L) {
      from.get(*l);
      if(from.eof() || (*l)=='\n') { *l=0; return; }
      if(isspace(*l)) {
	*l=0;
	SwallowRestofLine(from);
	return;
      }
      ++l;
    }
    error("GrowAccs: line longer than expected\n");
  }
  //----------------------------------------------------------------------------
  inline void read_line(std::istream& from, char*line, int const&n)
    // read line of size n
    // skip lines whose first non-space character is #
  {
    do {
      __read_line(from,line,n);
    } while(from.good() && *line == 0);
  }
  //////////////////////////////////////////////////////////////////////////////
  timer     TIMER;                          // our timer
  const int NAMAX=10;                       // max number of potentials
  int       NA;                             // actual number of potentials
  acc_pter  AC[NAMAX];                      // pointers to accelerations
  //----------------------------------------------------------------------------
  void iniacceleration(const double*pars,   // I:  array with parameters
		       int          npar,   // I:  number of parameters
		       const char  *file,   // I:  data file name
		       bool        *needM,  // O:  acceleration() needs masses?
		       bool        *needV)  // O:  acceleration() needs vel's?
  {
    if(npar < 3 || file==0)
      warning("GrowAccs acceleration field recognizes 3 parameters\n"
	      "and requires a datafile.\n"
	      "parameters:\n"
	      "  par[0] = controlling growth factor, see below         [9]\n"
	      "  par[1] = t0: start time for growth                    [0]\n"
	      "  par[2] = tau: time scale for growth                   [1]\n"
	      "  with par[0]=0: %s\n"
	      "       par[0]=1: %s\n"
	      "       par[0]=2: %s\n"
	      "       par[0]=3: %s\n"
	      "       par[0]=9: %s\n"
	      "the data file must contain up to %d entries of the form\n"
	      "  accname=ACCNAME\n[accpars=ACCPARS]\n[accfile=ACCFILE]\n"
	      "where [] indicates optional entries. Data between a '#' and\n"
	      "end-of-line are ignored (allowing comments)\n",
	      timer::describe(timer::adiabatic),
	      timer::describe(timer::saturate),
	      timer::describe(timer::quasi_linear),
	      timer::describe(timer::linear),
	      timer::describe(timer::constant),
	      NAMAX);
    if(file == 0) error("GrowAccs: no data file given");
    // initialize timer
    timer::index
      timin = (timer::index)(npar>0? int(pars[0]) : 9);
    double
      t0    = npar>1? pars[1] : 0.,
      tau   = npar>2? pars[2] : 1.;
    TIMER.init(timin,t0,tau);
    // now scan datafile and initialize accs
    const int size=200;
    std::ifstream in(file);
    if(!in.good())
      error("GrowAccs: couldn't open file \"%s\" for input\n",file);
    char Line[size], AccName[size], AccPars[size], AccFile[size];
    char*accname=0, *accpars=0, *accfile=0, unknown[9];
    read_line(in,Line,size);
    if(*Line == 0)
      error("GrowAccs: couldn't read data from file \"%s\"\n",file);
    if(strncmp(Line,"accname=",8))
      error("GrowAccs: first entry in file \"%s\" isn't \"accname=...\"\n",
	    file);
    NA     = 0;
    *needM = 0;
    *needV = 0;
    do {
      if(accname==0) {
	strcpy(AccName,Line+8);
	accname=AccName;
	accpars=0;
	accfile=0;
      }
      read_line(in,Line,size);
      if        (*Line == 0 || 0==strncmp(Line,"accname=",8)) {
	if(NA == NAMAX)
	  error("file \"%s\": more accnames than anticipated\n",file);
	//       // TEST
	//       std::cerr<<"loading acceleration "<<NA<<": \n"
	// 	       <<"  accname="<<accname<<'\n';
	//       if(accpars) std::cerr<<"  accpars="<<accpars<<'\n';
	//       if(accfile) std::cerr<<"  accfile="<<accfile<<'\n';
	//       // tensor_set
	bool nm,nv;
	AC[NA] = get_acceleration(accname,accpars,accfile,&nm,&nv);
	*needM = *needM || nm;
	*needV = *needV || nv;
	NA++;
	accname=0;
      } else if(0==strncmp(Line,"accpars=",8)) {
	if(accpars)
	  warning("GrowAccs: additional \"accpars=\""
		  " in file \"%s\" ignored\n",file);
	else {
	  strcpy(AccPars,Line+8);
	  accpars = AccPars;
	}
      } else if(0==strncmp(Line,"accfile=",8)) {
	if(accfile)
	  warning("GrowAccs: additional \"accfile=\""
		  " in file \"%s\" ignored\n",file);
	else {
	  strcpy(AccFile,Line+8);
	  accfile = AccFile;
	  if(0==strcmp(file,accfile))
	    error("GrowAccs: recursion accfile=datafile not allowed\n");
	}
      } else {
	strncpy(unknown,Line,8); unknown[9]=0;
	warning("GrowAccs: entry \"%s\" in file \"%s\" ignored\n",
		unknown,file);
      }
    } while(*Line != 0);
  }
  //----------------------------------------------------------------------------
  template <typename scalar> struct __type;
  template <> struct __type<float > { static const char t='f'; };
  template <> struct __type<double> { static const char t='d'; };
  //----------------------------------------------------------------------------
  template <int NDIM, typename scalar> inline
  void acc_T(double       time,              // I: simulation time              
	     int          nbod,              // I: number bodies =size of arrays
	     const scalar*mass,              // I: masses:         m[i]         
	     const scalar*pos,               // I: positions       (x,y,z)[i]   
	     const scalar*vel,               // I: velocities      (u,v,w)[i]   
	     const int   *flag,              // I: flags           f[i]         
	     scalar      *pot,               // O: potentials      p[i]         
	     scalar      *acc,               // O: accelerations   (ax,ay,az)[i]
	     int          add)               // I: add or assign pot & acc?     
  {
    // obtain amplitude of growth factor
    scalar ampl = TIMER(time);

    // if amplitude == 0, set pot / acc to zero if assigning required.
    if(ampl == 0.) {
      if(! (add&1))
	for(int n=0; n!=nbod; ++n)
	  if(flag==0 || flag[n] & 1) 
	    pot[n] = scalar(0);
      if(! (add&2))
	for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	  if(flag==0 || flag[n] & 1) 
	    v_set<NDIM>(acc+nn,scalar(0));
      return;
    }
    
    // if amplitude != 0, must compute accelerations:
    // if amplitude != 1, create arrays to write pot & acc into
    scalar*_pot = ampl!=1 ? new scalar[nbod]      : 0;
    scalar*_acc = ampl!=1 ? new scalar[nbod*NDIM] : 0;

    // define references to the arrays actually passed to accelerations
    scalar*&pots = ampl!=1 ? _pot : pot;
    scalar*&accs = ampl!=1 ? _acc : acc;

    // add/assign gravity from all the acceleration fields
    for(int iA=0; iA<NA; ++iA)
      (*(AC[iA]))(NDIM,time,nbod,
		  static_cast<const void*>(mass),
		  static_cast<const void*>(pos),
		  static_cast<const void*>(vel),
		  flag,
		  static_cast<      void*>(pots),
		  static_cast<      void*>(accs),
		  iA   != 0? 3 :            // further? add pot & acc           
		  ampl != 1? 0 : add,       // first? need to mul? ass:input    
		  __type<scalar>::t);

    // if amplitude != 1, multiply gravity by amplitude
    if(ampl!=1) {

      // add or assign potential times amplitude
      if(add & 1) {
	for(int n=0; n!=nbod; ++n)
	  if(flag==0 || flag[n] & 1) 
	    pot[n] += ampl * pots[n];
      } else {
	for(int n=0; n!=nbod; ++n)
	  if(flag==0 || flag[n] & 1) 
	    pot[n]  = ampl * pots[n];
      }
      delete[] _pot;

      // add or assign acceleration times amplitude
      if(add & 2) {
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	    if(flag==0 || flag[n] & 1)
	      v_addtimes<NDIM>(acc+nn, accs+nn, ampl);
      } else {
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	    if(flag==0 || flag[n] & 1)
	      v_asstimes<NDIM>(acc+nn, accs+nn, ampl);
      }
      delete[] _acc;
    }
  };
  //----------------------------------------------------------------------------
  void acceleration(int        ndim,         // I: number of dimensions         
		    double     time,         // I: simulation time              
		    int        nbod,         // I: number bodies =size of arrays
		    const void*mass,         // I: masses:         m[i]         
		    const void*pos,          // I: positions       (x,y,z)[i]   
		    const void*vel,          // I: velocities      (u,v,w)[i]   
		    const int *flag,         // I: flags           f[i]         
		    void      *pot,          // O: potentials      p[i]         
		    void      *acc,          // O: accelerations   (ax,ay,az)[i]
		    int        add,          // I: indicator (see note 6 above) 
		    char       type)         // I: type: 'f' or 'd'             
  {
    switch(type) {
    case 'f':
      switch(ndim) {
      case 2: return acc_T<2>(time,nbod,
			      static_cast<const float*>(mass),
			      static_cast<const float*>(pos),
			      static_cast<const float*>(vel),
			      flag,
			      static_cast<      float*>(pot),
			      static_cast<      float*>(acc),
			      add);
      case 3: return acc_T<3>(time,nbod,
			      static_cast<const float*>(mass),
			      static_cast<const float*>(pos),
			      static_cast<const float*>(vel),
			      flag,
			      static_cast<      float*>(pot),
			      static_cast<      float*>(acc),
			      add);
      default: error("GrowAccs: unsupported ndim: %d",ndim);
      }
      break;
    case 'd':
      switch(ndim) {
      case 2: return acc_T<2>(time,nbod,
			      static_cast<const double*>(mass),
			      static_cast<const double*>(pos),
			      static_cast<const double*>(vel),
			      flag,
			      static_cast<      double*>(pot),
			      static_cast<      double*>(acc),
			      add);
      case 3: return acc_T<3>(time,nbod,
			      static_cast<const double*>(mass),
			      static_cast<const double*>(pos),
			      static_cast<const double*>(vel),
			      flag,
			      static_cast<      double*>(pot),
			      static_cast<      double*>(acc),
			      add);
      default: error("GrowAccs: unsupported ndim: %d",ndim);
      }
      break;
    default: error("GrowAccs: unknown type \"%c\"",type);
    }
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
