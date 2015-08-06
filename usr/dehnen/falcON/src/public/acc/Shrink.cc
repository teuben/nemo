// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/acc/Shrink.cc
///
/// \author  Walter Dehnen
///
/// \date    2009,2011
/// 
/// \brief   contains code for a shrinking/expanding potential
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009,2011  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
// 
////////////////////////////////////////////////////////////////////////////////
///
/// \version 18-Nov-2009  Created, based on Combined.cc
///
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <acceleration.h>
#include <acc/timer.h>
#define  __NO_AUX_DEFACC
#include <defacc.h>      // $NMOINC/defacc.h
#include <utils/sse.h>
////////////////////////////////////////////////////////////////////////////////
namespace {

  using namespace WDutils;

  const int AccMax = 10;

  template <typename T> struct _type;
  template <> struct _type<float > { static const char t='f'; };
  template <> struct _type<double> { static const char t='d'; };
  template <typename T> char type(T = T(0)) { return _type<T>::t; }
  //
  /// a shrinking/growing external gravitational potential.
  /*!
    Given a set of up to ten gravitational potentials \f$\Phi_i\f$, we
    implement the potential
    /f[
    \Phi(\mathbf{x},t) = \alpha(t) \sum_i \Phi_i(\alpha\mathbf{x},t)
    /f]
    with
    \f[
    \alpha(t) = 1 + (\alpha_f-1) A(t)
    \f]
    where \f$\alpha_f\f$ is the final value for $\alpha$ and \f$A(t)\f$ one
    of the growth factors in file timer.h.
    \note the maximum number of \f$\Phi_i\f$ is set to ten.
    */
  class Shrink : private timer {
    /// swallow rest of line from istream
    static void SwallowRestofLine(std::istream& from)
    { char c; do from.get(c); while(from.good() && c !='\n'); }
    /// read line until character `#'; swallow rest of line
    static void _read_line(std::istream& from, char*line, int const&n)
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
      // 2. read line until EOL
      register char*l=line+1;
      const    char*L=line+n;
      while(l != L) {
	from.get(*l);
	if(from.eof() || (*l)=='\n') { *l=0; return; }
	++l;
      }
      error("Shrink: line longer than expected\n");
    }
    /// read line of size n; skip lines whose first non-space character is #
    static void read_line(std::istream& from, char*line, int const&n)
    {
      do {
	_read_line(from,line,n);
      } while(from.good() && *line == 0);
    }
    /// \name data
    //@{
    static const int NMAX=10;              ///< max number of potentials
    int              N;                    ///< actual number of potentials
    acc_pter         AC[NMAX];             ///< pointers to accelerations
    bool             NEEDM,NEEDV;          ///< need masses, velocities?
    double           AlfaI,AlfaF;          ///< alpha_initial, alpha_final
    mutable size_t   NPOS,NPOT,NACC;       ///< # bytes allocated: POS,POT,ACC
    mutable void    *POS;                  ///< memory for scaled positions
    mutable void    *POT;                  ///< memory for external potential
    mutable void    *ACC;                  ///< memory for external acceleration
    //@}
    /// given n positions, return array with n scaled positions
    /// \param[in] n  size of arrays modulo dimensionality
    /// \param[in] x  array with positions to scale
    /// \param[in] a  scale factor
    /// \return       array with n scaled positions
    template<int NDIM, typename scalar>
    const scalar* scaled(int n, const scalar*x, scalar a) const
    {
      size_t Nscl16 = SSE::Top<scalar>(n);
      size_t Nvec16 = SSE::Top<scalar>(n*NDIM);
      if(Nvec16*sizeof(scalar) > NPOS) {
	if(POS) WDutils_DEL16(POS);
	NPOS = Nvec16*sizeof(scalar);
	POS  = WDutils_NEW16(char,NPOS);
      }
      if(Nscl16*sizeof(scalar) > NPOT) {
	if(POT) WDutils_DEL16(POT);
	NPOT = Nvec16*sizeof(scalar);
	POT  = WDutils_NEW16(char,NPOT);
      }
      if(Nvec16*sizeof(scalar) > NACC) {
	if(ACC) WDutils_DEL16(ACC);
	NACC = Nvec16*sizeof(scalar);
	ACC  = WDutils_NEW16(char,NACC);
      }
      if(a == scalar(1)) return x;
      std::memcpy(POS,x,n*NDIM*sizeof(scalar));
      for(size_t i=0; i!=Nvec16; ++i)
	reinterpret_cast<scalar*>(POS)[i] *= a;
      return reinterpret_cast<const scalar*>(POS);
    }
    /// scale factor alpha(t)
    double Alpha(double time) const
    {
      return AlfaI + (AlfaF-AlfaI) * timer::operator()(time);
    }
    /// do the job
    /// \param[in]  time  current simulation time
    /// \param[in]  nbod  # bodies
    /// \param[in]  m     array: masses
    /// \param[in]  X     array: positions
    /// \param[in]  v     array: velocities
    /// \param[in]  f     array: flags
    /// \param[out] p     array: potentials
    /// \param[out] a     array: accelerations
    /// \param[in]  add   controller for adding/assigning
    template <int NDIM, typename scalar> inline
    void acc_T(double time, int nbod,
	       const scalar*m, const scalar*X, const scalar*v, const int*f,
	       scalar*p, scalar*a, int add)
    {
      // obtain scale factor
      scalar alpha = Alpha(time);
      // if alpha = 0: issue warning! and set zero velocities
      if(alpha == 0) {
	warning("Shrink: alpha(t=%f) = 0\n",time);
	if(! (add&1))
	  for(int n=0; n!=nbod; ++n)
	    if(f==0 || f[n] & 1) 
	      p[n] = scalar(0);
	if(! (add&2))
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	    if(f==0 || f[n] & 1) 
	      v_set<NDIM>(a+nn,scalar(0));
	return;
      }
      // if alpha != 0, must compute accelerations:
      // obtain scaled positions, ensure POT, ACC have enough memory
      const scalar*x = scaled<NDIM,scalar>(nbod,X,alpha);
      // add gravity from all the acceleration fields in POT, ACC
      for(int i=0; i<N; ++i)
	(*(AC[i]))(NDIM,time,nbod,
		   static_cast<const void*>(m),
		   static_cast<const void*>(x),
		   static_cast<const void*>(v),
		   f,POT,ACC, i? 3:0, type(*m));
      // add or assign potential times alpha
      if(add & 1) {
	for(int n=0; n!=nbod; ++n)
	  if(f==0 || f[n] & 1) 
	    p[n] += alpha * reinterpret_cast<const scalar*>(POT)[n];
      } else {
	for(int n=0; n!=nbod; ++n)
	  if(f==0 || f[n] & 1) 
	    p[n]  = alpha * reinterpret_cast<const scalar*>(POT)[n];
      }
      // add or assign acceleration times alpha^2
      scalar fc = alpha*alpha;
      if(add & 2) {
	for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	  if(f==0 || f[n] & 1)
	    v_addtimes<NDIM>(a+nn, reinterpret_cast<const scalar*>(ACC)+nn,fc);
      } else {
	for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	  if(f==0 || f[n] & 1)
	    v_asstimes<NDIM>(a+nn, reinterpret_cast<const scalar*>(ACC)+nn,fc);
      }
    }
    /// copy string up to and excluding white space
    static void copy2ws(char*dest, const char*srce)
    {
      while(srce && !isspace(*srce)) *(dest++) = *(srce++);
    }
    //
  public:
    static const char* name() { return "Shrink"; }
    bool const&NeedMass() const { return NEEDM; }
    bool const&NeedVels() const { return NEEDV; }
    /// ctor
    Shrink(const double*pars,
	   int          npar,
	   const char  *file) : NEEDM(0), NEEDV(0)
    {
      if(npar < 5 || file==0)
	warning("Shrink: recognizes 3 parameters and requires a data file.\n"
		"parameters:\n"
		"  par[0] = initial value for alpha                    [???]\n"
		"  par[1] = final value for alpha                      [???]\n"
		"  par[2] = controlling growth factor, see below         [9]\n"
		"  par[3] = t0: start time for growth                    [0]\n"
		"  par[4] = tau: time scale for growth                   [1]\n"
		"  with par[2]=0: %s\n"
		"       par[2]=1: %s\n"
		"       par[2]=2: %s\n"
		"       par[2]=3: %s\n"
		"       par[2]=9: %s\n"
		"the data file must contain up to %d entries of the form\n"
		"  accname=ACCNAME\n [accpars=ACCPARS]\n [accfile=ACCFILE]\n"
		"where [] indicates optional entries. Data between a '#' and\n"
		"end-of-line are ignored (allowing comments)\n\n",
		timer::describe(timer::adiabatic),
		timer::describe(timer::saturate),
		timer::describe(timer::quasi_linear),
		timer::describe(timer::linear),
		timer::describe(timer::constant),
		NMAX);
      if(file == 0) error("Shrink: not data file given");
      // initialize parameters
      AlfaI = pars[0];
      AlfaF = pars[1];
      timer::index
	timin = (timer::index)(npar>2? int(pars[2]) : 9);
      double
	_t0    = npar>3? pars[3] : 0.,
	_tau   = npar>4? pars[4] : 1.;
      timer::init(timin,_t0,_tau);
      if(npar>5) warning("Shrink: skipped parameters beyond 5");
      nemo_dprintf (1,
		    "initializing Shrink\n"
		    " parameters : alpha_i       = %f\n"
		    "              alpha_f       = %f\n"
		    "              timer::index  = %f -> %s\n"
		    "              t0            = %f\n"
		    "              tau           = %f\n",
		    AlfaI,AlfaF,timin,timer::describe(timin),_t0,_tau);
      // now scan datafile and initialize accs
      const int size=1024;
      std::ifstream inpt(file);
      if(!inpt.good())
	error("Shrink: couldn't open file \"%s\" for input\n",file);
      char Line[size], AccName[size], AccPars[size], AccFile[size];
      char*accname=0, *accpars=0, *accfile=0;
      nemo_dprintf (1,"sub-potentials:\n");
      for(N=0; N!=NMAX;) {
	read_line(inpt,Line,size);
	if(*Line==0) break;
	// read accname
	accname=strstr(Line,"accname=");
	if(accname) {
	  copy2ws(AccName,accname+8);
	  accname=AccName;
	} else {
	  warning("Shrink: entry \"%s\" in file \"%s\" ignored\n",
		  Line,file);
	  continue;
	}
	// read accpars
	accpars=strstr(Line,"accpars=");
	if(accpars) {
	  copy2ws(AccPars,accpars+8);
	  accpars=AccPars;
	}
	// read accfile
	accfile=strstr(Line,"accfile=");
	if(accfile) {
	  copy2ws(AccFile,accfile+8);
	  accfile=AccFile;
	}
	// load acceleration
	nemo_dprintf (1," accname=%s",accname);
	if(accpars) nemo_dprintf (1," accpars=%s",accpars);
	if(accfile) nemo_dprintf (1," accfile=%s",accfile);
	nemo_dprintf (1,"\n");
	bool nm,nv;
	AC[N] = get_acceleration(accname,accpars,accfile,&nm,&nv);
	if(nm) NEEDM = 1;
	if(nv) NEEDV = 1;
	N++;
      }
      if(N==0)
	error("Shrink: couldn't read data from file \"%s\"\n",file);
      else if(N==NMAX) {
	read_line(inpt,Line,size);
	if(*Line)
	  error("Shrink: file \"%s\" contains more \"accname=...\" "
		"than anticipated\n",file);
      }
      if(NEEDM) warning("Shrink: need_masses() = true:\n");
      if(NEEDV) warning("Shrink: need_velocities() = true:\n");
    }
    /// compute accelerations
    /// \param[in]  ndim  # dimensions  
    /// \param[in]  time  current simulation time
    /// \param[in]  nbod  # positions
    /// \param[in]  m     masses
    /// \param[in]  x     positions
    /// \param[in]  v     velocities
    /// \param[in]  f     flags
    /// \param[out] p     potentials
    /// \param[out] a     accelerations
    /// \param[in]  add   indicator for adding / assigning gravity
    /// \param[in]  type  'f' or 'd' for float or double.
    void acc(int ndim, double time, int nbod,
	     const void*m, const void*x, const void*v, const int*f,
	     void*p, void*a, int add, char type)
    {
      switch(type) {
      case 'f':
	switch(ndim) {
	case 2: return acc_T<2>(time,nbod,
				static_cast<const float*>(m),
				static_cast<const float*>(x),
				static_cast<const float*>(v),
				f,
				static_cast<      float*>(p),
				static_cast<      float*>(a),
				add);
	case 3: return acc_T<3>(time,nbod,
				static_cast<const float*>(m),
				static_cast<const float*>(x),
				static_cast<const float*>(v),
				f,
				static_cast<      float*>(p),
				static_cast<      float*>(a),
				add);
	default: error("Shrink: unsupported ndim: %d",ndim);
	}
	break;
      case 'd':
	switch(ndim) {
	case 2: return acc_T<2>(time,nbod,
				static_cast<const double*>(m),
				static_cast<const double*>(x),
				static_cast<const double*>(v),
				f,
				static_cast<      double*>(p),
				static_cast<      double*>(a),
				add);
	case 3: return acc_T<3>(time,nbod,
				static_cast<const double*>(m),
				static_cast<const double*>(x),
				static_cast<const double*>(v),
				f,
				static_cast<      double*>(p),
				static_cast<      double*>(a),
				add);
	default: error("Shrink: unsupported ndim: %d",ndim);
	}
	break;
      default: error("Shrink: unknown type \"%c\"",type);
      }
    }
    /// dtor
    ~Shrink()
    {
      if(POS) WDutils_DEL16(POS); POS=0;
      if(POT) WDutils_DEL16(POT); POT=0;
      if(ACC) WDutils_DEL16(ACC); ACC=0;
    }
  } *MyAcc[AccMax] = {0};
  int AccN = 0;
  
#undef  _DEF_ACC_NO
#define _DEF_ACC_NO(NUM)					\
void acceleration##NUM(int        d,				\
		       double     t,				\
		       int        n,				\
		       const void*m,				\
		       const void*x,				\
		       const void*v,				\
		       const int *f,				\
		       void      *p,				\
		       void      *a,				\
		       int        i,				\
		       char       y)				\
{ (MyAcc[NUM])->acc(d,t,n,m,x,v,f,p,a,i,y); }
_DEF_ACC_NO(0)
_DEF_ACC_NO(1)
_DEF_ACC_NO(2)
_DEF_ACC_NO(3)
_DEF_ACC_NO(4)
_DEF_ACC_NO(5)
_DEF_ACC_NO(6)
_DEF_ACC_NO(7)
_DEF_ACC_NO(8)
_DEF_ACC_NO(9)
  acc_pter Accs[AccMax] = {&acceleration0,
			   &acceleration1,
			   &acceleration2,
			   &acceleration3,
			   &acceleration4,
			   &acceleration5,
			   &acceleration6,
			   &acceleration7,
			   &acceleration8,
			   &acceleration9};
} // namespace {
////////////////////////////////////////////////////////////////////////////////
void iniacceleration(const double*pars,      // I:  array with parameters       
		     int          npar,      // I:  number of parameters        
		     const char  *file,      // I:  data file name              
		     acc_pter    *accel,     // O:  pter to acceleration()      
		     bool        *needM,     // O:  acceleration() needs masses?
		     bool        *needV)     // O:  acceleration() needs vel's? 
{
  if(AccN == AccMax) {
    warning("iniacceleration(): request to initialize "
	    "more than %d accelerations of type \"Shrink\"", AccMax);
    *accel = 0;
    return;
  }
  MyAcc[AccN] = new Shrink(pars,npar,file);
  if(needM) *needM = (MyAcc[AccN])->NeedMass();
  if(needV) *needV = (MyAcc[AccN])->NeedVels();
  *accel = Accs[AccN++];
}
////////////////////////////////////////////////////////////////////////////////
