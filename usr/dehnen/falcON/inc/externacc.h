// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// inc/externacc.h                                                             |
//                                                                             |
// Copyright (C) 2000-2005  Walter Dehnen                                      |
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
#ifndef falcON_included_externacc_h
#define falcON_included_externacc_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::acceleration                                               //
  //                                                                          //
  // is an abstract base class for a class that computes accelerations        //
  //                                                                          //
  // Notes:                                                                   //
  // 1 if the pointer to flags is NULL, all bodies are supposed to be active, //
  //   otherwise only those for which (f[i] & 1) is true.                     //
  // 2 the last argument of accelerations::set() indicates whether the        //
  //   accelerations and potential shall be assigned or added.                //
  //   If bit 0 is set, the potential    is added, otherwise assigned,        //
  //   If bit 1 is set, the acceleration is added, otherwise assigned.        //
  //   So, 0 means both are assigned.                                         //
  // 3 if masses and velocities are not required (as indicated by the         //
  //   routines need_masses() and need_velocities()), NULL pointers may be    //
  //   passed.                                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class acceleration {
    //--------------------------------------------------------------------------
    // types for vector in float or double                                      
    //--------------------------------------------------------------------------
  protected:
    typedef tupel<3,float>  vectf;
    typedef tupel<3,double> vectd;
    typedef bodies::block   block;
    //--------------------------------------------------------------------------
    // abstract methods to be satisfied by derived instances                    
    //--------------------------------------------------------------------------
  public:
    virtual bool is_empty()        const=0;        // R: no potential?          
    virtual bool need_masses()     const=0;        // R: masses needed?         
    virtual bool need_velocities() const=0;        // R: velocities needed?     
    //--------------------------------------------------------------------------
    // computing external gravity at a set of positions                         
    virtual
    void set(double      ,                         // I: simulation time        
	     int         ,                         // I: # bodies=size of tables
	     const float*,                         // I: masses                 
	     const vectf*,                         // I: positions              
	     const vectf*,                         // I: velocities             
	     const flags*,                         // I: flags                  
	     float      *,                         // O: potentials             
	     vectf      *,                         // O: accelerations          
	     int         )  const = 0;             // I: see note 2 above       
    //--------------------------------------------------------------------------
    virtual
    void set(double       ,                        // I: simulation time        
	     int          ,                        // I: # bodies=size of tables
	     const double*,                        // I: masses                 
	     const vectd *,                        // I: positions              
	     const vectd *,                        // I: velocities             
	     const flags *,                        // I: flags                  
	     double      *,                        // O: potentials             
	     vectd       *,                        // O: accelerations          
	     int          )  const = 0;            // I: see note 2 above       
    //--------------------------------------------------------------------------
    // non-virtual functions using the abstract methods above                   
    //--------------------------------------------------------------------------
    // adding external gravity to a snapshot                                    
    void set(const snapshot*snap,                  // I: bodies                 
	     bool           all,                   // I: all or active bodies?  
	     int            add = 2) const         //[I: see note 2 above]      
    {
      for(const block* b=snap->first_block(); b; b=b->next())
	set(snap->time(), b->N_bodies(),
	    need_masses()     ? b->const_data<fieldbit::m>() : 0,
	                        b->const_data<fieldbit::x>(),
	    need_velocities() ? b->const_data<fieldbit::v>() : 0,
	    all               ? 0 : b->const_data<fieldbit::f>(),
	                        b->data<fieldbit::q>(),
	                        b->data<fieldbit::a>(),
	    add);
    }
    // computing external gravity at a single position                          
    //                                                                          
    // WARNING  This routine makes sense only for static potentials.            
    void single_force(double       t,              // I: time                   
		      vectf  const&x,              // I: position               
		      float       &p,              // O: potential              
		      vectf       &a) const {      // O: acceleration           
      if(need_masses())
	warning("single_force(): accelertion requires m; force may be wrong");
      if(need_velocities())
	warning("single_force(): accelertion requires v; force may be wrong");
      const float m(one);
      const vectf v(zero);
      set(t,1,&m,&x,&v,0,&p,&a,0);
    }
    //                                                                          
    void single_force(double       t,              // I: time                   
		      vectd  const&x,              // I: position               
		      double      &p,              // O: potential              
		      vectd       &a) const {      // O: acceleration           
      if(need_masses())
	warning("single_force(): accelertion requires m; force may be wrong");
      if(need_velocities())
	warning("single_force(): accelertion requires v; force may be wrong");
      const double m(one);
      const vectd  v(zero);
      set(t,1,&m,&x,&v,0,&p,&a,0);
    }
    // emulating nemo_pot::pot_f()                                              
    // to enable the usage of class acceleration in place of class extpot       
    float pot_f      (vectf      &a,
		      vectf const&x,
		      double     t) const {
      float p;
      single_force(t,x,p,a);
      return p;
    }
    //                                                                          
    double pot_f     (vectd      &a,
		      vectd const&x,
		      double     t) const {
      double p;
      single_force(t,x,p,a);
      return p;
    }
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// forward declaration of some nemo stuff needed in construction of nemo_acc  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  void warning(char*, ...);
  typedef void (*pacc)(int,double,int,const void*,const void*,const void*,
		       const int*, void*, void*, int, char);
  pacc get_acceleration(const char*, const char*, const char*, bool*, bool*);
}
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::nemo_acc                                                   //
  //                                                                          //
  // non-abstract class, derived from acceleration                            //
  // loads external acceleration field at run time from .so file              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // A sensible usage of the constructor in a main() which has the options    //
  // accname, accpars, and accfile would be:                                  //
  // ...                                                                      //
  // acceleration *acc = hasvalue("accname")?                                 //
  //   new nemo_acc(getparam("accname"),                                      //
  //                hasvalue("accpars")? getparam("accpars") : 0,             //
  //                hasvalue("accfile")? getparam("accfile") : 0)             //
  //   : 0;                                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class nemo_acc: public acceleration {
  private:
    pacc ACC;
    bool NeedM, NeedV;
  public:
    //--------------------------------------------------------------------------
    bool is_empty() const { return ACC == 0; }
    bool need_masses() const { return ACC!=0 && NeedM; }
    bool need_velocities() const { return ACC!=0 && NeedV; }
    void set(double      t,
	     int         N,
	     const float*m,
	     const vectf*x,
	     const vectf*v,
	     const flags*f,
	     float      *p,
	     vectf      *a,
	     int         i) const
    {
      int ndim = Ndim;
      ACC(ndim,t,N,
	  static_cast<const void*>(m),
	  static_cast<const void*>(x),
	  static_cast<const void*>(v),
	  static_cast<const int*>(static_cast<const void*>(f)),
	  static_cast<void*>(p),
	  static_cast<void*>(a), i, 'f');
    }
    //--------------------------------------------------------------------------
    void set(double       t,
	     int          N,
	     const double*m,
	     const vectd *x,
	     const vectd *v,
	     const flags *f,
	     double      *p,
	     vectd       *a,
	     int          i) const
    {
      int ndim = Ndim;
      ACC(ndim,t,N,
	  static_cast<const void*>(m),
	  static_cast<const void*>(x),
	  static_cast<const void*>(v),
	  static_cast<const int*>(static_cast<const void*>(f)),
	  static_cast<void*>(p),
	  static_cast<void*>(a), i, 'd');
    }
    //--------------------------------------------------------------------------
    nemo_acc(const char*accname,
	     const char*accpars,
	     const char*accfile) : 
      ACC ( get_acceleration(accname,accpars,accfile,&NeedM,&NeedV) ) {}
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_externacc_h
