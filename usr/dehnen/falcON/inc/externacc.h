// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/externacc.h                                                     
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2007                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2007 Walter Dehnen                                        
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_externacc_h
#define falcON_included_externacc_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::acceleration                                                 
  //                                                                            
  /// is an abstract base class for a class that computes accelerations         
  ///                                                                           
  /// \note                                                                     
  ///  if the pointer to flags is NULL, all bodies are supposed to be active,   
  ///  otherwise only those for which (f[i] & 1) is true.                       
  /// \note                                                                     
  ///  the last argument of accelerations::set() indicates whether the          
  ///  accelerations and potential shall be assigned or added.\n                
  ///  If bit 0 is set, the potential    is added, otherwise assigned,\n        
  ///  If bit 1 is set, the acceleration is added, otherwise assigned.\n        
  ///  So, 0 means both are assigned.                                           
  /// \note                                                                     
  ///  if masses and velocities are not required (as indicated by the           
  ///  routines need_masses() and need_velocities()), NULL pointers may be      
  ///  passed.                                                                  
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class acceleration {
  protected:
    typedef tupel<3,float>  vectf;
    typedef tupel<3,double> vectd;
    typedef bodies::block   block;
    //--------------------------------------------------------------------------
    /// \name abstract methods to be satisfied by derived instances
    //@{
    //--------------------------------------------------------------------------
  public:
    virtual bool is_empty()        const=0;        ///< no potential?
    virtual bool need_masses()     const=0;        ///< masses needed?
    virtual bool need_velocities() const=0;        ///< velocities needed?
    //--------------------------------------------------------------------------
    /// computing external gravity at a set of positions
    ///
    /// The parameter \a add indicates whether the accelerations and potential 
    /// shall be assigned or added. If bit 0 is set, the potential is added,
    /// otherwise assigned; if bit 1 is set, the acceleration is added,
    /// otherwise assigned. So, 0 means both are assigned.
    /// \param t (input) simulation time
    /// \param n (input) # bodies=size of tables
    /// \param m (input) table with masses
    /// \param x (input) table with positions
    /// \param v (input) table with velocities
    /// \param f (input) table with flags
    /// \param m (output) table with potentials
    /// \param m (output) table with accelerations
    /// \param add (input) see detailed description
    virtual
    void set(double      t,
	     int         n,
	     const float*m,
	     const vectf*x,
	     const vectf*v,
	     const flags*f,
	     float      *p,
	     vectf      *a,
	     int         add)  const = 0;
    //--------------------------------------------------------------------------
    /// computing external gravity at a set of positions
    ///
    /// The parameter \a add indicates whether the accelerations and potential 
    /// shall be assigned or added. If bit 0 is set, the potential is added,
    /// otherwise assigned; if bit 1 is set, the acceleration is added,
    /// otherwise assigned. So, 0 means both are assigned.
    /// \param t (input) simulation time
    /// \param n (input) # bodies=size of tables
    /// \param m (input) table with masses
    /// \param x (input) table with positions
    /// \param v (input) table with velocities
    /// \param f (input) table with flags
    /// \param m (output) table with potentials
    /// \param m (output) table with accelerations
    /// \param add (input) see detailed description
    virtual
    void set(double       t,
	     int          n,
	     const double*m,
	     const vectd *x,
	     const vectd *v,
	     const flags *f,
	     double      *p,
	     vectd       *a,
	     int          add)  const = 0;
    //@}
    //--------------------------------------------------------------------------
    /// \name non-virtual functions using the abstract methods above
    //@{
    //--------------------------------------------------------------------------
    /// adding external gravity to a snapshot
    ///
    /// The parameter \a add indicates whether the accelerations and potential 
    /// shall be assigned or added. If bit 0 is set, the potential is added,
    /// otherwise assigned; if bit 1 is set, the acceleration is added,
    /// otherwise assigned. So, 0 means both are assigned, while the default is
    /// to add acceleration but assign potential to field pex (see fieldset).
    /// \param snap (input) snapshot = time + bodies
    /// \param all  (input) set gravity for all bodies or only active bodies?
    /// \param add  (optional input) see detailed description
    void set(const snapshot*snap,
	     bool           all,
	     int            add = 2) const
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
    //--------------------------------------------------------------------------
    /// computing external gravity at a single position
    ///
    /// \param t (input) simulation time
    /// \param x (input) position
    /// \param p (output) potential
    /// \param a (output) acceleration
    /// \warning  This routine makes sense only for static potentials.
    void single_force(double       t,
		      vectf  const&x,
		      float       &p,
		      vectf       &a) const {
      if(need_masses())
	falcON_Warning("single_force(): "
		       "accelertion requires m; force may be wrong");
      if(need_velocities())
	falcON_Warning("single_force(): "
		       "accelertion requires v; force may be wrong");
      const float m(one);
      const vectf v(zero);
      set(t,1,&m,&x,&v,0,&p,&a,0);
    }
    //--------------------------------------------------------------------------
    /// computing external gravity at a single position
    ///
    /// \param t (input) simulation time
    /// \param x (input) position
    /// \param p (output) potential
    /// \param a (output) acceleration
    /// \warning  This routine makes sense only for static potentials.
    void single_force(double       t,
		      vectd  const&x,
		      double      &p,
		      vectd       &a) const {
      if(need_masses())
	falcON_Warning("single_force(): "
		       "accelertion requires m; force may be wrong");
      if(need_velocities())
	falcON_Warning("single_force(): "
		       "accelertion requires v; force may be wrong");
      const double m(one);
      const vectd  v(zero);
      set(t,1,&m,&x,&v,0,&p,&a,0);
    }
    //--------------------------------------------------------------------------
    /// emulating nemo_pot::pot_f() to enable the usage of class acceleration
    /// in place of class extpot
    ///
    /// \return potential
    /// \param a (output) acceleration
    /// \param x (input) position
    /// \param t (input) simulation time
    float pot_f      (vectf      &a,
		      vectf const&x,
		      double     t) const {
      float p;
      single_force(t,x,p,a);
      return p;
    }
    //--------------------------------------------------------------------------
    /// emulating nemo_pot::pot_f() to enable the usage of class acceleration
    /// in place of class extpot
    ///
    /// \return potential
    /// \param a (output) acceleration
    /// \param x (input) position
    /// \param t (input) simulation time
    double pot_f     (vectd      &a,
		      vectd const&x,
		      double     t) const {
      double p;
      single_force(t,x,p,a);
      return p;
    }
    //@}
    /// noon dtor
    virtual~acceleration() {}
    //--------------------------------------------------------------------------
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::sum_acceleration                                             
  //                                                                            
  /// the sum of two acceleration fields.                                       
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class sum_acceleration : public acceleration {
  private:
    typedef const acceleration *cp_accel;
    cp_accel A1,A2;
  public:
    sum_acceleration(cp_accel a1, cp_accel a2) : A1(a1), A2(a2) {}
    //--------------------------------------------------------------------------
    /// \name implement abstract methods of base class
    //@{
  public:
    /// no potential?
    bool is_empty() const {
      return A1->is_empty() && A2->is_empty();
    }
    /// masses needed?
    bool need_masses() const {
      return A1->need_masses() || A2->need_masses();
    }
    /// velocities needed?
    bool need_velocities() const {
      return A1->need_velocities() || A2->need_velocities();
    }
    /// computing external gravity at a set of positions
    void set(double t, int n, const float*m, const vectf*x, const vectf*v,
	     const flags*f, float*p, vectf*a, int add) const
    {
      if(!A1->is_empty()) {
	A1->set(t,n,m,x,v,f,p,a,add);
	if(!A2->is_empty())
	  A2->set(t,n,m,x,v,f,p,a,3);
      } else {
	if(!A2->is_empty())
	  A2->set(t,n,m,x,v,f,p,a,add);
      }
    }
    /// computing external gravity at a set of positions
    void set(double t, int n, const double*m, const vectd*x, const vectd*v,
	     const flags*f, double*p, vectd*a, int add) const
    {
      if(!A1->is_empty()) {
	A1->set(t,n,m,x,v,f,p,a,add);
	if(!A2->is_empty())
	  A2->set(t,n,m,x,v,f,p,a,3);
      } else {
	if(!A2->is_empty())
	  A2->set(t,n,m,x,v,f,p,a,add);
      }
    }
    //@}
    /// noon dtor
    virtual~sum_acceleration() {}
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
//   void warning(char*, ...);   // nemo's warning
  typedef void (*pacc)(int,double,int,const void*,const void*,const void*,
		       const int*, void*, void*, int, char);
  pacc get_acceleration(const char*, const char*, const char*, bool*, bool*);
}
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::nemo_acc                                                     
  //                                                                            
  /// non-abstract class, derived from acceleration                             
  /// loads external acceleration field at run time from .so file               
  ///                                                                           
  /// A sensible usage of the constructor in a main() which has the options     
  /// accname, accpars, and accfile would be: \code                             
  /// acceleration *acc = hasvalue("accname")?                                  
  ///   new nemo_acc(getparam("accname"),                                       
  ///                hasvalue("accpars")? getparam("accpars") : 0,              
  ///                hasvalue("accfile")? getparam("accfile") : 0)              
  ///   : 0;                                                                    
  /// \endcode                                                                  
  //////////////////////////////////////////////////////////////////////////////
  class nemo_acc: public acceleration {
  private:
    pacc ACC;
    bool NeedM, NeedV;
  public:
    //--------------------------------------------------------------------------
    /// is there any external acceleration at all?
    bool is_empty() const { return ACC == 0; }
    //--------------------------------------------------------------------------
    /// do we need body masses?
    bool need_masses() const { return ACC!=0 && NeedM; }
    //--------------------------------------------------------------------------
    /// do we need body velocities?
    bool need_velocities() const { return ACC!=0 && NeedV; }
    //--------------------------------------------------------------------------
    /// computing external gravity at a set of positions
    ///
    /// The parameter \a i indicates whether the accelerations and potential 
    /// shall be assigned or added. If bit 0 is set, the potential is added,
    /// otherwise assigned; if bit 1 is set, the acceleration is added,
    /// otherwise assigned. So, 0 means both are assigned.
    /// \param t (input) simulation time
    /// \param n (input) # bodies=size of tables
    /// \param m (input) table with masses
    /// \param x (input) table with positions
    /// \param v (input) table with velocities
    /// \param f (input) table with flags
    /// \param m (output) table with potentials
    /// \param m (output) table with accelerations
    /// \param i (input) see detailed description
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
    /// computing external gravity at a set of positions
    ///
    /// The parameter \a i indicates whether the accelerations and potential 
    /// shall be assigned or added. If bit 0 is set, the potential is added,
    /// otherwise assigned; if bit 1 is set, the acceleration is added,
    /// otherwise assigned. So, 0 means both are assigned.
    /// \param t (input) simulation time
    /// \param n (input) # bodies=size of tables
    /// \param m (input) table with masses
    /// \param x (input) table with positions
    /// \param v (input) table with velocities
    /// \param f (input) table with flags
    /// \param m (output) table with potentials
    /// \param m (output) table with accelerations
    /// \param i (input) see detailed description
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
    /// ctor from name, parameters, and file
    ///
    /// The constructor will (try to) load a .so file of name given in \a
    /// accname and (try to) initialize an acceleration field using the 
    /// parameters (if any) provided in \a accpars and the data file (if any)
    /// provided in \a accfile.
    /// \param accname name of the external shared object file
    /// \param accpars parameters (comma-separated list of values)
    /// \param accfile data file possibly required by the external potential
    nemo_acc(const char*accname,
	     const char*accpars,
	     const char*accfile) : 
      ACC ( get_acceleration(accname,accpars,accfile,&NeedM,&NeedV) ) {}
    /// noon dtor
    virtual~nemo_acc() {}
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_externacc_h
