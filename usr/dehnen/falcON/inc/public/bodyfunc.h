// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/bodyfunc.h                                               
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2008 Walter Dehnen                                        
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
#ifndef falcON_included_bodyfunc_h
#define falcON_included_bodyfunc_h

#ifdef falcON_NEMO
#ifndef falcON_included_body_h
#  include <body.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary class bf_type<>                                                  
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct bf_type_base {};
  template<> struct bf_type_base<bool> {
    static const char type = 'b';
    static const char*name() { return "boolean"; }
    static  bool default_value() { return true; }
  };
  template<> struct bf_type_base<int> {
    static const char type = 'i';
    static const char*name() { return "integer"; }
    static  int default_value() { return 0; }
  };
  template<> struct bf_type_base<real> {
    static const char type = 'r';
    static const char*name() { return "scalar"; }
    static  real default_value() { return zero; }
  };
  template<> struct bf_type_base<vect> {
    static const char type = 'v';
    static const char*name() { return "vector"; }
    static  vect default_value() { return vect(zero); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct bf_type : public bf_type_base<T> {
    typedef T(*const bf_pter)(body   const&,double const&,const real*);
    static T func(void(*const f)(body const&,double const&,const real*),
		  body   const&b, double const&t, const real*p) {
      return f? (*(bf_pter)(f))(b,t,p) : bf_type_base<T>::default_value();
    }
    //--------------------------------------------------------------------------
    typedef T(*const Bf_pter)(bodies const&,double const&,const real*);
    static T func(void(*const f)(bodies const&,double const&,const real*),
		  bodies const&b, double const&t, const real*p) {
      return f? (*(Bf_pter)(f))(b,t,p) : bf_type_base<T>::default_value();
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class bodyfunc                                                             
  //                                                                            
  /// a function object taking a body, time, and possibly parameters, and       
  /// returning bool, int, real, or vect                                        
  ///                                                                           
  /// A bodyfunc is constructed from an expression (C-style string, see man     
  /// page (5bodyfunc)), which is used to generate a function (employing the    
  /// compiler) and loading it at run time. \n                                  
  /// A database is used to store expression and functions, which makes         
  /// repeated access to the same function very quick.                          
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class bodyfunc {
  protected:
    void      (*FUNC) (body const&, double const&, const real*);
    char        TYPE;
    int         NPAR;
    fieldset    NEED;
    char       *EXPR;
    //--------------------------------------------------------------------------
  public:
    /// maximum number of parameters allowed
    /// \note MAXPAR=10 is hardwired into bodyfunc.cc: only one digit is allowed
    static const int MAXPAR = 10;
    /// print info about bodyfuncs in database, if any
    /// \return true if something was printed out
    /// \param  out  ostream to print to
    static bool print_db(std::ostream&out);
    /// ctor from bodyfunc expression (see man pages)
    explicit bodyfunc(const char*) throw(falcON::exception);
    /// dtor: delete data
    ~bodyfunc() { if(EXPR) falcON_DEL_A(EXPR); EXPR=0; }
    /// return type: 'b', 'i', 'r', 'v' for bool, int, real, vect
    char const&type() const { return TYPE; }
    /// check return type
    template<typename T>
    bool is_type() const { return bf_type<T>::type == TYPE; }
    /// number of parameters
    int  const&npar() const { return NPAR; }
    /// return data needed
    fieldset const&need() const { return NEED; }
    /// do we need this datum?
    bool need(fieldbit b) const { return NEED.contain(b); }
    /// do we need this datum?
    bool need(fieldbit::bits b) const { return NEED.contain(b); }
    /// return original expression
    const char* expression() const { return EXPR; }
    /// function call, non-operator
    /// \param b body
    /// \param t time
    /// \param p parameters
    /// \return expression evaluated for body \a b, time \a t, parameters \a p
    template<typename T>
    T func (body const&b, double const&t, const real*p) const {
#if defined(DEBUG) || defined(EBUG)
      if(bf_type<T>::type != TYPE)
	falcON_THROW("bodyfunc::func_safe<%s>() called, but type is %s\n",
		     nameof(T), (
		     TYPE == 'b'? "bool" :
		     TYPE == 'i'? "int"  :
		     TYPE == 'r'? "real" :
		     TYPE == 'v'? "vect" : "unknown") );
      if(!b)
	falcON_THROW("bodyfunc::func_safe<%s>() called on invalid body\n",
		     nameof(T));
      if(!b.my_bodies()->have_all(NEED))
	falcON_THROW("bodyfunc::func_safe<%s>() not all data need (%s) "
		     "known for body\n", word(NEED));
#endif
      return bf_type<T>::func(FUNC,b,t,p);
    }
    /// function call operator (not very useful, since template)
    template<typename T>
    T operator() (body const&b, double const&t, const real*p) const {
      return func<T>(b,t,p);
    }
    /// is *this an empty bodyfunc?
    bool is_empty() const { return FUNC == 0; }
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class Bodyfunc                                                             
  //                                                                            
  /// a bodyfunc with parameters                                                
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class Bodyfunc : protected bodyfunc {
    real P[MAXPAR];
    char*PARS;
  public:
    /// construction from bodyfunc expression (can be empty)
    /// \param expr bodyfunc (5falcON) expression --- or NULL
    /// \param pars comma separated list of parameters
    Bodyfunc(const char*expr, const char*pars)
      throw(falcON::exception);
    /// construction from bodyfunc expression (can be empty)
    /// \param expr bodyfunc (5falcON) expression --- or NULL
    /// \param pars array with parameters 
    /// \param npar number of parameters
    /// \note there must be enough parameters given
    Bodyfunc(const char*expr, const real*pars, int npar)
      throw(falcON::exception);
    /// dtor: delete data
    ~Bodyfunc() { if(PARS) falcON_DEL_A(PARS); PARS=0; }
    /// return type: 'b', 'i', 'r', 'v' for bool, int, real, vect
    bodyfunc::type;
    /// check return type
    bodyfunc::is_type;
    /// return number of parameters used
    bodyfunc::npar;
    /// return nth parameter
    /// \param  n number of parameter asked
    real const&param(int n) const { return P[n]; }
    /// return fields required
    bodyfunc::need;
    /// return original expression
    bodyfunc::expression;
    /// return parameters
    const char*parameters() const { return PARS; }
    /// function call, non-operator
    /// \param b body
    /// \param t time
    /// \return expression evaluated for body \a b at time \a t
    template<typename T>
    T func (body const&b, double const&t) const {
      return bodyfunc::func<T>(b,t,P);
    }
    /// function call operator (not very useful, since template)
    template<typename T>
    T operator() (body const&b, double const&t) const {
      return func<T>(b,t);
    }
    /// is *this an empty bodyfunc?
    bodyfunc::is_empty;
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class BodyFunc<T>                                                          
  //                                                                            
  /// a function object taking body and time, returning T (template parameter)  
  ///                                                                           
  /// essentially this is a wrapper around class bodyfunc with the following    
  /// distinctions: \n                                                          
  /// - the return type of the bodyfunc expression \b must match T            \n
  /// - the parameters are read by the constructor, ie. cannot be varied later\n
  /// - if the bodyfunc is empty, we return: true, 0 , zero, vect(zero).        
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> class BodyFunc : private Bodyfunc {
  public:
    /// construction from bodyfunc expression (can be empty)
    /// \param expr body_func (5falcON) expression --- or NULL
    /// \param pars comma separated list of parameters
    /// \note the bodyfunc expression must return the type T
    BodyFunc(const char*expr, const char*pars)
      throw(falcON::exception);
    /// construction from bodyfunc expression (can be empty)
    /// \param expr body_func (5falcON) expression --- or NULL
    /// \param pars array with parameters 
    /// \param npar number of parameters
    /// \note there must be enough parameters given
    /// \note the bodyfunc expression must return the type T
    BodyFunc(const char*expr, const real*pars, int npar)
      throw(falcON::exception);
    /// return number of parameters used
    Bodyfunc::npar;
    /// return nth parameter
    /// \param  n number of parameter asked
    Bodyfunc::param;
    /// return fields required
    Bodyfunc::need;
    /// return original expression
    Bodyfunc::expression;
    /// return parameters
    Bodyfunc::parameters;
    /// function call
    /// \param b body
    /// \param t time
    /// \return expression evaluated for body \a b at time \a t
    T operator()(body const&b, double const&t) const {
      return Bodyfunc::func<T>(b,t);
    }
    /// is *this an empty bodyfunc?
    Bodyfunc::is_empty;
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  template<int BIT> struct BodyPropMap {
    typedef typename field_traits<BIT>::type type;
  };
  template<> struct BodyPropMap<fieldbit::f> { typedef int  type; };
  template<> struct BodyPropMap<fieldbit::k> { typedef int  type; };
  template<> struct BodyPropMap<fieldbit::s> { typedef real type; };
  template<> struct BodyPropMap<fieldbit::l> { typedef int  type; };
  template<> struct BodyPropMap<fieldbit::n> { typedef int  type; };
  template<> struct BodyPropMap<fieldbit::N> { typedef int  type; };
  template<> struct BodyPropMap<fieldbit::c> { typedef int  type; };
  template<> struct BodyPropMap<fieldbit::h> { typedef int  type; };
  // ---------------------------------------------------------------------------
  template<int BIT> struct BodyPropType {
    typedef typename field_traits<BIT>::type proptype;
    typedef typename BodyPropMap <BIT>::type functype;
    static proptype convert(functype const&x) {
      return static_cast<proptype const&>(x);
    }
  };
  template<> struct BodyPropType<fieldbit::f> {
    typedef field_traits<fieldbit::f>::type proptype;
    typedef BodyPropMap <fieldbit::f>::type functype;
    static proptype convert(functype const&x) {
      return *(static_cast<const flags*>(static_cast<const void*>(&x)));
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class BodyProp<BIT>                                                        
  //                                                                            
  /// wrapper around Bodyfunc, returns field_traits<BIT>::type                  
  ///                                                                           
  /// If the return type of bodyfunc is incompatible with our return type, an   
  /// exception is thrown.                                                      
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  template<int BIT> class BodyProp : private BodyPropType<BIT> {
    const Bodyfunc*F;   ///< Bodyfunc
    const double   T;   ///< time
  public:
    /// type returned by operator()
    typedef typename BodyPropType<BIT>::proptype proptype;
    /// type required by bodyfunc
    typedef typename BodyPropType<BIT>::functype functype;
    /// construction from bodyfunc and parameters
    /// \param f pter to bodyfunc
    /// \param p array with parameters (may be empty)
    /// \param n size of array p
    /// \note n must match b->npar() (we issue a warning if it exceeds)
    BodyProp(const Bodyfunc*f, double t) 
      throw(falcON::exception) : F(f), T(t)
    {
      if(!F->is_type<functype>())
	throw falcON::exception("BodyProp<%c>: type mismatch: "
				"bodyfunc returns %s, but need %s",
				field_traits<BIT>::word,
				F->type() == 'b'? "bool" :
				F->type() == 'i'? "int" :
				F->type() == 'r'? "real" : "vect",
				nameof(functype));
    }
    /// number of parameters
    int  const&npar() const { return F->npar(); }
    /// return nth parameter
    real const&param(int n) const { return F->param(n); }
    /// return fields required
    fieldset const&need() const { return F->need(); }
    /// return original expression
    const char* expression() const { return F->expression(); }
    /// return parameters
    const char*parameters() const { return F->parameters(); }
    /// function call
    /// \param b body
    /// \param t time
    /// \return expression evaluated for body \a b at time \a t
    proptype operator()(body const&b) const {
      return convert(F->func<functype>(b,T));
    }
    /// is *this an empty bodyfunc?
    bool is_empty() const { return F->is_empty(); }
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class BodyFilter                                                           
  //                                                                            
  /// a function object taking body, returning bool                             
  ///                                                                           
  /// essentially this is a wrapper around class BodyFunc<bool> with the        
  /// distinctions that simulation time is stored separately.                   
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class BodyFilter : private BodyFunc<bool> {
    double TIME;           ///< simulation time at which operator() is evaluated
  public:
    /// construction from bodyfunc expression (can be empty)
    /// \param expr body_func (5falcON) expression --- or NULL
    /// \param pars comma separated list of parameters
    /// \note the bodyfunc expression must return the type T
    BodyFilter(const char*expr, const char*pars)
      throw(falcON::exception);
    /// construction from bodyfunc expression (can be empty)
    /// \param expr body_func (5falcON) expression --- or NULL
    /// \param pars array with parameters 
    /// \param npar number of parameters
    /// \note there must be enough parameters given
    /// \note the bodyfunc expression must return the type T
    BodyFilter(const char*expr, const real*pars, int npar)
      throw(falcON::exception);
    /// set the time for future function calls
    /// \param t (input) simulation time
    void set_time(double t) { TIME = t; }
    /// return number of parameters used
    BodyFunc<bool>::npar;
    /// return nth parameter
    /// \param  n number of parameter asked
    BodyFunc<bool>::param;
    /// return fields required
    BodyFunc<bool>::need;
    /// return original expression
    BodyFunc<bool>::expression;
    /// return parameters
    BodyFunc<bool>::parameters;
    /// function call
    /// \param b body
    /// \return expression evaluated for body \a b at time set by set_time()
    bool operator()(body const&b) const {
      return BodyFunc<bool>::operator() (b,TIME);
    }
    /// is *this an empty bodyfunc?
    BodyFunc<bool>::is_empty;
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
    /// return number of bodies passing the filter
    /// \param B (input) bodies
    unsigned N_bodies(const bodies*B) const {
      unsigned n = 0;
      LoopAllBodies(B,b) if(operator()(b)) ++n;
      return n;
    }
    /// return first filtered body
    body first(const bodies*B) const {
      LoopAllBodies(B,b)
	if(operator()(b)) return b;
      return bodies::bodyNil();
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class bodiesfunc                                                           
  //                                                                            
  /// a function object taking bodies, time, and possibly parameters, and       
  /// returning bool, int, real, or vect                                        
  ///                                                                           
  /// A bodyfunc is constucted from an expression (C-style string), which is    
  /// used to generate a function (employing the compiler) and loading it at    
  /// run time. \n                                                              
  /// A database is used to store expression and functions, which makes         
  /// repeated access to the same function very quick.                          
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class bodiesfunc {
    void   (*FUNC) (bodies const&, double const&, const real*);
    char     TYPE;
    int      NPAR;
    fieldset NEED;
    //--------------------------------------------------------------------------
  public:
    /// print info about bodyfuncs in database, if any
    /// \return true if something was printed out
    /// \param  out  ostream to print to
    static bool print_db(std::ostream&out);
    /// construction from bodyfunc expression
    explicit bodiesfunc(const char*) throw(falcON::exception);
    /// return type: 'b', 'i', 'r', 'v' for bool, int, real, vect
    char const&type() const { return TYPE; }
    /// check return type
    template<typename T>
    bool is_type() const { return bf_type<T>::type == TYPE; }
    /// number of parameters
    int  const&npar() const { return NPAR; }
    /// return data needed
    fieldset const&need() const { return NEED; }
    /// non-operator function call with bodies, time, and parameters
    template<typename T>
    T func (bodies const&b, double const&t, const real*p) const {
#if defined(DEBUG) || defined(EBUG)
      if(bf_type<T>::type != TYPE)
	falcON_THROW("bodiesfunc::func<%s>() called, but type is %s\n",
		     nameof(T), (
		     TYPE == 'b'? "bool" :
		     TYPE == 'i'? "int"  :
		     TYPE == 'r'? "real" :
		     TYPE == 'v'? "vect" : "unknown") );
      if(!b.have_all(NEED))
	falcON_THROW("bodiesfunc::func<%s>() not all data need (%s) "
		     "known for body\n", word(NEED));
#endif
      return bf_type<T>::func(FUNC,b,t,p);
    }
    /// non-operator function call with snapshot and parameters
    template<typename T>
    T func (snapshot const&s, const real*p) const {
      return func<T>(static_cast<const bodies&>(s),s.time(),p);
    }
    /// function call operator (not very useful, since template)
    template<typename T>
    T operator() (bodies const&b, double const&t, const real*p) const {
      return func<T>(b,t,p);
    }
    /// function call operator (not very useful, since template)
    template<typename T>
    T operator() (snapshot const&s, const real*p) const {
      return func<T>(s,p);
    }
    /// is *this an empty bodiesfunc?
    bool is_empty() const { return FUNC == 0; }
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class bodiesmethod                                                         
  //                                                                            
  /// a function object taking a pointer to anything (void*), bodies, time, and 
  /// possibly parameters (for example reading the body angular momenta into an 
  /// array).                                                                   
  ///                                                                           
  /// A bodiesmethod is constucted from an expression (C-style string), which   
  /// is used to generate a function (employing the compiler) and loading it at 
  /// run time. \n                                                              
  /// A database is used to store expression and functions, which makes         
  /// repeated access to the same function very quick.                          
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class bodiesmethod {
    void   (*FUNC) (void*, bodies const&, double const&, const real*);
    char     TYPE;
    int      NPAR;
    fieldset NEED;
    //--------------------------------------------------------------------------
    template<typename T>
    void func (T*d, bodies const&b, double const&t, const real*p) const {
#if defined(DEBUG) || defined(EBUG)
      if(bf_type<T>::type != TYPE)
	falcON_THROW("bodiesmethod::func<%s>() called, but type is %s\n",
		     nameof(T), (
		     TYPE == 'b'? "bool" :
		     TYPE == 'i'? "int"  :
		     TYPE == 'r'? "real" :
		     TYPE == 'v'? "vect" : "unknown") );
      if(!b.have_all(NEED))
	falcON_THROW("bodiesfunc::func<%s>() not all data need (%s) "
		     "known for body\n", word(NEED));
#endif
      return FUNC(static_cast<void*>(d),b,t,p);
    }
    //--------------------------------------------------------------------------
  public:
    /// construction from bodyfunc expression
    explicit bodiesmethod(const char*) falcON_THROWING;
    /// return type: 'b', 'i', 'r', 'v' for bool, int, real, vect
    char const&type() const { return TYPE; }
    /// check return type
    template<typename T>
    bool is_type() const { return bf_type<T>::type == TYPE; }
    /// number of parameters
    int  const&npar() const { return NPAR; }
    /// return data needed
    fieldset const&need() const { return NEED; }
    /// function call
    template<typename T>
    void operator() (T*d, bodies const&b, double const&t, const real*p) const {
      return func(d,b,t,p);
    }
    /// is *this an empty bodiesmethod?
    bool is_empty() const { return FUNC == 0; }
    /// is *this valid (non-empty)?
    operator bool() const { return !is_empty(); }
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  template<> struct traits<falcON::bodiesmethod > {
    static const char  *name () { return "bodiesmethod"; }
  };
}
#endif // falcON_PROPER
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  template<> struct traits<falcON::bodyfunc > {
    static const char  *name () { return "bodyfunc"; }
  };
  template<typename T> struct traits<falcON::BodyFunc<T> > {
    static const char  *name () {
      return message("BodyFunc<%s>",traits<T>::name()); }
  };
  template<> struct traits<falcON::bodiesfunc > {
    static const char  *name () { return "bodiesfunc"; }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
#endif // falcON_included_bodyfunc_h

 
