/*
 *****************************************************************************
 *
 * defacc.h
 *
 *  declaration of functions that any implementation of an external,
 *  dynamically loadable acceleration must define. See also acceleration.h
 *
 *  to ensure that your implementation of external acceleration is consistent,
 *  include this header file into your C or C++ source code.
 *
 *  Notes:
 *  1 the 4th argument to iniacceleration() is boolean and returns whether
 *    masses are required as input for acceleration().
 *  2 the 5th argument to iniacceleration() is boolean and returns whether
 *    velocities are required as input for acceleration().
 *    Velocities may be used to compute friction forces, such as the drag
 *    a gaseous disk is generating on stars crossing it.
 *  3 arrays are passed to acceleration() as pointer to void. They must be
 *    either all of type float or all of type double as indicated by the last
 *    argument being 'f' or 'd', respectively.
 *  4 arrays of vector quantities are in the order x0,y0,z0, x1,y1,z1, ...
 *  5 if the pointer to flags is NULL, all bodies are supposed to be active,
 *    otherwise only those for which (f[i] & 1) is true.
 *  6 the argument "indicator" of acceleration() indicates whether the
 *    accelerations and potential shall be assigned or added.
 *    If bit 0 is set, the potential    is added, otherwise assigned,
 *    If bit 1 is set, the acceleration is added, otherwise assigned.
 *    So, 0 means both are assigned.
 *
 ******************************************************************************
 * version 0.0  17/06/2004  WD
 * version 1.0  22/06/2004  WD  support for acc & pot via C++
 *
 */
#ifndef included_defacc_h
#define included_defacc_h

#ifdef __cplusplus
# include <cstring>
  extern "C" {
#else
# include <string.h>
#endif

/*
 * functions that may be used in the implementation of the above.
 * they are resolved in the NEMO library
 * (this way, we avoid including the nemo.h)
 */

void warning(char *, ...);
void error  (char *, ...);
int  nemo_dprintf(int, const char *, ...);

/*
 * functions that to be defined in #including file.
 */

void iniacceleration(                 /* return: void                         */
		     const double*,   /* input:  array with parameters        */
		     int,             /* input:  number of parameters         */
		     const char*,     /* input:  data file name (may be NULL) */
		     bool*,           /* output: acceleration() needs masses? */
		     bool*);          /* output: acceleration() needs vel's?  */

void acceleration   (                 /* return: void                         */
		     int,             /* input:  number of dimensions         */
		     double,          /* input:  simulation time              */
		     int,             /* input:  number bodies =size of arrays*/
		     const void*,     /* input:  masses:         m[i]         */
		     const void*,     /* input:  positions       (x,y,z)[i]   */
		     const void*,     /* input:  velocities      (u,v,w)[i]   */
		     const int *,     /* input:  flags           f[i]         */
		     void*,           /* output: potentials      p[i]         */
		     void*,           /* output: accelerations   (ax,ay,az)[i]*/
		     int,             /* input:  indicator (see note 6 above) */
		     char);           /* input:  type: 'f' or 'd'             */

#ifdef POT_DEF
/*
 * functions that may be defined in #including file.
 */
void inipotential    (                /* return: void                         */
		      const int*,     /* input:  number of parameters         */
		      const double*,  /* input:  array with parameters        */
		      const char*);   /* input:  data file name (may be NULL) */
void potential_double(const int*,     /* input:  number of dimensions         */
		      const double*,  /* input:  position                     */
		      double*,        /* output: acceleration                 */
		      double*,        /* output: potential                    */
		      const double*); /* input:  time                         */
void potential_float (const int*,     /* input:  number of dimensions         */
		      const float*,   /* input:  position                     */
		      float*,         /* output: acceleration                 */
		      float*,         /* output: potential                    */
		      const float*);  /* input:  time                         */
#endif

#ifdef __cplusplus
}
////////////////////////////////////////////////////////////////////////////////
//
// template and macros that ease the implementation of the above routines
// in C++
//
// needed is a class with the following properties
//
// class NAME {
//   public:
//   static const char* name();
//   NAME(const double*pars,
//        int          npar,
//        const char  *file);
//   bool NeedMass() const;
//   bool NeedVels() const;
//   inline void set_time(double time) const
//   template<int NDIM, typename scalar>
//   inline void acc(const scalar*mas,   // pter may be 0 if NeedMass()==0
// 	             const scalar*vel,   // pter may be 0 if NeedVels()==0
// 	             const scalar*pos,
// 	             scalar      &pot,
// 	             scalar      *acc) const;
// };
//

namespace {
  template<class Acceleration> class AccInstall {
  private:
    Acceleration *Acc;         // our acceleration
    double        Time;        // last/actual simulation time
    bool          First;       // true only after initialization
    double       *Pars;        // parameters of last call to init()
    int           Npar;        // # ---
    char         *File;        // data file of  last call to init()

    bool differ(const double*pars, int npar, const char*file)
      // check whether actual parameters & data file differ from last
    {
      if(Acc ==0) return true;
      if(File==0 && file!=0) return true;
      if(File!=0 && file==0) return true;
      if(File && file && strcmp(File,file)) return true;
      if(npar != Npar) return true;
      for(int n=0; n!=npar; ++n)
	if(pars[n] != Pars[n]) return true;
      return false;
    }

    void remember(const double*pars, int npar, const char*file)
      // remember parameters and data file
    {
      if(Pars) delete[] Pars;
      if(pars && npar>0) {
	Pars = new double[npar];
	for(int n=0; n!=npar; ++n) Pars[n] = pars[n];
      } else
	Pars = 0;
      Npar = npar;
      if(File) delete[] File;
      if(file) {
	size_t n = strlen(file)+1;
	File = new char[n];
	strncpy(File,file,n);
      } else
	File = 0;
    }

    void set_time(double time)
      // if time has changed or this is the first call, Acc->set_time()
    {
      if(First || Time != time) {
	Time  = time;
	First = false;
	if(Acc) Acc->set_time(Time);
      }
    }

    template<typename scalar>
    void acc_T(int          ndim,
	       double       time,
	       int          nbod,
	       const scalar*mas,
	       const scalar*pos,
	       const scalar*vel,
	       const int   *flag,
	       scalar      *pot,
	       scalar      *acc,
	       int          add)
      // compute the pot&acc and add/set pot/acc
    {

      // in case of adding, create arrays to write pot & acc into
      scalar*_pot = add&1 ? new scalar[nbod]      : 0;
      scalar*_acc = add&2 ? new scalar[nbod*ndim] : 0;

      // define references to the arrays actually passed to grav::get()
      scalar*&pots = add&1 ? _pot : pot;
      scalar*&accs = add&2 ? _acc : acc;

      // next get pot & acc for all or all flagged bodies
      if(Acc) {
	set_time(time);
	switch(ndim) {
	case 2:
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	    if(flag==0 || flag[n] & 1)
	      Acc->template acc<2>(mas+n, pos+nn, vel+nn, pots[n], accs+nn);
	  break;
	case 3:
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=3)
	    if(flag==0 || flag[n] & 1)
	      Acc->template acc<3>(mas+n, pos+nn, vel+nn, pots[n], accs+nn);
	  break;
	default:
	  error("acceleration \"%s\": ndim=%d not supported",
		Acceleration::name(),ndim);
	}
      } else {
	switch(ndim) {
	case 2:
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	    if(flag==0 || flag[n] & 1) {
	      pots[n   ] = scalar(0);
	      accs[nn+0] = scalar(0);
	      accs[nn+1] = scalar(0);
	    }
	  break;
	case 3:
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=3)
	    if(flag==0 || flag[n] & 1) {
	      pots[n   ] = scalar(0);
	      accs[nn+0] = scalar(0);
	      accs[nn+1] = scalar(0);
	      accs[nn+2] = scalar(0);
	    }
	  break;
	default:
	  error("acceleration \"%s\": ndim=%d not supported",
		Acceleration::name(),ndim);
	}
      }
      // in case of potential to be added, add it now & delete array
      if(_pot) {
	for(int n=0; n!=nbod; ++n)
	  if(flag==0 || flag[n] & 1)
	    pot[n] += _pot[n];
	delete[] _pot;
      }
 
      // in case of acceleration to be added, add it now & delete array
      if(_acc) {
	switch(ndim) {
	case 2:
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	    if(flag==0 || flag[n] & 1) {
	      acc[nn]   += _acc[nn];
	      acc[nn+1] += _acc[nn+1];
	    }
	  break;
	case 3:
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=3)
	    if(flag==0 || flag[n] & 1) {
	      acc[nn]   += _acc[nn];
	      acc[nn+1] += _acc[nn+1];
	      acc[nn+2] += _acc[nn+2];
	    }
	}
	delete[] _acc;
      }
    }

  public:
    AccInstall() : Acc(0), First(true), Pars(0), Npar(0) {}

   ~AccInstall()
    {
      if(Acc)  delete   Acc;
      if(Pars) delete[] Pars;
      if(File) delete[] File;
    }

    void init(const double*pars,
	      int          npar,
	      const char  *file,
	      bool        *need_mass,
	      bool        *need_vels)
    {
      if(differ(pars,npar,file)) {
	if(Acc) {
	  if(need_mass)
	    warning("iniacceleration: re-initializing \"%s\" "
		    "with different parameters or data file",
		    Acceleration::name());
	  else
	    warning("inipotential: re-initializing \"%s\" "
		    "with different parameters or data file",
		    Acceleration::name());
	  delete Acc;
	}
	Acc = new Acceleration(pars,npar,file);
	remember(pars,npar,file);
	First = true;
      } else if(Acc) {
	if(need_mass)
	  warning("iniacceleration: re-initializing \"%s\" "
		  "with idential parameters: no action taken",
		  Acceleration::name());
	else
	  warning("inipotential: re-initializing \"%s\" "
		  "with idential parameters: no action taken",
		  Acceleration::name());
      }
      if(need_mass)
	*need_mass = Acc->NeedMass();
      else if(Acc->NeedMass())
	error("inipotential: cannot use \"%s\", since masses are needed",
	      Acceleration::name());
      if(need_vels)
	*need_vels = Acc->NeedVels();
      else if(Acc->NeedVels())
	error("inipotential: cannot use \"%s\", since velocities are needed",
	      Acceleration::name());
    }

    void acc(int        nd,
	     double     t,
	     int        nb,
	     const void*m,
	     const void*x,
	     const void*v,
	     const int *f,
	     void      *p,
	     void      *a,
	     int        add,
	     char       type)
    {
      if(!Acc)
	warning("acceleration \"%s\": not initialized", Acceleration::name());
      switch(type) {
      case 'f': return acc_T(nd,t,nb,
			     static_cast<const float*>(m),
			     static_cast<const float*>(x),
			     static_cast<const float*>(v),
			     f,
			     static_cast<float*>(p),
			     static_cast<float*>(a),
			     add);
      case 'd': return acc_T(nd,t,nb,
			     static_cast<const double*>(m),
			     static_cast<const double*>(x),
			     static_cast<const double*>(v),
			     f,
			     static_cast<double*>(p),
			     static_cast<double*>(a),
			     add);
      default: error("acceleration \"%s\": unknown type ('%s')",
		     Acceleration::name(),type);
      }
    }
#ifdef POT_DEF
    template<typename scalar>
    void pot_T(const int   *ndim,
	       const scalar*pos,
	       scalar      *acc,
	       scalar      *pot,
	       const scalar*time)
    {
      if(!Acc)
	warning("potential \"%s\": not initialized",Acceleration::name());
      set_time(*time);
      switch(*ndim) {
      case 2: Acc->template acc<2>(static_cast<const scalar*>(0), pos, 
				   static_cast<const scalar*>(0), *pot, acc);
	break;
      case 3: Acc->template acc<3>(static_cast<const scalar*>(0), pos, 
				   static_cast<const scalar*>(0), *pot, acc);
	break;
      default: error("potential \"%s\": ndim=%d not supported",
		     Acceleration::name(),ndim);
      }
    }
#endif
  };
  
  template<int N, int I=0> struct Dims {
    typedef Dims<N,I+1> D;
    template<typename X>
    static X norm(const X*a) { return a[I]*a[I] + D::norm(a); }
    template<typename X>
    static void ass (X*a, const X*b) { a[I] =b[I]; D::ass(a,b); }
    template<typename X>
    static void add (X*a, const X*b) { a[I]+=b[I]; D::add(a,b); }
    template<typename X>
    static void asstimes(X*a,const X*b,X x) { a[I] =x*b[I]; D::asstimes(a,b,x);}
    template<typename X>
    static void addtimes(X*a,const X*b,X x) { a[I]+=x*b[I]; D::addtimes(a,b,x);}
  };
  template<int I> struct Dims<I,I> {
    template<typename X>
    static X norm(const X*a) { return a[I]*a[I]; }
    template<typename X>
    static void ass (X*a, const X*b) { a[I] =b[I]; }
    template<typename X>
    static void add (X*a, const X*b) { a[I]+=b[I]; }
    template<typename X>
    static void asstimes(X*a,const X*b,X x) { a[I] =x*b[I]; }
    template<typename X>
    static void addtimes(X*a,const X*b,X x) { a[I]+=x*b[I]; }
  };
  template<int N, typename X>
  inline X norm(const X*a) { return Dims<N-1>::norm(a); }
  template<int N, typename X>
  inline void ass(X*a, const X*b) { Dims<N-1>::ass(a,b); }
  template<int N, typename X>
  inline void add(X*a, const X*b) { Dims<N-1>::add(a,b); }
  template<int N, typename X>
  inline void asstimes(X*a, const X*b,X x) { Dims<N-1>::asstimes(a,b,x); }
  template<int N, typename X>
  inline void addtimes(X*a, const X*b,X x) { Dims<N-1>::addtimes(a,b,x); }
}

#define __DEF__ACC(NAME)			\
namespace {					\
  AccInstall<NAME> MyAcc;			\
}						\
void iniacceleration(const double*pars,		\
		     int          npar,		\
		     const char  *file,		\
		     bool        *need_m,	\
		     bool        *need_v)	\
{						\
  MyAcc.init(pars,npar,file,need_m,need_v);	\
}						\
void acceleration(int        nd,		\
		  double     t,			\
		  int        nb,		\
		  const void*m,			\
		  const void*x,			\
		  const void*v,			\
		  const int *f,			\
		  void      *p,			\
		  void      *a,			\
		  int        add,		\
		  char       type)		\
{						\
  MyAcc.acc(nd,t,nb,m,x,v,f,p,a,add,type);	\
}

#ifdef POT_DEF

# define __DEF__POT				\
void inipotential(const int   *npar,		\
		  const double*pars,		\
		  const char  *file)		\
{						\
  MyAcc.init(pars,*npar,file,0,0);		\
}						\
void potential_double(const int   *ndim,	\
		      const double*pos,		\
		      double      *acc,		\
		      double      *pot,		\
		      const double*time)	\
{						\
  MyAcc.pot_T(ndim,pos,acc,pot,time);		\
}						\
void potential_float(const int  *ndim,		\
		     const float*pos,		\
		     float      *acc,		\
		     float      *pot,		\
		     const float*time)		\
{						\
  MyAcc.pot_T(ndim,pos,acc,pot,time);		\
}

#else  // POT_DEF
# define __DEF__POT
#endif // POT_DEF

#endif // __cplusplus

#endif // included_defacc_h
