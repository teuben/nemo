//------------------------------------------------------------------------------
//                                                                              
// acceleration.cc                                                              
//                                                                              
// implements acceleration.h                                                    
//                                                                              
// Note that acceleration.h is a C header file, while this implementation       
// is in C++.                                                                   
//------------------------------------------------------------------------------
//                                                                              
// version 0.0  17/06/2004  WD   It works!                                      
// version 1.0  22/06/2004  WD   allow for up to 10 fallbacks, no last_...      
// version 2.0  24/08/2004  WD   have iniacceleration() give the acc_pter       
//                                                                              
//------------------------------------------------------------------------------

extern "C" {
#include  <acceleration.h>      // the header we are implementing
}
#include  <cstring>             // C type string manipultions
extern "C" {
#include  <stdinc.h>            // nemo's string (used in potential.h)
#include  <potential.h>         // the external potential for falling back
#include  <getparam.h>          // getting name of main()
#include  <loadobj.h>           // loading shared object files
#include  <filefn.h>            // finding a function in a loaded file
}

//==============================================================================
//                                                                              
// 1 Fall-back mechanism.                                                       
//                                                                              
// We shall implement a fall-back mechanism to catch cases for which            
// there is an inipotential() and potential(), but no iniacceleration()         
// yet. In this case, we return a simple acceleration() which just calls the    
// functions defined in potential.h.                                            
//                                                                              
// implemented via templating the data type.                                    
// That's why this is implemented in C++ rather than C.                         
//                                                                              

namespace {

  //
  // auxiliary template classes used by the helper functions below
  //
  template<int N, int I=0> struct Dims {
    typedef Dims<N,I+1> D;
    template<typename X>
    static X norm(const X*a) { return a[I]*a[I] + D::norm(a); }
    template<typename X>
    static void set (X*a, X x) { a[I] =x; D::set(a,x); }
    template<typename X, typename Y>
    static void ass (X*a, const Y*b) { a[I] =b[I]; D::ass(a,b); }
    template<typename X, typename Y>
    static void add (X*a, const Y*b) { a[I]+=b[I]; D::add(a,b); }
    template<typename X, typename Y>
    static void asstimes(X*a,const Y*b,Y x) { a[I] =x*b[I]; D::asstimes(a,b,x);}
    template<typename X, typename Y>
    static void addtimes(X*a,const Y*b,Y x) { a[I]+=x*b[I]; D::addtimes(a,b,x);}
  };
  template<int I> struct Dims<I,I> {
    template<typename X>
    static X norm(const X*a) { return a[I]*a[I]; }
    template<typename X>
    static void set (X*a, X x) { a[I] =x; }
    template<typename X, typename Y>
    static void ass (X*a, const Y*b) { a[I] =b[I]; }
    template<typename X, typename Y>
    static void add (X*a, const Y*b) { a[I]+=b[I]; }
    template<typename X, typename Y>
    static void asstimes(X*a,const Y*b,Y x) { a[I] =x*b[I]; }
    template<typename X, typename Y>
    static void addtimes(X*a,const Y*b,Y x) { a[I]+=x*b[I]; }
  };
  //
  // helper functions useful for manipulating pointers to scalar
  // as N dimensional vectors.
  //
  template<int N, typename X>
  inline X v_norm(const X*a) { return Dims<N-1>::norm(a); }
  template<int N, typename X>
  inline void v_set(X*a, X x) { Dims<N-1>::set(a,x); }
  template<int N, typename X, typename Y>
  inline void v_ass(X*a, const Y*b) { Dims<N-1>::ass(a,b); }
  template<int N, typename X, typename Y>
  inline void v_add(X*a, const Y*b) { Dims<N-1>::add(a,b); }
  template<int N, typename X, typename Y>
  inline void v_asstimes(X*a, const Y*b,Y x) { Dims<N-1>::asstimes(a,b,x); }
  template<int N, typename X, typename Y>
  inline void v_addtimes(X*a, const Y*b,Y x) { Dims<N-1>::addtimes(a,b,x); }

  //
  // global variables associated with fallback mechanism
  //

  const int MaxFB = 10;  // max number of fall backs
  int       ndim;        // number of dimensions, used in call to potential()
  double    t_double;    // simulation time, used in call to potential()
  float     t_float;     // simulation time, used in call to potential()

  //----------------------------------------------------------------------------
  // grav<scalar>::get() auxiliary for defining fallback::acc_T<>()
  template<typename scalar, int NDIM> class grav;

  template<int NDIM> struct grav<float,NDIM> {
    static void get(potproc_double pd,
		    potproc_float  pf,
		    const float   *pos,
		    float         *pot,
		    float         *acc)
    {
      if(pf)
	pf(&ndim,pos,acc,pot,&t_float);
      else if(pd) {
	double d__pot;
	double d__pos[NDIM];
	double d__acc[NDIM];
	v_ass<NDIM>(d__pos,pos);
	pd(&ndim,d__pos,d__acc,&d__pot,&t_double);
	pot[0] = d__pot;
	v_ass<NDIM>(acc,d__acc);
      } else {
	pot[0] = 0.f;
	v_set<NDIM>(acc,0.f);
      }
    }
  };

  template<int NDIM> struct grav<double,NDIM> {
    static void get(potproc_double pd,
		    potproc_float  pf,
		    const double  *pos,
		    double        *pot,
		    double        *acc)
    {
      if(pd)
	pd(&ndim,pos,acc,pot,&t_double);
      else if(pf) {
	float f__pot;
	float f__pos[NDIM];
	float f__acc[NDIM];
	v_ass<NDIM>(f__pos,pos);
	pf(&ndim,f__pos,f__acc,&f__pot,&t_float);
	pot[0] = f__pot;
	v_ass<NDIM>(acc,f__acc);
      } else {
	pot[0] = 0.;
	v_set<NDIM>(acc,0.);
      }
    }
  };

  //----------------------------------------------------------------------------
  // class fallback
  class fallback {
    potproc_double pd;
    potproc_float  pf;

    template<typename scalar>
    void acc_T(int          nbod,
	       const scalar*pos,
	       const int   *flag,
	       scalar      *pot,
	       scalar      *acc,
	       int          add)
    {

      // in case of adding, create arrays to write pot & acc into
      scalar*_pot = add&1 ? new scalar[nbod]      : 0;
      scalar*_acc = add&2 ? new scalar[nbod*ndim] : 0;

      // define references to the arrays actually passed to grav::get()
      scalar*&pots = add&1 ? _pot : pot;
      scalar*&accs = add&2 ? _acc : acc;

      // next get pot & acc for all or all flagged bodies
      if(ndim==2) {
	for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	  if(flag==0 || flag[n] & 1)
	    grav<scalar,2>::get(pd,pf,pos+nn, pots+n, accs+nn);
      } else {
	for(int n=0,nn=0; n!=nbod; ++n,nn+=3)
	  if(flag==0 || flag[n] & 1)
	    grav<scalar,3>::get(pd,pf,pos+nn, pots+n, accs+nn);
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
	if(ndim==2) {
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=2)
	    if(flag==0 || flag[n] & 1)
	      v_add<2>(acc+nn,_acc+nn);
	} else {
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=3)
	    if(flag==0 || flag[n] & 1)
	      v_add<3>(acc+nn,_acc+nn);
	}
	delete[] _acc;
      }
    }
    
  public:

    fallback() {}

    void set (potproc_double d, potproc_float f)
    {
      pd = d;
      pf = f;
    }

    bool is_set() const
    {
      return pd!=0 || pf!=0;
    }

    void acceleration(int        nd,
		      double     t,
		      int        n,
		      const void*m,
		      const void*x,
		      const void*v,
		      const int *f,
		      void      *p,
		      void      *a,
		      int        add,
		      char       type)
    {
      if(!is_set())
	warning("fallback::acceleration: potential not set");
      
      // set ndim and time
      ndim     = nd;
      t_float  = t;
      t_double = t;

      // now call acceleration_t ()
      if     (type=='f')
	acc_T(n, static_cast<const float*>(x),
	      f, static_cast<float*>(p), static_cast<float*>(a), add);
      else if(type=='d')
	acc_T(n, static_cast<const double*>(x),
	      f, static_cast<double*>(p), static_cast<double*>(a), add);
      else error("fallback::acceleration: unknown type ('%s')",type);
    }

  } 
  // define array with MaxFB=10 FallBacks
  FallBack[MaxFB];

  // define the acceleration0() to acceleration9() which employ
  // FallBack[i]. This means that we have function pointers for
  // up to 10 fall backs
#define ACC_FALLBACK(NUM)					\
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
{ (FallBack[NUM]).acceleration(d,t,n,m,x,v,f,p,a,i,y); }

  ACC_FALLBACK(0);
  ACC_FALLBACK(1);
  ACC_FALLBACK(2);
  ACC_FALLBACK(3);
  ACC_FALLBACK(4);
  ACC_FALLBACK(5);
  ACC_FALLBACK(6);
  ACC_FALLBACK(7);
  ACC_FALLBACK(8);
  ACC_FALLBACK(9);

#undef ACC_FALLBACK

  // define array with function pointers
  acc_pter AccFallBack[MaxFB] = {&acceleration0,
				 &acceleration1,
				 &acceleration2,
				 &acceleration3,
				 &acceleration4,
				 &acceleration5,
				 &acceleration6,
				 &acceleration7,
				 &acceleration8,
				 &acceleration9};
  // index of next fallback to use; must not be >= MaxFB
  int fb_index = 0;

} // anonymous namespace
// end of fall-back mechanism

//==============================================================================
//                                                                              
// 2 get_acceleration()                                                         
//                                                                              
// 1. we try to load the routines iniacceleration() and acceleration() and,     
//    if successful, call the first and return the latter.                      
//                                                                              
// 2. Failing that, we try to load the routines inipotential(),                 
//    potential_float(), and potential_double() and restore to the              
//    fall-back mechanism.                                                      
//    NOTE: we cannot just call get_potential(), for we want BOTH               
//          potential_float and potential_double, but must call                 
//          inipotential() only once                                            
//                                                                              
//==============================================================================


// declare more auxiliary stuff
namespace {

  // declare type of pointer to iniacceleration()
  typedef void(*iniacc_pter)       // return: void
    (const double*,                // input:  array with parameters
     int,                          // input:  number of parameters
     const char*,                  // input:  data file name (may be NULL)
     acc_pter*,                    // output: pter to acc_pter to acceleration()
     bool*,                        // output: acceleration() needs masses?
     bool*);                       // output: acceleration() needs vel's?

  // declare type of pointer to inipotential()
  typedef void(*inipot_pter)       // return: void
    (int*,                         // input:  number of parameters
     double*,                      // input:  array with parameters
     char*);                       // input:  data file name (may be NULL)

  // boolean indicating whether we have loaded local symbols
  bool first = true;

  typedef void(*proc)();
  // routine for finding a function  "NAME", "NAME_" or "NAME__"
  inline proc findfunc(const char*func)
  {
    char fname[256];
    strcpy(fname,func);
    mapsys(fname);
    proc f = findfn(fname);
    for(int i=0; f==0 && i!=2; ++i) {
      strcat(fname,"_");
      f= findfn(fname);
    }
    return f;
  }

} // anonymous namespace

//
// finally define get_acceleration()
//
acc_pter get_acceleration(
			  const char*accname,
			  const char*accpars,
			  const char*accfile,
			  bool      *need_mass,
			  bool      *need_vels)
{
  //
  // NOTE: the present implementation will NOT try to compile a source code
  //       but abort if no .so file is found.
  //

  // 1. parse the parameters (if any)
  const int MAXPAR = 64;
  double pars[MAXPAR];
  int    npar;
  if(accpars && *accpars) {
    npar = nemoinpd(const_cast<char*>(accpars),pars,MAXPAR);
    if(npar > MAXPAR)
      error("get_acceleration: too many parameters (%d > %d)",npar,MAXPAR);
    if(npar < 0)
      warning("get_acceleration: parsing error in parameters: \"%s\"",accpars);
  } else
    npar = 0;
  
  // 2. load local symbols
  if(first) {
    mysymbols(getparam("argv0"));
    first = false;
  }

  // 3. seek file accname.so and load it
  char path[256] = ".";
  char*accpath = getenv("ACCPATH");               // try $ACCPATH
  if(accpath == 0) accpath = getenv("POTPATH");   // try $POTPATH
  if(accpath == 0) {
    accpath = path;                               // try .
    char* nemopath = getenv("NEMO");              // IF $NEMO defined
    if(nemopath) {
      strcat(accpath,":");
      strcat(accpath,nemopath);
      strcat(accpath,"/obj/potential");           //   try .:$NEMO/obj/potential
    }
  }      
  char name[256];
  strcpy(name,accname);
  strcat(name,".so");
  char*fullname = pathfind(accpath,name);         // seek for file in path
  if(fullname == 0)
    error("get_acceleration: cannot find file \"%s\" in path \"%s\"",
	  name,accpath);
  loadobj(fullname);

  // 4. try to get iniacceleration()
  //    if found, call it and return acceleration() given by it
  iniacc_pter ia = (iniacc_pter) findfunc("iniacceleration");
  if(ia) {
    acc_pter ac;
    ia(pars,npar,accfile,&ac,need_mass,need_vels);
    return ac;
  }
  warning("get_acceleration: no acceleration found in file \"%s\", "
	  "trying potential instead",fullname);

  // 5. fall-back: use potential.
  //    NOTE: we cannot just call get_potential(), for we want BOTH
  //          potential_float and potential_double, but must call
  //          inipotential() only once

  if(fb_index >= MaxFB)
    error("get_acceleration: cannot fallback more than %d times",MaxFB);

  // 5.1 try to get inipotential and potential
  inipot_pter    ip = (inipot_pter)    findfunc("inipotential");
  potproc_double pd = (potproc_double) findfunc("potential_double");
  if(pd == 0)    pd = (potproc_double) findfunc("potential");
  potproc_float  pf = (potproc_float)  findfunc("potential_float");
  if((pf==0 && pd==0) || ip==0)
    error("get_acceleration: no potential found either");

  // 5.2 call inipotential() once
  ip(&npar,pars,const_cast<char*>(accfile));
  *need_mass = false;
  *need_vels = false;

  // 5.3 initialize FallBack[fb_index] and return AccFallBack[fb_index++]
  (FallBack[fb_index]).set(pd,pf);
  return AccFallBack[fb_index++];
}

void dummy() {
  ::warning     (  "fool the linker to include 'warning'");
  ::error       (  "fool the linker to include 'error'");
  ::nemo_dprintf(2,"fool the linker to include 'nemo_dprintf'");
}
 
////////////////////////////////////////////////////////////////////////////////
