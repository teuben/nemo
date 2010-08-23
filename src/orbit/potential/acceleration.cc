// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///
/// \file  acceleration.cc
///
/// \brief implements acceleration.h
///
/// \author Walter Dehnen
/// \date   2004,2010
/// \note   acceleration.h is a C header file, while this implementation is C++
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004,2010 Walter Dehnen
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
//------------------------------------------------------------------------------
//
// version 0.0  17/06/2004  WD  It works!
// version 1.0  22/06/2004  WD  allow for up to 10 fallbacks, no last_...
// version 2.0  24/08/2004  WD  have iniacceleration() give the acc_pter
// version 2.1  26/08/2004  WD  check for 2nd use of same potential
//                              16 fallbacks
// version 3.0  26/08/2004  WD  allow for multiple accnames
// version 3.1  26/08/2004  WD  avoid array allocation in fallback
// version 3.2  26/08/2004  WD  check for 2nd use of same acceleration
// version 3.3  17/02/2010  WD  allow for nemo string bug
// version 3.4  23/08/2010  WD  empty accnames: return 0 (rather than trigger
//                              Segmentation fault)
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

  const int FbMax = 16;  // max number of fall backs
                         // changes here must be reflect in changes below!!!
  int       FbInd = 0;   // index of first free fallback
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
		    float         &pot,
		    float         *acc)
    {
      if(pf)
	pf(&ndim,pos,acc,&pot,&t_float);
      else if(pd) {
	double d__pot;
	double d__pos[NDIM];
	double d__acc[NDIM];
	v_ass<NDIM>(d__pos,pos);
	pd(&ndim,d__pos,d__acc,&d__pot,&t_double);
	pot = d__pot;
	v_ass<NDIM>(acc,d__acc);
      } else {
	pot = 0.f;
	v_set<NDIM>(acc,0.f);
      }
    }
  };

  template<int NDIM> struct grav<double,NDIM> {
    static void get(potproc_double pd,
		    potproc_float  pf,
		    const double  *pos,
		    double        &pot,
		    double        *acc)
    {
      if(pd)
	pd(&ndim,pos,acc,&pot,&t_double);
      else if(pf) {
	float f__pot;
	float f__pos[NDIM];
	float f__acc[NDIM];
	v_ass<NDIM>(f__pos,pos);
	pf(&ndim,f__pos,f__acc,&f__pot,&t_float);
	pot = f__pot;
	v_ass<NDIM>(acc,f__acc);
      } else {
	pot = 0.;
	v_set<NDIM>(acc,0.);
      }
    }
  };

  //----------------------------------------------------------------------------
  // declare type of pointer to inipotential()
  typedef void(*inipot_pter)       // return: void
    (int*,                         // input:  number of parameters
     double*,                      // input:  array with parameters
     char*);                       // input:  data file name (may be NULL)
  //----------------------------------------------------------------------------
  // class fallback
  class fallback {
    inipot_pter    ip;
    potproc_double pd;
    potproc_float  pf;
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void acc_TT(int          nbod,
		const scalar*pos,
		const int   *flag,
		scalar      *pot,
		scalar      *acc,
		int          add)
    {
      if(add & 1) {
	if(add & 2) {
	  // add both potential and acceleration
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	    if(flag==0 || flag[n] & 1) {
	      register scalar P, A[NDIM];
	      grav<scalar,NDIM>::get(pd,pf,pos+nn, P, A);
	      pot[n] += P;
	      v_add<NDIM>(acc+nn,A);
	    }
	} else {
	  // add potential, assign acceleration
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	    if(flag==0 || flag[n] & 1) {
	      register scalar P;
	      grav<scalar,NDIM>::get(pd,pf,pos+nn, P, acc+nn);
	      pot[n] += P;
	    }
	}
      } else {
	if(add & 2) {
	  // assign potential, add acceleration
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	    if(flag==0 || flag[n] & 1) {
	      register scalar A[NDIM];
	      grav<scalar,NDIM>::get(pd,pf,pos+nn, pot[n], A);
	      v_add<NDIM>(acc+nn,A);
	    }
	} else {
	  // assign both potential and acceleration
	  for(int n=0,nn=0; n!=nbod; ++n,nn+=NDIM)
	    if(flag==0 || flag[n] & 1)
	      grav<scalar,NDIM>::get(pd,pf,pos+nn, pot[n], acc+nn);
	}
      }
    }
    //--------------------------------------------------------------------------

    template<typename scalar>
    void acc_T(int          nbod,
	       const scalar*pos,
	       const int   *flag,
	       scalar      *pot,
	       scalar      *acc,
	       int          add)
    {
      switch(ndim) {
      case 2: acc_TT<2,scalar>(nbod,pos,flag,pot,acc,add); break;
      case 3: acc_TT<3,scalar>(nbod,pos,flag,pot,acc,add); break;
      default: error(const_cast<char*>("acceleration: unsupported ndim (%d)"),ndim);
      }
    }
    
  public:

    fallback() : ip(0), pd(0), pf(0) {}

    void set (inipot_pter i, potproc_double d, potproc_float f)
    {
      ip = i;
      pd = d;
      pf = f;
    }

    inipot_pter const&IniPotPter() const { return ip; }

    bool is_set() const
    {
      return pd!=0 || pf!=0;
    }

    void accel(int        nd,
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
	warning(const_cast<char*>("fallback::acceleration: potential not set"));
      
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
      else error(const_cast<char*>("fallback::acceleration: unknown type ('%s')"),type);
    }

  } 
  // define array with FbMax FallBacks
  FallBack[FbMax];

  // define the functions AccFallBack0() to AccFallBack15() which employ
  // FallBack[i]. This means that we have function pointers for up to 
  // FbMax fall backs
#define ACC_FALLBACK(NUM)					\
  void AccFallBack##NUM(int        d,				\
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
  { (FallBack[NUM]).accel(d,t,n,m,x,v,f,p,a,i,y); }

  ACC_FALLBACK(0)
  ACC_FALLBACK(1)
  ACC_FALLBACK(2)
  ACC_FALLBACK(3)
  ACC_FALLBACK(4)
  ACC_FALLBACK(5)
  ACC_FALLBACK(6)
  ACC_FALLBACK(7)
  ACC_FALLBACK(8)
  ACC_FALLBACK(9)
  ACC_FALLBACK(10)
  ACC_FALLBACK(11)
  ACC_FALLBACK(12)
  ACC_FALLBACK(13)
  ACC_FALLBACK(14)
  ACC_FALLBACK(15)

#undef ACC_FALLBACK

  // define array with function pointers
  acc_pter AccFallBack[FbMax] = {&AccFallBack0,
				 &AccFallBack1,
				 &AccFallBack2,
				 &AccFallBack3,
				 &AccFallBack4,
				 &AccFallBack5,
				 &AccFallBack6,
				 &AccFallBack7,
				 &AccFallBack8,
				 &AccFallBack9,
				 &AccFallBack10,
				 &AccFallBack11,
				 &AccFallBack12,
				 &AccFallBack13,
				 &AccFallBack14,
				 &AccFallBack15};
  // table of accnames used for the fallbacks sofar
  char acc_names[FbMax][256] = {0};

  //                                                                            
  // 2 single_acceleration()                                                    
  //                                                                            
  // 0. we try to match accname against those already successfully done;        
  //    if match found, we simply call the corresponding iniacceleration() and  
  //    return.                                                                 
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
  //----------------------------------------------------------------------------
  const int IniAcMax=256;
  const int AcNamMax=128;
  int       IniAcInd=0;
  char AcNames[IniAcMax][AcNamMax];
  iniacc_pter IniAc[IniAcMax] = {0};
  //----------------------------------------------------------------------------
  // this routine used to be get_acceleration()
  // we now allow for several accelerations added, the names, parameters, and
  // files of which are seperated by ',', ';', and ';', respectively.
  acc_pter single_acceleration(const char*accname,
			       const char*accpars,
			       const char*accfile,
			       bool      *need_mass,
			       bool      *need_vels)
  {
    //
    // NOTE: the present implementation will NOT try to compile a source code
    //       but abort if no .so file is found.
    //

    nemo_dprintf(2,"single_acceleration(\"%s\", \"%s\", \"%s\")\n",
		 accname,accpars,accfile);
    // 1. parse the parameters (if any)
    const int MAXPAR = 64;
    double pars[MAXPAR];
    int    npar;
    if(accpars && *accpars) {
      npar = nemoinpd(const_cast<char*>(accpars),pars,MAXPAR);
      if(npar > MAXPAR)
	error(const_cast<char*>("get_acceleration: too many parameters (%d > %d)"),npar,MAXPAR);
      if(npar < 0)
	warning(const_cast<char*>("get_acceleration: parsing error in parameters: \"%s\""),
		accpars);
    } else
      npar = 0;
  
    // 2. load local symbols
    if(first) {
      mysymbols(getparam(const_cast<char*>("argv0")));
      first = false;
    }

    // 3. try to find accname in list of accnames already done
    for(int i=0; i!=IniAcInd; ++i)
      if(0 == strcmp(accname, AcNames[i])) {
	nemo_dprintf(2,"single_acceleration: accname=\"%s\": "
		     "use iniacc_pter known already\n",accname);
	acc_pter ac;
	(*IniAc[i])(pars,npar,accfile,&ac,need_mass,need_vels);
	return ac;
      }

    // 4. seek file accname.so and load it
    char path[256] = ".";
    char*accpath = getenv("ACCPATH");             // try $ACCPATH
    if(accpath == 0) accpath = getenv("POTPATH"); // try $POTPATH
    if(accpath == 0) {
      accpath = path;                             // try .
      char* nemopath = getenv("NEMO");            // IF $NEMO defined
      if(nemopath) {
	strcat(accpath,":");
	strcat(accpath,nemopath);
	strcat(accpath,"/obj/potential");         //   try .:$NEMO/obj/potential
      }
    }      
    char name[256];
    strcpy(name,accname);
    strcat(name,".so");
    char*fullname = pathfind(accpath,name);       // seek for file in path
    if(fullname == 0)
      error(const_cast<char*>("get_acceleration: cannot find file \"%s\" in path \"%s\""),
	    name,accpath);
    loadobj(fullname);

    // 5. try to get iniacceleration()
    //    if found, remember it, call it and return acceleration() given by it
    iniacc_pter ia = (iniacc_pter) findfunc("iniacceleration");
    if(ia) {
      if(IniAcInd < IniAcMax && strlen(accname) < AcNamMax) {
	strcpy(AcNames[IniAcInd],accname);
	IniAc[IniAcInd] = ia;
	IniAcInd++;
      }
      acc_pter ac;
      ia(pars,npar,accfile,&ac,need_mass,need_vels);
      return ac;
    }
    warning(const_cast<char*>("get_acceleration: no acceleration found in file \"%s\", "
			      "trying potential instead"),fullname);

    // 6. fall-back: use potential.
    //    NOTE: 1 we cannot just call get_potential(), for we want BOTH
    //            potential_float and potential_double, but must call
    //            inipotential() only once
    //          2 only one instance of any external potentials can exist
    //            (there is only one function potential() to load), so if
    //            the same accname occurs a second time, we simply warn
    //            and call inipotential() again (no loading needed).

    if(need_mass) *need_mass = false;
    if(need_vels) *need_vels = false;
    // 6.0 check for second call of some accname
    int fb = 0;
    for(; fb!=FbInd; ++fb)
      if(0 == strcmp(accname, acc_names[fb])) break;
    if(fb < FbInd) {
      // 6.A second call for accname:
      // 6.A.1 warn about re-initialisation
      warning(const_cast<char*>("get_acceleration: re-initializing potential \"%s\""),accname);
      // 6.A.2 call inipotential() again
      (FallBack[fb]).IniPotPter()(&npar,pars,const_cast<char*>(accfile));
    } else {
      // 6.B first call for accname:
      // 6.B.0 remember accname, error
      if(fb >= FbMax)
	error(const_cast<char*>("get_acceleration: cannot fallback more than %d times"),FbMax);
      strcpy(acc_names[fb],accname);
      ++FbInd;
      // 6.B.1 try to get inipotential and potential
      inipot_pter    ip = (inipot_pter)    findfunc("inipotential");
      potproc_double pd = (potproc_double) findfunc("potential_double");
      if(pd == 0)    pd = (potproc_double) findfunc("potential");
      potproc_float  pf = (potproc_float)  findfunc("potential_float");
      if((pf==0 && pd==0) || ip==0)
	error(const_cast<char*>("get_acceleration: no potential found either"));
      // 6.B.2 call inipotential() once
      ip(&npar,pars,const_cast<char*>(accfile));
      // 6.B.3 initialize FallBack[fb] and return AccFallBack[fb]
      (FallBack[fb]).set(ip,pd,pf);
    }
    return AccFallBack[fb];
  }
  //                                                                            
  // 3 added accelerations                                                      
  //                                                                            
  const int AcMax = 10;                      // max number of added accs        
  int       AcInd = 0;                       // next free added acc             
  //----------------------------------------------------------------------------
  class __AddedAccs {
    static const int NMAX=10;                // max number of potentials        
    int              N;                      // actual number of potentials     
    acc_pter         AC[NMAX];               // pointers to accelerations       
    //--------------------------------------------------------------------------
  public:
    void set(int         nacc,
	     const char**accname,
	     const char**accpars,
	     const char**accfile,
	     bool       *need_mass,
	     bool       *need_vels)
    {
      if(nacc > NMAX)
	error(const_cast<char*>("get_acceleration: more accnames than expected (%d)"),NMAX);
      N = nacc;
      if(need_mass) *need_mass = 0;
      if(need_vels) *need_vels = 0;
      for(int n=0; n!=N; ++n) {
	if(accname[n]==0 || accname[n][0] == 0)
	  error(const_cast<char*>("get_acceleration: accname #%d empty "
				  "(parse error in \"accname=...\"?)"),n);
	bool nm,nv;
	AC[n] = single_acceleration(accname[n],
				    (accpars[n] && accpars[n][0])? accpars[n]:0,
				    (accfile[n] && accfile[n][0])? accfile[n]:0,
				    &nm,&nv);
	if(need_mass && nm) *need_mass = 1;
	if(need_vels && nv) *need_vels = 1;
      }
    }
    //--------------------------------------------------------------------------
    void accel(int        ndim,              // I: number of dimensions         
	       double     time,              // I: simulation time              
	       int        nbod,              // I: number bodies =size of arrays
	       const void*mass,              // I: masses:         m[i]         
	       const void*pos,               // I: positions       (x,y,z)[i]   
	       const void*vel,               // I: velocities      (u,v,w)[i]   
	       const int *flag,              // I: flags           f[i]         
	       void      *pot,               // O: potentials      p[i]         
	       void      *acc,               // O: accelerations   (ax,ay,az)[i]
	       int        add,               // I: indicator (see note 6 above) 
	       char       type)              // I: type: 'f' or 'd'             
    {
      for(int n=0; n<N; ++n)
	(*(AC[n]))(ndim,time,nbod,mass,pos,vel,flag,pot,acc,n?3:add,type);
    }

  } // class __AddedAccs
  // define array with AcMax __AddedAccs
  Added[AcMax];

  // define the functions AddedAcc0() to AddedAcc9() which employ
  // Added[i]. This means that we have function pointers for up to AcMax
  // sums of acceleration fields
#define ACC_ADDED(NUM)					\
  void AddedAcc##NUM(int        d,			\
	  	     double     t,			\
		     int        n,			\
		     const void*m,			\
		     const void*x,			\
		     const void*v,			\
		     const int *f,			\
	             void      *p,			\
	             void      *a,			\
	             int        i,			\
	             char       y)			\
  { (Added[NUM]).accel(d,t,n,m,x,v,f,p,a,i,y); }

  ACC_ADDED(0)
  ACC_ADDED(1)
  ACC_ADDED(2)
  ACC_ADDED(3)
  ACC_ADDED(4)
  ACC_ADDED(5)
  ACC_ADDED(6)
  ACC_ADDED(7)
  ACC_ADDED(8)
  ACC_ADDED(9)

#undef ACC_ADDED

  // define array with function pointers
  acc_pter AddedAcc[AcMax] = {&AddedAcc0,
			      &AddedAcc1,
			      &AddedAcc2,
			      &AddedAcc3,
			      &AddedAcc4,
			      &AddedAcc5,
			      &AddedAcc6,
			      &AddedAcc7,
			      &AddedAcc8,
			      &AddedAcc9};

  //----------------------------------------------------------------------------
  typedef const char* c_string;

  inline bool is_sep(char const&c, const char*seps) {
    for(const char*s=seps; *s; ++s)
      if( c == *s ) return true;
    return false;
  }

  c_string *splitstring(c_string list, c_string seps)
  {
    char    *words = new char[strlen(list)+1];
    c_string*wlist = new c_string[128];
    int   n    = 0;
    char *w    = words;
    c_string l = list;
    wlist[n]   = w;
    do {                                      // LOOP list
      if(*l == 0 || is_sep(*l,seps)) {        //   IF end or seperator found
	*w++ = 0;                             //     close word
	if(++n == 128) error(const_cast<char*>("too many words in list"));
	wlist[n] = w;                         //     get begin of new word
      } else                                  //   ELSE
	*w++ = *l;                            //     copy character
    } while (*l++ != 0);                      // until list ends
    wlist[n] = 0;                             // close word list
//     for(int i=0; i!=n; ++i)
//       nemo_dprintf(2,"wlist[%d]=%d:\"%s\"\n",i,int(wlist[i]-words),wlist[i]);
    return wlist;
  }

  inline void freestrings(c_string*wlist) {
    delete[] wlist[0];
    delete[] wlist;
  }
  //////////////////////////////////////////////////////////////////////////////
  char NameSeps[3] = {',','+',0};
  char ParsSeps[3] = {';','#',0};
  char FileSeps[3] = {';','#',0};

} // namespace {
////////////////////////////////////////////////////////////////////////////////
//
// finally, after all this, we can define get_acceleration() itself!
//
////////////////////////////////////////////////////////////////////////////////
acc_pter get_acceleration(const char*accnames,
			  const char*accparss,
			  const char*accfiles,
			  bool      *need_mass,
			  bool      *need_vels)
{
  // 0. split accnames and count number of sub acceleration fields
  nemo_dprintf(2,"get_acceleration(\"%s\",\"%s\",\"%s\")\n",
	       accnames,accparss,accfiles);
  
  int nacc=0;
  c_string*accname;
  if(accnames) {
    accname = splitstring(accnames,NameSeps);
    while(accname[nacc]) nacc++;
  }
  if(nacc == 0) return 0;
  // 1. just one accname, then return single acceleration()
  if(nacc == 1) {
    freestrings(accname);
    return single_acceleration(accnames, accparss, accfiles,
			       need_mass, need_vels);
  }
  // 2. several accnames
  // 2.1 split accparss, allow for empty accparss -> all accpars=0
  c_string accpars_[1024] = {0}, *accpars = accpars_;
  if(accparss) {
    accpars = splitstring(accparss,ParsSeps);
    int n=0; while(accpars[n]) n++;
    if(n!=nacc)
      error(const_cast<char*>("get_acceleration: %d names but %d parameter sets"),nacc,n);
  }
  // 2.2 split accfiles, allow for empty accfiles -> all accfile=0
  c_string accfile_[1024] = {0}, *accfile = accfile_;
  if(accfiles) {
    accfile = splitstring(accfiles,FileSeps);
    int n=0; while(accfile[n]) n++;
    if(n!=nacc)
      error(const_cast<char*>("get_acceleration: %d names but %d data files"),nacc,n);
  }
  // 2.3 check for availability of another added acceleration
  if(AcInd >= AcMax)
    error(const_cast<char*>("get_acceleration: called more than %d times with multiple accnames"),
	  AcMax);
  // 2.4 initialize added acceleration and return function using it
  (Added[AcInd]).set(nacc,accname,accpars,accfile,need_mass,need_vels);
  freestrings(accname);
  if(accparss) freestrings(accpars);
  if(accfiles) freestrings(accfile);
  return AddedAcc[AcInd++];
}
////////////////////////////////////////////////////////////////////////////////
