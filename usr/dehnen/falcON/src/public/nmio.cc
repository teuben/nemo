//-----------------------------------------------------------------------------+
//                                                                             |
// nmio.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/nmio.h>
#ifdef   falcON_NEMO                               // empty if NEMO not running 
#include <public/exit.h>
#include <iostream>
#include <cstdlib>
#include <cstring>

// nemo header files
extern "C" {
# include <stdinc.h>
# include <filestruct.h>
# include <history.h>
# include <snapshot/snapshot.h>
# ifndef GasDensTag
#   define GasDensTag "GasDensity"
# endif
# ifndef NumberTag
#   define NumberTag "NPartners"
# endif
# ifndef SPHNumberTag
#   define SPHNumberTag "NSPHPartners"
# endif
}
#ifndef DensityTag
#  error
#  error DensityTag not #defined
#  error    you are presumably and old version of NEMO
#  error    please update to a more recent version
#  error
#endif
#ifndef FlagTag
#  define FlagTag "Flag"
#endif
#ifndef LevelTag
#  define LevelTag "Level"
#endif

#ifdef falcON_REAL_IS_FLOAT
#  define FalconType FloatType
#else
#  define FalconType DoubleType
#endif


#define mystream static_cast<stream>(STREAM)

using std::cerr;
using namespace nbdy;
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::Set X)
{
  switch(X) {
  case nemo_io::snap:   return SnapShotTag;
  case nemo_io::param:  return ParametersTag;
  case nemo_io::bodies: return ParticlesTag;
  case nemo_io::diags:  return DiagnosticsTag;
  default: nbdy::error("Unknown nemo_io::Set");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline int tag(const nemo_io::CoSys X)
{
  switch(X) {
  case nemo_io::cart:  return Cartesian;
  case nemo_io::spher: return Spherical;
  case nemo_io::scat:  return Scattering;
  default: nbdy::error("Unknown nemo_io::CoSys");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::SingleScalar X)
{
  switch(X) {
  case nemo_io::time:    return TimeTag;
  case nemo_io::cputime: return CPUTimeTag;
  default: nbdy::error("Unknown nemo_io::SingleScalar");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::SingleVector X)
{
  switch(X) {
  case nemo_io::energy:  return EnergyTag;
  default: nbdy::error("Unknown nemo_io::SingleVector");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::SingleMatrix X)
{
  switch(X) {
  case nemo_io::KinT: return KETensorTag;
  case nemo_io::PotT: return PETensorTag;
  case nemo_io::AmT:  return AMTensorTag;
  default: nbdy::error("Unknown nemo_io::SingleMatrix");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::SinglePhases X)
{
  switch(X) {
  case nemo_io::cofm: return CMPhaseSpaceTag;
  default: nbdy::error("Unknown nemo_io::SinglePhases");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::BodiesScalar X)
{
  switch(X) {
  case nemo_io::mass: return MassTag;
  case nemo_io::pot:  return PotentialTag;
  case nemo_io::rho:  return DensityTag;
  case nemo_io::aux:  return AuxTag;
  case nemo_io::eps:  return EpsTag;
  default: nbdy::error("Unknown nemo_io::BodiesScalar");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::SPHScalar X)
{
  switch(X) {
  case nemo_io::uin:  return UinternTag;
  case nemo_io::udin: return UdotIntTag;
  case nemo_io::udex: return UdotRadTag;
  case nemo_io::entr: return EntFuncTag;
  case nemo_io::srho: return GasDensTag;
  case nemo_io::h   : return SmoothTag;
  default: nbdy::error("Unknown nemo_io::SPHScalar");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::BodiesVector X)
{
  switch(X) {
  case nemo_io::acc:    return AccelerationTag;
  case nemo_io::pos:    return PosTag;
  case nemo_io::vel:    return VelTag;
  default: nbdy::error("Unknown nemo_io::BodiesVector");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::BodiesPhases X)
{
  switch(X) {
  case nemo_io::posvel: return PhaseSpaceTag;
  default: nbdy::error("Unknown nemo_io::BodiesPhases");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::BodiesInteger X)
{
  switch(X) {
  case nemo_io::key:     return KeyTag;
  case nemo_io::flag:    return FlagTag;
  case nemo_io::numb:    return NumberTag;
  case nemo_io::numbSPH: return SPHNumberTag;
  default: nbdy::error("Unknown nemo_io::BodiesInteger");
  }
  return 0;
}
//------------------------------------------------------------------------------
inline char* tag(const nemo_io::BodiesShort X)
{
  switch(X) {
  case nemo_io::level:  return LevelTag;
  default: nbdy::error("Unknown nemo_io::BodiesShort");
  }
  return 0;
}
//------------------------------------------------------------------------------
void nemo_io::open(const char* file, const char* control)
{
  if(STREAM) strclose(mystream);
  STREAM = static_cast<void*>(stropen(const_cast<char*>(file),
				      const_cast<char*>(control)));
  if     (0 == std::strcmp(const_cast<char*>(control), "r"))
    get_history(mystream);
  else if(0 == std::strcmp(const_cast<char*>(control), "w")  ||
	  0 == std::strcmp(const_cast<char*>(control), "w!") ||
	  0 == std::strcmp(const_cast<char*>(control), "s") )
    put_history(mystream);
}  
//------------------------------------------------------------------------------
void nemo_io::close()
{
  if(STREAM) strclose(mystream);
  STREAM = 0;
}  
//------------------------------------------------------------------------------
#define ISPRESENT(NAME)						\
bool nemo_io::is_present(const NAME X) const			\
{ return get_tag_ok(mystream,tag(X)); }

ISPRESENT(Set)
ISPRESENT(SingleScalar)
ISPRESENT(SingleVector)
ISPRESENT(SingleMatrix)
ISPRESENT(SinglePhases)
ISPRESENT(BodiesScalar)
ISPRESENT(SPHScalar)
ISPRESENT(BodiesVector)
ISPRESENT(BodiesPhases)
ISPRESENT(BodiesInteger)
ISPRESENT(BodiesShort)
#undef ISPRESENT
//------------------------------------------------------------------------------
nemo_io::nemo_io(const char* file, const char* control) :
  STREAM ( static_cast<void*>(stropen(const_cast<char*>(file),
				      const_cast<char*>(control))) ),
  N             (0),
  NS            (0),
  BODIESSCALAR  (0),
  BODIESARRAYS  (0),
  BODIESINTEGER (0),
  BODIESSHORT   (0)
{
  OPEN[0]= OPEN[1] = OPEN[2] = OPEN[3] = false;
  if     (0 == std::strcmp(const_cast<char*>(control), "r"))
    get_history(mystream);
  else if(0 == std::strcmp(const_cast<char*>(control), "w")  ||
	  0 == std::strcmp(const_cast<char*>(control), "w!") ||
	  0 == std::strcmp(const_cast<char*>(control), "s") )
    put_history(mystream);
}
//------------------------------------------------------------------------------
nemo_io::~nemo_io()
{
  close();
  reset();
}
//------------------------------------------------------------------------------
void nemo_io::open_set(const Set S, const bool get) const
{
  if(get) get_set(mystream,tag(S));
  else    put_set(mystream,tag(S));
  OPEN[S] = true;
}
//------------------------------------------------------------------------------
void nemo_io::close_set(const Set S, const bool get) const
{
  if(get) get_tes(mystream,tag(S));
  else    put_tes(mystream,tag(S));
  OPEN[S] = false;
}
//------------------------------------------------------------------------------
void nemo_io::read_N() const
{
  reset();
  get_data(mystream,NobjTag,IntType,&N,0);
  if(get_tag_ok(mystream,NGasTag))
    get_data(mystream,NGasTag,IntType,&NS,0);
  else NS=0;
}
//------------------------------------------------------------------------------
double nemo_io::read(const SingleScalar X) const {
  double scal;
  get_data_coerced(mystream,tag(X),DoubleType,&scal,0);
  return scal;
}
void nemo_io::read(const SingleScalar X, double*Y) const {
  get_data_coerced(mystream,tag(X),DoubleType,Y,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const SingleVector X) const {
  get_data_coerced(mystream,tag(X),FalconType,SINGLEVECTOR,NDM,0); 
}
void nemo_io::read(const SingleVector X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,NDM,0); 
}
//------------------------------------------------------------------------------
void nemo_io::read(const SinglePhases X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,2,NDM,0);
}
void nemo_io::read(const SinglePhases X) const {
  get_data_coerced(mystream,tag(X),FalconType,SINGLEPHASES,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const SingleMatrix X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,NDM,NDM,0);
}
void nemo_io::read(const SingleMatrix X) const {
  get_data_coerced(mystream,tag(X),FalconType,SINGLEMATRIX,NDM,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesScalar X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,N,0);
}
void nemo_io::read(const BodiesScalar X) const {
  allocscalar();
  get_data_coerced(mystream,tag(X),FalconType,BODIESSCALAR,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const SPHScalar X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,NS,0);
}
void nemo_io::read(const SPHScalar X) const {
  allocscalar();
  get_data_coerced(mystream,tag(X),FalconType,BODIESSCALAR,NS,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesVector X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,N,NDM,0);
}
void nemo_io::read(const BodiesVector X) const {
  allocarrays();
  if(X==vel)
    get_data_coerced(mystream,tag(X),FalconType,BODIESARRAYS+N*NDM,N,NDM,0);
  else
    get_data_coerced(mystream,tag(X),FalconType,BODIESARRAYS,N,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesPhases X, real*Y) const {
  get_data_coerced(mystream,tag(X),FalconType,Y,N,2,NDM,0);
}
void nemo_io::read(const BodiesPhases X) const {
  allocarrays();
  get_data_coerced(mystream,tag(X),FalconType,BODIESARRAYS,N,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesInteger X, int*Y) const {
  get_data_coerced(mystream,tag(X),IntType,Y,N,0);
}
void nemo_io::read(const BodiesInteger X) const {
  allocinteger();
  get_data_coerced(mystream,tag(X),IntType,BODIESINTEGER,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesShort X, short*Y) const {
  get_data_coerced(mystream,tag(X),ShortType,Y,N,0);
}
void nemo_io::read(const BodiesShort X) const {
  allocshort();
  get_data_coerced(mystream,tag(X),ShortType,BODIESSHORT,N,0);
}
//==============================================================================
void nemo_io::write_N(int const&n, int const&ns) const
{
  reset();
  N  = n;
  NS = ns;
  put_data(mystream,NobjTag,IntType,&N,0);
  put_data(mystream,NGasTag,IntType,&NS,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const CoSys X) const
{
  CS = CSCode(tag(X),NDM,2);
  put_data(mystream,CoordSystemTag,IntType,&CS,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleScalar X, real const&scal) const
{
  double s = scal;
  put_data(mystream,tag(X),FalconType,&s,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleScalar X, double const&scal) const
{
  put_data(mystream,tag(X),DoubleType,const_cast<double*>(&scal),0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleVector X, const real* Y) const {
  put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),NDM,0);
}
void nemo_io::write(const SingleVector X) const {
  put_data(mystream,tag(X),FalconType,SINGLEVECTOR,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SinglePhases X, const real* Y) const {
  put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),2,NDM,0);
}
void nemo_io::write(const SinglePhases X) const {
  put_data(mystream,tag(X),FalconType,SINGLEPHASES,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleMatrix X, const real* Y) const {
  put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),NDM,NDM,0);
}
void nemo_io::write(const SingleMatrix X) const {
  put_data(mystream,tag(X),FalconType,SINGLEMATRIX,NDM,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesScalar X, const real* Y) const {
  put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),N,0);
}
void nemo_io::write(const BodiesScalar X) const {
  if(BODIESSCALAR==0)
    nbdy::error("[nemo_io::write(BodiesScalar)]: no memory allocated");
  else  put_data(mystream,tag(X),FalconType,BODIESSCALAR,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SPHScalar X, const real* Y) const {
  put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),N,0);
}
void nemo_io::write(const SPHScalar X) const {
  if(BODIESSCALAR==0)
    nbdy::error("[nemo_io::write(BodiesScalar)]: no memory allocated");
  else  put_data(mystream,tag(X),FalconType,BODIESSCALAR,NS,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesVector X, const real* Y) const {
  if(Y) put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),N,NDM,0);
}
void nemo_io::write(const BodiesVector X) const {
  if(BODIESARRAYS==0)
    nbdy::error("[nemo_io::write(BodiesVector)]: no memory allocated");
  if(X==vel) put_data(mystream,tag(X),FalconType,BODIESARRAYS+N*NDM,
		      N,NDM,0);
  else       put_data(mystream,tag(X),FalconType,BODIESARRAYS,N,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesPhases X, const real* Y) const {
  put_data(mystream,tag(X),FalconType,const_cast<real*>(Y),N,2,NDM,0);
}
void nemo_io::write(const BodiesPhases X) const {
  if(BODIESARRAYS==0)
    nbdy::error("[nemo_io::write(BodiesPhases)]: no memory allocated");
  else  put_data(mystream,tag(X),FalconType,BODIESARRAYS,N,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesInteger X, const int* Y) const {
  put_data(mystream,tag(X),IntType,const_cast<int*>(Y),N,0);
}
void nemo_io::write(const BodiesInteger X) const {
  if(BODIESINTEGER==0)
    nbdy::error("[nemo_io::write(BodiesInteger)]: no memory allocated");
  else  put_data(mystream,tag(X),IntType,BODIESINTEGER,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesShort X, const short* Y) const {
  put_data(mystream,tag(X),ShortType,const_cast<short*>(Y),N,0);
}
void nemo_io::write(const BodiesShort X) const {
  if(BODIESSHORT==0)
    nbdy::error("[nemo_io::write(BodiesShort)]: no memory allocated");
  else  put_data(mystream,tag(X),ShortType,BODIESSHORT,N,0);
}
#undef mystream
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  static int within_count = 0;
  inline bool within(real const&val, char*range, real const&fuzz)
    // almost identical to nemo::within
  {
    char*endptr, *subptr, *sepptr, *colptr;
    real sublow, subhi;
    int  count;
    subptr = range;
    if (*subptr++ == '#') {                        // special select Nth        
      within_count++;                              // (first=1) occurence       
      count = atoi(subptr);
      return count == within_count;
    }
    endptr = range + strlen(range);                // point to term. NULL       
    for(subptr = range; subptr != endptr; ) {      // for each subrange         
        sepptr = strchr(subptr, ',');              //   pnt to subrange end     
        if (sepptr == 0)                           //   last subrange listed?   
            sepptr = endptr;                       //     fix up subend ptr     
        colptr = strchr(subptr, ':');              //   scan subrange for :     
        if (colptr > sepptr)                       //   in another subrange?    
            colptr = 0;                            //     then dont use it      
        sublow = atof(subptr) - fuzz/2.0;          //   set low end of range    
        if (colptr != 0)                           //   high end specified?     
            subhi = atof(colptr+1) + fuzz/2.0;     //     set high end          
        else
            subhi = sublow + fuzz;                 //     just use low end      
        if (sublow <= val && val <= subhi)         //   within subrange?        
            return true;
        subptr = sepptr;                           //   advance subrange ptr    
        if (*subptr == ',')                        //   more ranges to do?      
            subptr++;                              //     move on to next       
    }
    return false;
  }
}
////////////////////////////////////////////////////////////////////////////////
bool nbdy::time_in_range(nbdy::real const &t, const char*times)
{
  return  times == 0
    ||    std::strcmp(const_cast<char*>(times),"all") == 0
    ||    nbdy::within(t,const_cast<char*>(times),0.0005);
}
////////////////////////////////////////////////////////////////////////////////
#endif                                               // falcON_NEMO             
