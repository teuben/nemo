//-----------------------------------------------------------------------------+
//                                                                             |
// nmio.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2001                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/nmio.h>
#ifdef   ALLOW_NEMO                                // empty if NEMO not running 
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
  bool within(::real, char*, ::real);              // NEMO range checker        
}
#ifndef DensityTag
#  define DensityTag "Density"
#endif
#ifndef EpsTag
#  define EpsTag "Eps"
#endif
#ifndef CPUTimeTag
#  define CPUTimeTag "cputime"
#endif
#ifndef PosTag
#  define PosTag "Position"
#endif
#ifndef VelTag
#  define VelTag "Velocity"
#endif
#ifndef FlagTag
#  define FlagTag "Flag"
#endif
#ifndef LevelTag
#  define LevelTag "Level"
#endif

#ifdef TWODIMENSIONAL
#  define NDM 2
#else
#  define NDM 3
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
  case nemo_io::key:  return KeyTag;
  case nemo_io::flag: return FlagTag;
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
ISPRESENT(BodiesVector)
ISPRESENT(BodiesPhases)
ISPRESENT(BodiesInteger)
ISPRESENT(BodiesShort)
#undef ISPRESENT
//------------------------------------------------------------------------------
nemo_io::nemo_io(const char* file, const char* control) :
  STREAM ( static_cast<void*>(stropen(const_cast<char*>(file),
				      const_cast<char*>(control))) ),
  N ( 0 ),
  BODIESSCALAR  (0),
  BODIESARRAYS  (0),
  BODIESINTEGER (0),
  BODIESSHORT   (0)
{
  OPEN[0]= OPEN[1] = OPEN[2] = OPEN[3] = false;
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
int nemo_io::read_N() const
{
  reset();
  get_data(mystream,NobjTag,IntType,&N,0);
  return N;
}
//------------------------------------------------------------------------------
float nemo_io::read(const SingleScalar X, float*Y) const
{
  if(Y) {
    get_data_coerced(mystream,tag(X),FloatType,Y,0);
    return *Y;
  } else {
    float scal;
    get_data_coerced(mystream,tag(X),FloatType,&scal,0);
    return scal;
  }
}
//------------------------------------------------------------------------------
void nemo_io::read(const SingleVector X, float*Y) const
{
  if(Y) get_data_coerced(mystream,tag(X),FloatType,Y,NDM,0);
  else  get_data_coerced(mystream,tag(X),FloatType,SINGLEVECTOR,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const SinglePhases X, float*Y) const
{
  if(Y) get_data_coerced(mystream,tag(X),FloatType,Y,2,NDM,0);
  else  get_data_coerced(mystream,tag(X),FloatType,SINGLEPHASES,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const SingleMatrix X, float*Y) const
{
  if(Y) get_data_coerced(mystream,tag(X),FloatType,Y,NDM,NDM,0);
  else  get_data_coerced(mystream,tag(X),FloatType,SINGLEMATRIX,NDM,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesScalar X, float*Y) const
{
  if(Y)
    get_data_coerced(mystream,tag(X),FloatType,Y,N,0);
  else {
    allocscalar();
    get_data_coerced(mystream,tag(X),FloatType,BODIESSCALAR,N,0);
  }
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesVector X, float*Y) const
{
  if(Y)
    get_data_coerced(mystream,tag(X),FloatType,Y,N,NDM,0);
  else {
    allocarrays();
    if(X==vel)
      get_data_coerced(mystream,tag(X),FloatType,BODIESARRAYS+N*NDM,N,NDM,0);
    else
      get_data_coerced(mystream,tag(X),FloatType,BODIESARRAYS,N,NDM,0);
  }
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesPhases X, float*Y) const
{
  if(Y)
    get_data_coerced(mystream,tag(X),FloatType,Y,N,2,NDM,0);
  else {
    allocarrays();
    get_data_coerced(mystream,tag(X),FloatType,BODIESARRAYS,N,2,NDM,0);
  }
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesInteger X, int*Y) const
{
  if(Y)
    get_data_coerced(mystream,tag(X),IntType,Y,N,0);
  else {
    allocinteger();
    get_data_coerced(mystream,tag(X),IntType,BODIESINTEGER,N,0);
  }
}
//------------------------------------------------------------------------------
void nemo_io::read(const BodiesShort X, short*Y) const
{
  if(Y)
    get_data_coerced(mystream,tag(X),ShortType,Y,N,0);
  else {
    allocshort();
    get_data_coerced(mystream,tag(X),ShortType,BODIESSHORT,N,0);
  }
}
//------------------------------------------------------------------------------
void nemo_io::read_history() const
{
  get_history(mystream);
}
//------------------------------------------------------------------------------
void nemo_io::write_N(const int n) const
{
  reset();
  N = n;
  put_data(mystream,NobjTag,IntType,&N,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const CoSys X) const
{
  CS = CSCode(tag(X), NDM, 2);
  put_data(mystream,CoordSystemTag,IntType,&CS,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleScalar X, const float scal) const
{
  float s = scal;
  put_data(mystream,tag(X),FloatType,&s,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleVector X, float* Y) const
{
  if(Y) put_data(mystream,tag(X),FloatType,Y,NDM,0);
  else  put_data(mystream,tag(X),FloatType,SINGLEVECTOR,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SinglePhases X, float* Y) const
{
  if(Y) put_data(mystream,tag(X),FloatType,Y,2,NDM,0);
  else  put_data(mystream,tag(X),FloatType,SINGLEPHASES,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const SingleMatrix X, float* Y) const
{
  if(Y) put_data(mystream,tag(X),FloatType,Y,NDM,NDM,0);
  else  put_data(mystream,tag(X),FloatType,SINGLEMATRIX,NDM,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesScalar X, float* Y) const
{
  if(Y) put_data(mystream,tag(X),FloatType,Y,N,0);
  else if(BODIESSCALAR==0)
    nbdy::error("[nemo_io::write(BodiesScalar)]: no memory allocated");
  else  put_data(mystream,tag(X),FloatType,BODIESSCALAR,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesVector X, float* Y) const
{
  if(Y)           put_data(mystream,tag(X),FloatType,Y,N,NDM,0);
  else if(BODIESARRAYS==0)
    nbdy::error("[nemo_io::write(BodiesVector)]: no memory allocated");
  else if(X==vel) put_data(mystream,tag(X),FloatType,BODIESARRAYS+N*NDM,N,NDM,0);
  else            put_data(mystream,tag(X),FloatType,BODIESARRAYS,N,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesPhases X, float* Y) const
{
  if(Y) put_data(mystream,tag(X),FloatType,Y,N,2,NDM,0);
  else if(BODIESARRAYS==0)
    nbdy::error("[nemo_io::write(BodiesPhases)]: no memory allocated");
  else  put_data(mystream,tag(X),FloatType,BODIESARRAYS,N,2,NDM,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesInteger X, int* Y) const
{
  if(Y) put_data(mystream,tag(X),IntType,Y,N,0);
  else if(BODIESINTEGER==0)
    nbdy::error("[nemo_io::write(BodiesInteger)]: no memory allocated");
  else  put_data(mystream,tag(X),IntType,BODIESINTEGER,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::write(const BodiesShort X, short* Y) const
{
  if(Y) put_data(mystream,tag(X),ShortType,Y,N,0);
  else if(BODIESSHORT==0)
    nbdy::error("[nemo_io::write(BodiesShort)]: no memory allocated");
  else  put_data(mystream,tag(X),ShortType,BODIESSHORT,N,0);
}
//------------------------------------------------------------------------------
void nemo_io::write_history() const
{
  put_history(mystream);
}
//------------------------------------------------------------------------------
#undef mystream
#undef NDM
////////////////////////////////////////////////////////////////////////////////
bool nbdy::time_in_range(const nbdy::real& t, const char*times)
{
  return  times == 0
    ||    std::strcmp(const_cast<char*>(times),"all") == 0
    ||    within(t,const_cast<char*>(times),0.0005);
}
////////////////////////////////////////////////////////////////////////////////
#endif                                               // ALLOW_NEMO              
