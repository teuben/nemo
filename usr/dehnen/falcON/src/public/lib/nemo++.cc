// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/nemo++.cc
///
/// \author  Walter Dehnen
///
/// \date    2008-2009
///
/// \brief   implements methods declared in inc/nemo++.h
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008-2009  Walter Dehnen
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
#include <public/nemo++.h>
#include <utils/io.h>

#ifdef falcON_NEMO

extern "C" {
#  include <public/basic.h>
#  include <stdinc.h>
#  include <loadobj.h>
#  include <getparam.h>
#  include <history.h>
#  include <filestruct.h>
#  include <filefn.h>
#  include <snapshot/snapshot.h>
}

#undef getrparam
#undef nemoinpr
#undef getargv0
#undef getversion

////////////////////////////////////////////////////////////////////////////////
//
// NOTE
//
// There is a good reason for not making namespace falcON available in the
// global or anonymous namespace (via "using namespace falcON;"):
//
// If we did and forget to include all relevant nemo header files, then the
// definition of falcON routines via nemo routines (which are in the global
// namespace) with identical name and syntax will result in an infinite
// recursive loop (because the nemo routine is accidentally not declared, so
// instead it uses the routine to be defined itself), instead of a compile
// time error!
//
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  extern int debug_level;         // NEMO's debug_level, see nemo's dprintf.c
}

int falcON::nemo_debug_level() {
  return debug_level;
}

////////////////////////////////////////////////////////////////////////////////

namespace {
  template<typename Type> struct __inpA;
  template<> struct __inpA<double> {
    static int inp(const char*p, double*a, int m) {
      return ::nemoinpd(const_cast<char*>(p),a,m);
    } };
  template<> struct __inpA<float> {
    static int inp(const char*p, float*a, int m) {
      return ::nemoinpf(const_cast<char*>(p),a,m);
    } };
  template<> struct __inpA<int> {
    static int inp(const char*p, int*a, int m) {
      return ::nemoinpi(const_cast<char*>(p),a,m);
    } };
  template<> struct __inpA<long> {
    static int inp(const char*p, long*a, int m) {
      return ::nemoinpl(const_cast<char*>(p),a,m);
    } };
  template<> struct __inpA<unsigned> {
    static int inp(const char*p, unsigned*a, int m) {
      return ::nemoinpi(const_cast<char*>(p),
			static_cast<int*>(static_cast<void*>(a)),m);
    } };
  template<> struct __inpA<bool> {
    static int inp(const char*p, bool*a, int m) {
      return ::nemoinpb(const_cast<char*>(p),a,m);
    } };
}

template<typename T>
int falcON::nemoinp(const char*e, T*x, int n) {
  return __inpA<T>::inp(e,x,n);
}

template int falcON::nemoinp(const char*, bool    *, int);
template int falcON::nemoinp(const char*, int     *, int);
template int falcON::nemoinp(const char*, long    *, int);
template int falcON::nemoinp(const char*, unsigned*, int);
template int falcON::nemoinp(const char*, float   *, int);
template int falcON::nemoinp(const char*, double  *, int);

int nemoinpx(const char*e, double*x, int n) {
  return ::nemoinpx(const_cast<char*>(e),x,n);
}

////////////////////////////////////////////////////////////////////////////////

void falcON::initparam(const char**argv, const char**defv) {
  ::initparam(const_cast<char**>(argv),const_cast<char**>(defv));
}

void falcON::finiparam() {
  ::finiparam();
}

bool falcON::hasvalue(const char*p) {
  return ::hasvalue(const_cast<char*>(p));
}

const char*falcON::getparam(const char*p) {
  return ::getparam(const_cast<char*>(p));
}

bool falcON::getbparam(const char*p) {
  return ::getbparam(const_cast<char*>(p));
}

int falcON::getiparam(const char*p) {
  return ::getiparam(const_cast<char*>(p));
}

long falcON::getlparam(const char*p) {
  return ::getlparam(const_cast<char*>(p));
}

double falcON::getdparam(const char*p) {
  return ::getdparam(const_cast<char*>(p));
}

falcON::vect falcON::getvparam(const char*p) falcON_THROWING {
  vect X;
  int  N = nemoinp(getparam(p), static_cast<real*>(X), Ndim);
  if(N!=Ndim) {
    if(N<0) falcON_THROW("parse error: processing parameter \"%s\"\n",p);
    else    falcON_THROW("parameter \"%s\" requires %d values, but %d given\n",
			 p,Ndim,N);
  }
  return X;
}

falcON::vect falcON::getvrparam(const char* p) falcON_THROWING {
  vect X;
  int  N = nemoinp(getparam(p), static_cast<real*>(X), Ndim);
  if(N==1) for(int d=1; d!=Ndim; ++d) X[d]=X[0];
  else if(N!=Ndim) {
    if(N<0) falcON_THROW("parse error: processing parameter \"%s\"\n",p);
    else    falcON_THROW("parameter \"%s\" requires %d values or 1, "
			 "but %d given\n", p,Ndim,N);
  }
  return X;
}

falcON::vect* falcON::getvparam_z(const char* p, vect&X) falcON_THROWING {
  if(!hasvalue(p)) return 0;
  int N = nemoinp(getparam(p), static_cast<real*>(X), Ndim);
  if(Ndim != N) {
    if(N<0) falcON_THROW("parse error: processing parameter \"%s\"\n",p);
    else  falcON_Warning("parameter \"%s\" requires %d values, but %d given\n",
			 p,Ndim,N);
    return 0;
  }
  return &X;
}

falcON::vect* falcON::getvrparam_z(const char* p, vect&X) falcON_THROWING {
  if(!hasvalue(p)) return 0;
  int  N = nemoinp(getparam(p), static_cast<real*>(X), Ndim);
  if(N==1) for(int d=1; d!=Ndim; ++d) X[d]=X[0];
  else if(N!=Ndim) {
    if(N<0) falcON_THROW("parse error: processing parameter \"%s\"\n",p);
    else  falcON_Warning("parameter \"%s\" requires %d values or 1, "
			 "but %d given\n", p,Ndim,N);
    return 0;
  }
  return &X;
}

////////////////////////////////////////////////////////////////////////////////

const char*falcON::fullname(const char*file) {
  return ::fullname(const_cast<char*>(file));
}

const char*falcON::pathfind(const char*path, const char*file) {
  return ::pathfind(const_cast<char*>(path),const_cast<char*>(file));
}

void falcON::mysymbols(const char*p) {
  ::mysymbols(const_cast<char*>(p));
}

void falcON::loadobj(const char*f) {
  ::loadobj(const_cast<char*>(f));
}

void(*falcON::findfn(const char*f))() {
  return ::findfn(const_cast<char*>(f));
}

void falcON::mapsys(char*f) {
  ::mapsys(f);
}

////////////////////////////////////////////////////////////////////////////////
//
// nemo I/O
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NSinkTag
# define NSinkTag "NSink"
#endif
#ifndef SoundSpeedTag
# define SoundSpeedTag "SoundSpeed"
#endif
#ifndef TimeStepTag
# define TimeStepTag "TimeStep"
#endif
#ifndef OrbitalPeriodTag
# define OrbitalPeriodTag "OrbitalPeriod"
#endif
#ifndef LevelTag
# define LevelTag "Level"
#endif
#ifndef JerkTag
# define JerkTag "Jerk"
#endif
#ifndef PhaseSpaceDensityTag
# define PhaseSpaceDensityTag "PhaseSpaceDensity"
#endif
#ifndef GasHdotTag
# define GasHdotTag "Hdot"
#endif
#ifndef GasFactTag
# define GasFactTag "SPHFactor"
#endif
#ifndef MolWeightTag
# define MolWeightTag "MolecularWeight"
#endif
#ifndef ArtViscTag
# define ArtViscTag "ArtificialViscosity"
#endif
#ifndef GasDivVTag
# define GasDivVTag "Divergence(Velocity)"
#endif
#ifndef SpinTag
# define SpinTag "SpinVector"
#endif

namespace {
  inline const char* NemoTag(falcON::nemo_io::Field f) falcON_THROWING
  {
    switch(f) {
    case falcON::nemo_io::mass   : return MassTag;
    case falcON::nemo_io::pos    : return PosTag;
    case falcON::nemo_io::vel    : return VelTag;
    case falcON::nemo_io::eps    : return EpsTag;
    case falcON::nemo_io::key    : return KeyTag;
    case falcON::nemo_io::step   : return TimeStepTag;
    case falcON::nemo_io::pot    : return PotentialTag;
    case falcON::nemo_io::acc    : return AccelerationTag;
    case falcON::nemo_io::jerk   : return JerkTag;
    case falcON::nemo_io::dens   : return DensityTag;
    case falcON::nemo_io::aux    : return AuxTag;
    case falcON::nemo_io::zet    : return AuxVecTag;
    case falcON::nemo_io::lev    : return LevelTag;
    case falcON::nemo_io::num    : return NumberTag;
    case falcON::nemo_io::phden  : return PhaseSpaceDensityTag;
    case falcON::nemo_io::torb   : return OrbitalPeriodTag;
    case falcON::nemo_io::posvel : return PhaseSpaceTag;
    case falcON::nemo_io::Size   : return SmoothTag;
    case falcON::nemo_io::Gasnum : return SPHNumberTag;
    case falcON::nemo_io::Uin    : return UinternTag;
    case falcON::nemo_io::Uindot : return UdotIntTag;
    case falcON::nemo_io::Uinrad : return UdotRadTag;
    case falcON::nemo_io::Entr   : return EntFuncTag;
    case falcON::nemo_io::Gasdens: return GasDensTag;
    case falcON::nemo_io::Sizedot: return GasHdotTag;
    case falcON::nemo_io::Sphfact: return GasFactTag;
    case falcON::nemo_io::Csound : return SoundSpeedTag;
    case falcON::nemo_io::AlphaAV: return ArtViscTag;
    case falcON::nemo_io::DivV   : return GasDivVTag;
    case falcON::nemo_io::MolWght: return MolWeightTag;
    case falcON::nemo_io::Spin   : return SpinTag;
    case falcON::nemo_io::null:
      falcON_THROW("nemo I/O: nemo_io::null not I/O able");
    default:
      falcON_THROW("nemo I/O: unknown nemo_io::Field '%d'",f);
    }
  }
  inline falcON::nemo_io::DataType Type(const char* NemoType)
  {
    if(0==std::strcmp(NemoType, ByteType  ) ) return falcON::nemo_io::Byte;
    if(0==std::strcmp(NemoType, ShortType ) ) return falcON::nemo_io::Short;
    if(0==std::strcmp(NemoType, IntType   ) ) return falcON::nemo_io::Integer;
    if(0==std::strcmp(NemoType, LongType  ) ) return falcON::nemo_io::Long;
    if(0==std::strcmp(NemoType, FloatType ) ) return falcON::nemo_io::Single;
    if(0==std::strcmp(NemoType, DoubleType) ) return falcON::nemo_io::Double;
    return falcON::nemo_io::Null;
  }
  inline const char* NemoType(falcON::nemo_io::DataType t)
  {
    switch(t) {
    case falcON::nemo_io::Byte   : return ByteType;
    case falcON::nemo_io::Short  : return ShortType;
    case falcON::nemo_io::Integer: return IntType;
    case falcON::nemo_io::Long   : return LongType;
    case falcON::nemo_io::Single : return FloatType;
    case falcON::nemo_io::Double : return DoubleType;
    default              : return AnyType;
    }
  }
  inline bool is_scalar(falcON::nemo_io::Field f) {
    switch(f) {
    case falcON::nemo_io::mass:
    case falcON::nemo_io::eps:
    case falcON::nemo_io::key:
    case falcON::nemo_io::step:
    case falcON::nemo_io::pot:
    case falcON::nemo_io::dens:
    case falcON::nemo_io::aux:
    case falcON::nemo_io::lev:
    case falcON::nemo_io::num:
    case falcON::nemo_io::phden:
    case falcON::nemo_io::torb:
    case falcON::nemo_io::Size:
    case falcON::nemo_io::Gasnum:
    case falcON::nemo_io::Uin:
    case falcON::nemo_io::Uindot:
    case falcON::nemo_io::Uinrad:
    case falcON::nemo_io::Entr:
    case falcON::nemo_io::Gasdens:
    case falcON::nemo_io::Sizedot:
    case falcON::nemo_io::Sphfact:
    case falcON::nemo_io::Csound:  
    case falcON::nemo_io::AlphaAV:  
    case falcON::nemo_io::DivV:  
    case falcON::nemo_io::MolWght: return true;
    default:                       return false;
    }
  }
  inline bool is_vector(falcON::nemo_io::Field f) {
    switch(f) {
    case falcON::nemo_io::pos:
    case falcON::nemo_io::vel:
    case falcON::nemo_io::acc:
    case falcON::nemo_io::zet:
    case falcON::nemo_io::jerk: 
    case falcON::nemo_io::Spin: return true;
    default:                    return false;
    }
  }
  inline bool is_phases(falcON::nemo_io::Field f) {
    return f == falcON::nemo_io::posvel;
  }
} // namespace {
//------------------------------------------------------------------------------
/// define some auxiliary functionality
/// \note most of this can be scrapped once NEMO is fixed
namespace Aux {
  std::FILE* stropen(const char*file, const char*mode) {
    return ::stropen(file,const_cast<char*>(mode));
  }
  inline void put_set(std::FILE* file, const char*set) {
    ::put_set(file, const_cast<char*>(set));
  }
  inline void put_tes(std::FILE* file, const char*set) {
    ::put_tes(file, const_cast<char*>(set));
  }
  inline void put_data(std::FILE* file, const char*tag,
		       const char*type, const void*X,
		       int N0, int N1=0, int N2=0, int N3=0) {
    ::put_data(file,const_cast<char*>(tag),const_cast<char*>(type),
		 const_cast<void*>(X),N0,N1,N2,N3,0);
  }
  inline void put_data_set(std::FILE* file, const char*tag,
		       const char*type, int N0, int N1=0, int N2=0, int N3=0) {
    ::put_data_set(file,const_cast<char*>(tag),const_cast<char*>(type),
		   N0,N1,N2,N3,0);
  }
  inline void put_data_blocked(std::FILE* file, const char*tag, const void*X,
			       int N) {
    ::put_data_blocked(file,const_cast<char*>(tag),const_cast<void*>(X),N);
  }
  inline void put_data_tes(std::FILE* file, const char*tag) {
    ::put_data_tes(file,const_cast<char*>(tag));
  }
  inline bool get_tag_ok(std::FILE* file, const char*tag) {
    return ::get_tag_ok(file,const_cast<char*>(tag));
  }
  inline const char*get_type(std::FILE* file, const char*tag) {
    return ::get_type(file,const_cast<char*>(tag));
  }
  inline const int*get_dims(std::FILE* file, const char*tag) {
    return ::get_dims(file,const_cast<char*>(tag));
  }
  inline void get_set(std::FILE* file, const char*set) {
    ::get_set(file,const_cast<char*>(set));
  }
  inline void get_tes(std::FILE* file, const char*set) {
    ::get_tes(file,const_cast<char*>(set));
  }
  inline void get_data(std::FILE* file, const char*tag, const char*type,
		       void*X, int N0, int N1=0, int N2=0, int N3=0) {
    ::get_data(file,const_cast<char*>(tag),const_cast<char*>(type),X,
	       N0,N1,N2,N3,0);
  }
  inline void get_data_set(std::FILE* file, const char*tag, const char*type,
			   int N0, int N1=0, int N2=0, int N3=0) {
    ::get_data_set(file,const_cast<char*>(tag),const_cast<char*>(type),
		   N0,N1,N2,N3,0);
  }
  inline void get_data_blocked(std::FILE* file, const char*tag, void*X, int N) {
    ::get_data_blocked(file,const_cast<char*>(tag),X,N);
  }
  inline void get_data_tes(std::FILE* file, const char*tag) {
    ::get_data_tes(file,const_cast<char*>(tag));
  }
}
using namespace Aux;
//------------------------------------------------------------------------------
// class falcON::nemo_io
//------------------------------------------------------------------------------
void falcON::nemo_io::open(const char*file, const char*mode) falcON_THROWING
{
  close();
  if(file && file[0]) {
    PIPE = !std::strcmp(file,"-");
    SINK = !std::strcmp(file,".");
    if       (0 == std::strcmp(mode, "r")) {
      if(SINK) falcON_THROW("nemo_io::open(): "
			    "will not open file '%s' for input\n",file);
      if(PIPE) input::open_std();
      STREAM = stropen(file,mode);
      get_history(STREAM);
      IN = 1;
      DebugInfo(4,"nemo_io: opened file '%s' for input\n",file);
    } else if(0 == std::strcmp(mode, "w")  ||
	      0 == std::strcmp(mode, "w!") ||
	      0 == std::strcmp(mode, "a")  ||
	      0 == std::strcmp(mode, "a!") ||
	      0 == std::strcmp(mode, "s") ) {
      if(PIPE) output::open_std();
      STREAM = stropen(file,mode);
      put_history(STREAM);
      OUT  = 1;
      DebugInfo(4,"nemo_io: opened file '%s' for output with mode '%s'\n",
		file,mode);
    } else
      falcON_THROW("nemo_io::open(): unknown mode '%s'\n",mode);
  }
}
//------------------------------------------------------------------------------
void falcON::nemo_io::close()
{
  if(STREAM) ::strclose(STREAM);
  if(PIPE) {
    if(IN)  input ::close_std();
    if(OUT) output::close_std();
  }
  if(STREAM) DebugInfo(4,"nemo_io: closed %s stream\n", IN? "input" : "output");
  STREAM=0, IN=0, OUT=0, PIPE=0, SINK=0;
}
//------------------------------------------------------------------------------
// class falcON::nemo_in
//------------------------------------------------------------------------------
void falcON::nemo_in::close() falcON_THROWING
{
  if(SNAP) {
    DebugInfo(4,"nemo_in::close(): closing open snap_in first ... \n");
    SNAP->~snap_in();
    SNAP = 0;
  }
  nemo_io::close();
}
//------------------------------------------------------------------------------
bool falcON::nemo_in::has_snapshot() const
{
  // added to read any history padded between snapshots (note that gyrfalcON
  // does that upon appending to an existing file when resuming a simulation)
  get_history(STREAM);
  return STREAM && get_tag_ok(STREAM,SnapShotTag);
}
//------------------------------------------------------------------------------
// class falcON::snap_in
//------------------------------------------------------------------------------
falcON::snap_in::snap_in(nemo_in const&in) falcON_THROWING : 
  INPUT(in), DATA(0), FIELDS_READ(0), HAS_TIME(0), NTOT(0u), TIME(0.)
{
  DebugInfo(4,"snap_in::snap_in() ...\n");
  for(bodytype t; t; ++t) NBOD[t] = 0u;
  if(! INPUT.has_snapshot())
    falcON_THROW("cannot open snapshot from nemo input stream");
  if(INPUT.SNAP)
    falcON_THROW("trying to open 2nd snapshot from nemo input stream");
  // 1 open snapshot set
  get_set(INPUT.STREAM,SnapShotTag);
  INPUT.SNAP = this;
  DebugInfo(5,"  snap_in::snap_in(): snapshot opened\n");
  // 2 open parameter set
  if(!get_tag_ok(INPUT.STREAM,ParametersTag)) {
    get_tes(INPUT.STREAM,SnapShotTag);
    INPUT.SNAP = 0;
    falcON_THROW("cannot read parameters from nemo input stream");
  }
  get_set(INPUT.STREAM,ParametersTag);
  DebugInfo(5,"  snap_in::snap_in(): parameter set opened\n");
  // 3 read parameter set
  // 3.1 read total # bodies 
  if(!get_tag_ok(INPUT.STREAM,NobjTag)) {
    get_tes(INPUT.STREAM,ParametersTag);
    get_tes(INPUT.STREAM,SnapShotTag);
    INPUT.SNAP = 0;
    falcON_THROW("cannot read # bodies from nemo input stream");
  }
  get_data(INPUT.STREAM,NobjTag,IntType,&NTOT,0);
  DebugInfo(5,"  snap_in::snap_in(): read Nobj = %u\n",NTOT);
  // 3.2 try to read # SINK bodies
  if(get_tag_ok(INPUT.STREAM,NSinkTag)) {
    get_data(INPUT.STREAM,NSinkTag,IntType,&(NBOD[bodytype::sink]),0);
    DebugInfo(5,"  snap_in::snap_in(): read Nsink = %u\n",
	      NBOD[bodytype::sink]);
  }
  // 3.3 try to read # SPH bodies
  if(get_tag_ok(INPUT.STREAM,NGasTag)) {
    get_data(INPUT.STREAM,NGasTag,IntType,&(NBOD[bodytype::gas]),0);
    DebugInfo(5,"  snap_in::snap_in(): read Nsph = %u\n",
	      NBOD[bodytype::gas]);
  }
  // 3.4 set # STD bodies
  unsigned n(0u);
  for(bodytype t; t; ++t) n += NBOD[t];
  if(n > NTOT)
    falcON_THROW("read nemo data: more non-STD bodies than total");
  NBOD[bodytype::std] = NTOT - n;
  // 3.5 try to read simulation time
  if(get_tag_ok(INPUT.STREAM,TimeTag)) {
    HAS_TIME = true;
    const char* time_type = get_type(INPUT.STREAM,TimeTag);
    if(0 == std::strcmp(time_type,DoubleType))
	get_data(INPUT.STREAM,TimeTag,DoubleType,&TIME,0);
    else if(0 == std::strcmp(time_type,FloatType)) {
      float __TIME;
      get_data(INPUT.STREAM,TimeTag,FloatType,&__TIME,0);
      TIME = __TIME;
    } else
      falcON_Warning("nemo input: unknown type '%s' for time\n",time_type);
  }
  if(HAS_TIME)
    DebugInfo(5,"  read time = %f\n",TIME);
  // 4 close parameter set
  get_tes(INPUT.STREAM,ParametersTag);
  DebugInfo(5,"  snap_in::snap_in(): parameter set read & closed\n");
  // 5 open particle set
  if(!get_tag_ok(INPUT.STREAM,ParticlesTag)) {
    get_tes(INPUT.STREAM,SnapShotTag);
    INPUT.SNAP = 0;
    falcON_THROW("cannot open particle set from nemo input stream");
  }
  get_set(INPUT.STREAM,ParticlesTag);
  DebugInfo(5,"  snap_in::snap_in(): particles set opened\n");
}
//------------------------------------------------------------------------------
falcON::snap_in::~snap_in() falcON_THROWING
{
  if(DATA) {
    DebugInfo(4,"snap_in::~snap_in(): closing open data_in first\n");
    DATA->~data_in();
  }
  HAS_TIME = false;
  NTOT = 0;
  for(bodytype t; t; ++t) NBOD[t] = 0u;
  get_tes(INPUT.STREAM,ParticlesTag);
  get_tes(INPUT.STREAM,SnapShotTag);
  // added to read any history padded between snapshots (note that gyrfalcON
  // does that upon appending to an existing file when resuming a simulation)
  get_history(INPUT.STREAM);
  INPUT.SNAP = 0;
  DebugInfo(4,"snap_in: closed\n");
}
//------------------------------------------------------------------------------
bool falcON::snap_in::has(nemo_io::Field f) const
{
  return 
    !has_been_read(f) &&
    get_tag_ok(INPUT.STREAM,NemoTag(f));
}
//------------------------------------------------------------------------------
// class falcON::data_in
//------------------------------------------------------------------------------
namespace {
  using namespace falcON;
  typedef tupel<Ndim,notreal> Vect;
}
falcON::data_in::data_in(snap_in const&snap, nemo_io::Field f) falcON_THROWING :
  INPUT(snap), FIELD(f), NREAD(0), NTOT(0), TYPE(nemo_io::Null), SUBN(0)
{
  DebugInfo(5,"data_in::data_in(%s) ...\n",NemoTag(FIELD));
  if( INPUT.DATA )
    falcON_THROW("cannot read %s: nemo input still engaged",NemoTag(FIELD));
  if( !INPUT.has(f) )
    falcON_THROW("cannot read %s: not given with nemo input",NemoTag(FIELD));
  if( INPUT.has_been_read(f) )
    falcON_THROW("cannot read %s: already read from nemo input",
		    NemoTag(FIELD));
  // 1 get type of data on file and check for mismatch
  const char*TypeTag=get_type(INPUT.stream(),NemoTag(FIELD));
  TYPE = Type(TypeTag);
  if(nemo_io::mismatch(TYPE, nemo_io::type(FIELD)))
    falcON_THROW("cannot read %s: type mismatch (got %s, expect %s)",
		 NemoTag(FIELD), NemoType(TYPE), 
		 NemoType(nemo_io::type(FIELD)));
  DebugInfo(6,"  data type: %s\n",nemo_io::type_name(TYPE));
  // 2 get dimensions of data on file
  const int*dim=get_dims(INPUT.stream(),NemoTag(FIELD));
  if(!dim)
    falcON_THROW("cannot read # %s data",NemoTag(FIELD));
  NTOT = dim[0];
  // 3 check for dimensions to match expectations AND open data set
  if(NTOT != snap.N(FIELD) ) {
    // 3.1 dim[0] != snap_in::N: error out
    falcON_THROW("nemo input of %s: found %d data, expected %d",
		 NemoTag(FIELD), dim[0], snap.N(FIELD));
  } else if(dim[1] == 0) {
    // 3.2 array of scalars
    if(!::is_scalar(FIELD) )
      falcON_THROW("nemo input of %s: found scalars",NemoTag(FIELD));
    DebugInfo(6,"  opening data set for %d scalars\n",NTOT);
    get_data_set(INPUT.stream(),NemoTag(FIELD),TypeTag,NTOT);
    SUBN = 1;
  } else if(dim[2] == 0) {
    // 3.3 array of vectors
    if(!::is_vector(FIELD) )
      falcON_THROW("nemo input of %s: found vectors",NemoTag(FIELD));
    if(dim[1] != NDIM)
      falcON_THROW("nemo input of %s: Ndim mismatch",NemoTag(FIELD));
    DebugInfo(6,"  opening data set for %d vectors\n",NTOT);
    get_data_set(INPUT.stream(),NemoTag(FIELD),TypeTag,NTOT,NDIM);
    SUBN = NDIM;
  } else if(dim[3] == 0) {
    // 3.4 array of phases
    if(!is_phases(FIELD) )
      falcON_THROW("nemo input of %s: found phases",NemoTag(FIELD));
    if(dim[1] != 2 && dim[2] != NDIM)
      falcON_THROW("nemo input of %s: Ndim mismatch",NemoTag(FIELD));
    DebugInfo(6,"  opening data set for %d phases\n",NTOT);
    get_data_set(INPUT.stream(),NemoTag(FIELD),TypeTag,NTOT,2,NDIM);
    SUBN = 2*NDIM;
  } else {
    // 3.5 array of yet higher rank
    falcON_THROW("nemo input of %s: found high-rank data",NemoTag(FIELD));
  }
  INPUT.DATA = this;
}
//------------------------------------------------------------------------------
falcON::data_in::~data_in()
{
  get_data_tes(INPUT.stream(),NemoTag(FIELD));
  INPUT.DATA = 0;
  INPUT.FIELDS_READ |= FIELD;
  DebugInfo(5,"data_in(%s) closed\n",NemoTag(FIELD));
}
//------------------------------------------------------------------------------
void falcON::data_in::read(void*data, unsigned N)
{
  if(NREAD >= NTOT) {
    falcON_Warning("nemo input of %s: cannot read any more (all %d read)\n",
		   NemoTag(FIELD),NREAD); 
    return;
  }
  if(N == 0)
    N = NTOT - NREAD;
  else if(NREAD + N > NTOT) {
    falcON_Warning("nemo input of %s: cannot read %d, only %d data left",
		   NemoTag(FIELD), N, NTOT-NREAD);
    N = NTOT - NREAD;
  }
  if(nemo_io::coercing(TYPE, nemo_io::type(FIELD))) {
    DebugInfo(1,"data_in::read(%s): must coerce\n",NemoTag(FIELD));
    unsigned n = N*SUBN;
    notreal*buf= falcON_NEW(notreal, N);
    get_data_blocked(INPUT.stream(),NemoTag(FIELD), buf, n);
    for(unsigned i=0; i!=n; ++i) static_cast<real*>(data)[i] = buf[i];
    falcON_DEL_A(buf);
  } else 
    get_data_blocked(INPUT.stream(),NemoTag(FIELD), data, N*SUBN);
  DebugInfo(5,"data_in::read(): %d %s read\n",N,NemoTag(FIELD));
  NREAD += N;
}
//------------------------------------------------------------------------------
void falcON::data_in::read_phases(void*pos, void*vel, unsigned N)
{
  if(FIELD != nemo_io::posvel)
    falcON_THROW("data_in::read_phases(%s)\n",NemoTag(FIELD));
  if(pos == 0 && vel == 0) {
    falcON_Warning("data_in::read_phases(): pos=%p, vel=%p\n",pos,vel);
    return;
  }
  if(NREAD >= NTOT) {
    falcON_Warning("data_in::read_phases() cannot read any more "
		   "(all %d read)\n", NREAD); 
    return;
  }
  if(N == 0)
    N = NTOT - NREAD;
  else if(NREAD + N > NTOT) {
    falcON_Warning("nemo input of %s: cannot read %d, only %d data left",
		   NemoTag(FIELD), N, NTOT-NREAD);
    N = NTOT - NREAD;
  }
  const bool coerce = nemo_io::coercing(TYPE, nemo_io::type(FIELD));
  if(coerce)
    DebugInfo(1,"data_in::read_phases(): must coerce\n");
  void* phases = coerce?
    static_cast<void*>(falcON_NEW(Vect,2*N)) : 
    static_cast<void*>(falcON_NEW(vect,2*N)) ;
  get_data_blocked(INPUT.stream(),NemoTag(FIELD), phases, N*SUBN);
  if(pos) {
    vect*to = static_cast<vect*>(pos);
    if(coerce) {
      const Vect*ph = static_cast<const Vect*>(phases);
      for(unsigned n=0; n!=N; ++n,++to,ph+=2) *to = *ph;
    } else {
      const vect*ph = static_cast<const vect*>(phases);
      for(unsigned n=0; n!=N; ++n,++to,ph+=2) *to = *ph;
    }
  }
  if(vel) {
    vect*to = static_cast<vect*>(vel);
    if(coerce) {
      const Vect*ph = static_cast<const Vect*>(phases) + 1;
      for(unsigned n=0; n!=N; ++n,++to,ph+=2) *to = *ph;
    } else {
      const vect*ph = static_cast<const vect*>(phases) + 1;
      for(unsigned n=0; n!=N; ++n,++to,ph+=2) *to = *ph;
    }
  }
  if(coerce) falcON_DEL_A(static_cast<Vect*>(phases));
  else       falcON_DEL_A(static_cast<vect*>(phases));
  if(pos)
    if(vel)
      DebugInfo(5,"data_in::read_phases(): %d %s & %s read\n",N,PosTag,VelTag);
    else
      DebugInfo(5,"data_in::read_phases(): %d %s read\n",N,PosTag);
  else
    DebugInfo(5,"data_in::read_phases(): %d %s read\n",N,VelTag);
  NREAD += N;
}
//------------------------------------------------------------------------------
// class falcON::nemo_out
//------------------------------------------------------------------------------
void falcON::nemo_out::close() falcON_THROWING
{
  if(SNAP) {
    DebugInfo(4,"nemo_out::close(): closing open snap_out first\n");
    SNAP->~snap_out();
    SNAP = 0;
  }
  nemo_io::close();
}
//------------------------------------------------------------------------------
falcON::nemo_out& falcON::nemo_out::open(const char*file, bool app)
falcON_THROWING
{
  // NOTE: we MUST call nemo_io::open() even if file==0
  // If the trailing letter of file is '!' or '@', we temporarily remove it
  // (setting it to nil), then call nemo_io::open(), and finally add it back on.
  close();
  char*last = (file&&file[0])? const_cast<char*>(file)+std::strlen(file)-1 : 0;
  char lett = last? *last : 0;            // remember last letter in filename
  if     (lett=='!') { *last=0; nemo_io::open(file,app? "a!":"w!"); }
  else if(lett=='@') { *last=0; nemo_io::open(file, "a"); }
  else                          nemo_io::open(file,app? "a":"w");
  if(last) *last = lett;                  // ensure file name remains unchanged
  return *this;
}
//------------------------------------------------------------------------------
// class falcON::snap_out
//------------------------------------------------------------------------------
falcON::snap_out::snap_out(nemo_out const&out, const unsigned nbod[BT_NUM],
			   double time) falcON_THROWING :
  OUTPUT(out), DATA(0), FIELDS_WRITTEN(0), NTOT(0u)
{
  DebugInfo(4,"snap_out::snap_out() ...\n");
  // 0 set # bodies
  for(bodytype t; t; ++t) NTOT += NBOD[t] = nbod[t];
  if(OUTPUT.SNAP)
    falcON_THROW("cannot open 2nd snapshot from nemo output stream");
  // 1 open snapshot set
  put_set(OUTPUT.STREAM,SnapShotTag);
  OUTPUT.SNAP = this;
  DebugInfo(5,"  snapshot opened\n");
  // 2 write parameter set
  put_set      (OUTPUT.STREAM,ParametersTag);
  Aux::put_data(OUTPUT.STREAM,NobjTag, IntType,&NTOT,0);
  Aux::put_data(OUTPUT.STREAM,NGasTag, IntType,&(NBOD[bodytype::gas]),0);
  Aux::put_data(OUTPUT.STREAM,NSinkTag,IntType,&(NBOD[bodytype::sink]),0);
  Aux::put_data(OUTPUT.STREAM,TimeTag ,DoubleType,&time,0);
  put_tes      (OUTPUT.STREAM,ParametersTag);
  DebugInfo(5,"  snap_out::snap_out(): parameter written:"
	     " Nbod=%d, Nsph=%d, Nsink=%d, time=%f\n",
	     NTOT, NBOD[bodytype::gas], NBOD[bodytype::sink], time);
  // 3 open particle set
  put_set (OUTPUT.STREAM,ParticlesTag);
  int CS = CSCode(Cartesian,Ndim,2);
  Aux::put_data(OUTPUT.STREAM,CoordSystemTag,IntType,&CS,0);
}
//------------------------------------------------------------------------------
falcON::snap_out::~snap_out() falcON_THROWING
{
  if(DATA) {
    DebugInfo(4,"snap_out::~snap_out(): closing open data_out first\n");
    DATA->~data_out();
  }
  NTOT = 0;
  for(bodytype t; t; ++t) NBOD[t] = 0u;
  put_tes(OUTPUT.STREAM,ParticlesTag);
  put_tes(OUTPUT.STREAM,SnapShotTag);
  OUTPUT.SNAP = 0;
  DebugInfo(4,"snap_out closed\n");
}
//------------------------------------------------------------------------------
// class falcON::data_out
//------------------------------------------------------------------------------
falcON::data_out::data_out(snap_out const&snap, nemo_io::Field f)
falcON_THROWING
: OUTPUT(snap), FIELD(f), NWRITTEN(0),
  NTOT(OUTPUT.N(FIELD)), TYPE(nemo_io::type(FIELD)),
  SUBN(::is_scalar(FIELD)? 1: ::is_vector(FIELD)? NDIM : 2*NDIM)
{
  DebugInfo(5,"data_out::data_out(%s) ...\n",NemoTag(FIELD));
  if( OUTPUT.DATA )
    falcON_THROW("cannot write %s: nemo output still engaged",
		 NemoTag(FIELD));
  if( OUTPUT.has_been_written(FIELD) )
    falcON_THROW("cannot write %s: has already been written",
		 NemoTag(FIELD));
  if(::is_scalar(FIELD)) {
    put_data_set(OUTPUT.stream(),NemoTag(FIELD),NemoType(TYPE),NTOT);
    DebugInfo(6,"  opening data set for %d scalars\n",NTOT);
  } else if(::is_vector(FIELD)) {
    put_data_set(OUTPUT.stream(),NemoTag(FIELD),NemoType(TYPE),NTOT,NDIM);
    DebugInfo(6,"  opening data set for %d vectors\n",NTOT);
  } else {
    put_data_set(OUTPUT.stream(),NemoTag(FIELD),NemoType(TYPE),
		 NTOT,2,NDIM);
    DebugInfo(6,"  opening data set for %d phases\n",NTOT);
  }
  OUTPUT.DATA = this;
}
//------------------------------------------------------------------------------
falcON::data_out::~data_out()
{
  if(NWRITTEN != NTOT)
    falcON_Warning("nemo output of %s: assigned %d, written only %d bodies\n",
		   NemoTag(FIELD), NTOT, NWRITTEN);
  put_data_tes(OUTPUT.stream(),NemoTag(FIELD));
  OUTPUT.DATA = 0;
  OUTPUT.FIELDS_WRITTEN |= FIELD;
  DebugInfo(5,"data_out(%s) closed\n",NemoTag(FIELD));
}
//------------------------------------------------------------------------------
void falcON::data_out::write(const void*data, unsigned n)
{
  if(NWRITTEN + n > NTOT) {
    falcON_Warning("nemo output of %s: "
		   "cannot write %d, only %d free spaces left\n",
		   NemoTag(FIELD), n, NTOT-NWRITTEN);
    n = NTOT - NWRITTEN;
  }
  put_data_blocked(OUTPUT.stream(),NemoTag(FIELD),data, n*SUBN);
  DebugInfo(6,"  %d %s written\n",n,NemoTag(FIELD));
  NWRITTEN += n;
}
//------------------------------------------------------------------------------
void falcON::data_out::write(const void*data)
{
  if(NWRITTEN < NTOT) {
    unsigned n = NTOT - NWRITTEN;
    put_data_blocked(OUTPUT.stream(),NemoTag(FIELD),data, n*SUBN);
    DebugInfo(6,"  %d %s written\n",n,NemoTag(FIELD));
    NWRITTEN += n;
  }
}
//------------------------------------------------------------------------------
namespace {
  int within_count = 0;
  inline bool Within(double const&val, char*range, falcON::real const&fuzz)
    // almost identical to nemo::within
  {
    char*endptr, *subptr, *sepptr, *colptr;
    falcON::real sublow, subhi;
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
//------------------------------------------------------------------------------
bool falcON::time_in_range(double t, const char*times)
{
  return  times == 0
    ||    std::strcmp(const_cast<char*>(times),"all") == 0
    ||    Within(t,const_cast<char*>(times),0.0005);
}

////////////////////////////////////////////////////////////////////////////////

#endif // falcON_NEMO
