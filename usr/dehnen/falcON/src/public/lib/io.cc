// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// io.cc                                                                       |
//                                                                             |
// Copyright (C) 2000-2007  Walter Dehnen                                      |
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
#include <public/io.h>
#include <fstream>
#include <cstdlib>
#include <cstring>

#ifdef   falcON_NEMO
// nemo header files
extern "C" {
# include <stdinc.h>
# include <filestruct.h>
# include <history.h>
# include <snapshot/snapshot.h>
}
#ifndef SoundSpeedTag
#  define SoundSpeedTag "SoundSpeed"
#endif
#ifndef TimeStepTag
#  define TimeStepTag "TimeStep"
#endif
#ifndef AuxVectorTag
#  define AuxVectorTag "AuxiliaryVector"
#endif
#ifndef LevelTag
#  define LevelTag "Level"
#endif
#ifndef JerkTag
#  define JerkTag "Jerk"
#endif
#ifndef GasHdotTag
#  define GasHdotTag "Hdot"
#endif
#ifndef GasFactTag
#  define GasFactTag "SPHFactor"
#endif
#ifndef MolWeightTag
#  define MolWeightTag "MolecularWeight"
#endif
#ifndef GasDensTag
#  warning
#  warning GasDensTag not #defined by NEMO
#  warning    please update to a more recent version of NEMO
#  warning
#  define GasDensTag "GasDensity"
#endif
#ifndef NumberTag
#  warning
#  warning NumberTag not #defined by NEMO
#  warning    please update to a more recent version of NEMO
#  warning
#  define NumberTag "NPartners"
#endif
#ifndef SPHNumberTag
#  warning
#  warning SPHNumberTag not #defined by NEMO
#  warning    please update to a more recent version of NEMO
#  warning
#  define SPHNumberTag "NSPHPartners"
#endif
#ifndef DensityTag
#  warning
#  warning DensityTag not #defined by NEMO
#  warning    please update to a more recent version of NEMO
#  warning
#  define DensityTag "Density"
#endif

#endif // falcON_NEMO

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
namespace {
  //----------------------------------------------------------------------------
  int openstdout = 0;
  int openstdin  = 0;
  //----------------------------------------------------------------------------
  inline void open_stdout() falcON_THROWING
  {
    if( ++openstdout > 1 )
      falcON_THROW("trying to open more than one output to stdout");
  }
  //----------------------------------------------------------------------------
  inline void close_stdout()
  {
    if(openstdout) --openstdout;
  }
  //----------------------------------------------------------------------------
  inline void open_stdin() falcON_THROWING
  {
    if( ++openstdin > 1 )
      falcON_THROW("trying to open more than one input from stdin");
  }
  //----------------------------------------------------------------------------
  inline void close_stdin()
  {
    if(openstdin) --openstdin;
  }
#ifdef   falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
  inline char* NemoTag(nemo_io::Field f) falcON_THROWING
  {
    switch(f) {
    case nemo_io::mass   : return MassTag;
    case nemo_io::pos    : return PosTag;
    case nemo_io::vel    : return VelTag;
    case nemo_io::eps    : return EpsTag;
    case nemo_io::key    : return KeyTag;
    case nemo_io::tau    : return TimeStepTag;
    case nemo_io::pot    : return PotentialTag;
    case nemo_io::acc    : return AccelerationTag;
    case nemo_io::jerk   : return JerkTag;
    case nemo_io::dens   : return DensityTag;
    case nemo_io::aux    : return AuxTag;
    case nemo_io::zet    : return AuxVectorTag;
    case nemo_io::lev    : return LevelTag;
    case nemo_io::num    : return NumberTag;
    case nemo_io::posvel : return PhaseSpaceTag;
    case nemo_io::SPHh   : return SmoothTag;
    case nemo_io::SPHnum : return SPHNumberTag;
    case nemo_io::SPHu   : return UinternTag;
    case nemo_io::SPHudot: return UdotIntTag;
    case nemo_io::SPHurad: return UdotRadTag;
    case nemo_io::SPHentr: return EntFuncTag;
    case nemo_io::SPHdens: return GasDensTag;
    case nemo_io::SPHhdot: return GasHdotTag;
    case nemo_io::SPHfact: return GasFactTag;
    case nemo_io::SPHcs  : return SoundSpeedTag;
    case nemo_io::SPHmu  : return MolWeightTag;
    case nemo_io::null:
      falcON_THROW("nemo I/O: nemo_io::null not I/O able");
    default:
      falcON_THROW("nemo I/O: unknown nemo_io::Field");
    }
  }
  //----------------------------------------------------------------------------
  inline nemo_io::DataType Type(const char* NemoType)
  {
    if(0 == std::strcmp(NemoType, ShortType ) ) return nemo_io::Short;
    if(0 == std::strcmp(NemoType, IntType   ) ) return nemo_io::Integer;
    if(0 == std::strcmp(NemoType, LongType  ) ) return nemo_io::Long;
    if(0 == std::strcmp(NemoType, FloatType ) ) return nemo_io::Single;
    if(0 == std::strcmp(NemoType, DoubleType) ) return nemo_io::Double;
    return nemo_io::Null;
  }
  //----------------------------------------------------------------------------
  inline char* NemoType(nemo_io::DataType t)
  {
    switch(t) {
    case nemo_io::Short  : return ShortType;
    case nemo_io::Integer: return IntType;
    case nemo_io::Long   : return LongType;
    case nemo_io::Single : return FloatType;
    case nemo_io::Double : return DoubleType;
    default              : return AnyType;
    }
  }
  //----------------------------------------------------------------------------
  inline bool is_scalar(nemo_io::Field f) {
    switch(f) {
    case nemo_io::mass:
    case nemo_io::eps:
    case nemo_io::key:
    case nemo_io::tau:
    case nemo_io::pot:
    case nemo_io::dens:
    case nemo_io::aux:
    case nemo_io::lev:
    case nemo_io::num:
    case nemo_io::SPHh:
    case nemo_io::SPHnum:
    case nemo_io::SPHu:
    case nemo_io::SPHudot:
    case nemo_io::SPHurad:
    case nemo_io::SPHentr:
    case nemo_io::SPHdens:
    case nemo_io::SPHhdot:
    case nemo_io::SPHfact:
    case nemo_io::SPHcs:  
    case nemo_io::SPHmu:   return true;
    default:               return false;
    }
  }
  //----------------------------------------------------------------------------
  inline bool is_vector(nemo_io::Field f) {
    switch(f) {
    case nemo_io::pos:
    case nemo_io::vel:
    case nemo_io::acc:
    case nemo_io::zet:
    case nemo_io::jerk: return true;
    default:            return false;
    }
  }
  //----------------------------------------------------------------------------
  inline bool is_phases(nemo_io::Field f) {
    return f == nemo_io::posvel;
  }
#endif // falcON_NEMO
  //////////////////////////////////////////////////////////////////////////////
} // namespace {
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::output                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void output::__open(bool append)
{
  APPENDING = false;
  if     (0 == FILE    ||
	  0 == FILE[0] ||
	  0 == std::strcmp(FILE,".") ) {
    OUT = 0;
    debug_info(2,"output: open sink\n");
  } else if(0 == std::strcmp(FILE,"-") ) {
    open_stdout();
    OUT = &std::cout;
    debug_info(2,"output: open stdout\n");
  } else {
    std::ofstream *FOUT = new std::ofstream();
    if(append) {
      FOUT->open(FILE,std::ios::out | std::ios::app);
      if(FOUT->is_open()) {
	APPENDING = true;
	debug_info(2,"output: append to file \"%s\"\n",FILE);
      }
    }
    if(!FOUT->is_open() )
      FOUT->open(FILE,std::ios::out);
    if( FOUT->is_open() ) {
      OUT = FOUT;
      debug_info(2,"output: open file \"%s\"\n",FILE);
    } else {
      debug_info(2,"output: could not open file \"%s\"\n",FILE);
      OUT = 0;
      falcON_DEL_O(FOUT);
    }
  }
}
//------------------------------------------------------------------------------
void output::close() {
  if(FREC) {
    if(FILE)
      warning("closing FortranORec before output from file \"%s\"\n",FILE);
    else
      warning("closing FortranORec before output\n");
    FREC->close();
  }
  debug_info(2,"output: closing\n");
  if(OUT == &std::cout) close_stdout();
  else if(OUT) falcON_DEL_O(OUT);
  APPENDING = false;
  OUT = 0;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::input                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void input::__open() {
  if     (0 == FILE || FILE[0] == 0) {
    IN = 0;
    debug_info(2,"input: empty file\n");
  } else if(0 == std::strcmp(FILE,"-") ) {
    open_stdin();
    IN= &std::cin;
    debug_info(2,"input: stdin\n");
  } else {
    std::ifstream *FIN = new std::ifstream(FILE);
    if( FIN->is_open() ) {
      IN = FIN;
      debug_info(2,"input: open file \"%s\"\n",FILE);
    } else {
      debug_info(2,"input: could not open file \"%s\"\n",FILE);
      IN = 0;
      falcON_DEL_O(FIN);
    }
  }
}
//------------------------------------------------------------------------------
void input::close() {
  if(FREC) {
    if(FILE)
      warning("closing FortranIRec before input from file \"%s\"\n",FILE);
    else
      warning("closing FortranIRec before input\n");
    FREC->close();
  }
  debug_info(2,"input: closing\n");
  if(IN == &std::cin) close_stdin();
  else if(IN) falcON_DEL_O(IN);
  IN = 0;
}
#ifdef   falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::nemo_io                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
nemo_io&nemo_io::open(const char* file, const char* control)
{
  if(STREAM) strclose(static_cast<stream>(STREAM));
  STREAM = static_cast<void*>(stropen(const_cast<char*>(file),
				      const_cast<char*>(control)));
  if     (0 == std::strcmp(const_cast<char*>(control), "r"))
    get_history(static_cast<stream>(STREAM));
  else if(0 == std::strcmp(const_cast<char*>(control), "w")  ||
	  0 == std::strcmp(const_cast<char*>(control), "w!") ||
	  0 == std::strcmp(const_cast<char*>(control), "s") )
    put_history(static_cast<stream>(STREAM));
  return *this;
}
//------------------------------------------------------------------------------
void nemo_io::close()
{
  if(STREAM) strclose(static_cast<stream>(STREAM));
  STREAM = 0;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::nemo_in                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void nemo_in::close() falcON_THROWING
{
  if(STREAM == 0) return;
  if(SNAP_IN) {
    debug_info(4,"nemo_in::close(): closing open snap_in first\n");
    SNAP_IN->~snap_in();
  }
  if(IS_PIPE) {
    close_stdin();
    IS_PIPE = false;
  }
  nemo_io::close();
  debug_info(4,"nemo_in: closed stream\n");
}
//------------------------------------------------------------------------------
nemo_in& nemo_in::open(const char* file) falcON_THROWING
{
  close();
  if(file == 0 || file[0] == 0) return *this;
  IS_PIPE = !std::strcmp(const_cast<char*>(file),"-");
  if(IS_PIPE) open_stdin();
  SNAP_IN = 0;
  nemo_io::open(file,"r");
  debug_info(4,"nemo_in: opened file '%s'\n",file);
  return *this;
}
//------------------------------------------------------------------------------
nemo_in::nemo_in(const char* file, const char* mode) :
  IS_PIPE(0), SNAP_IN(0)
{
  if(file == 0 || file[0] == 0) return;
  IS_PIPE = !std::strcmp(const_cast<char*>(file),"-");
  if(IS_PIPE) open_stdin();
  nemo_io::open(file,mode);
  debug_info(4,"nemo_in: opened file '%s'\n",file);
}
//------------------------------------------------------------------------------
nemo_in::~nemo_in() falcON_THROWING
{
  if(IS_PIPE) close_stdin();
  if(SNAP_IN) {
    debug_info(4,"nemo_in::~nemo_in(): closing open snap_in first\n");
    SNAP_IN->~snap_in();
  }
  debug_info(4,"nemo_in: closed stream\n");
}
//------------------------------------------------------------------------------
bool nemo_in::has_snapshot() const
{
  return STREAM && get_tag_ok(static_cast< ::stream >(STREAM),SnapShotTag);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::snap_in                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
snap_in::snap_in(nemo_in const&in) falcON_THROWING : 
  INPUT(in), DATA_IN(0), FIELDS_READ(0),
  HAS_NSPH(0), HAS_TIME(0), NBOD(0), NSPH(0), TIME(0.)
{
  debug_info(4,"snap_in::snap_in() ...\n");
  if(! INPUT.has_snapshot())
    falcON_THROW("cannot open snapshot from nemo input stream");
  if(INPUT.SNAP_IN)
    falcON_THROW("trying to open 2nd snapshot from nemo input stream");
  // 1 open snapshot set
  get_set(static_cast< ::stream >(INPUT.stream()),SnapShotTag);
  INPUT.SNAP_IN = this;
  debug_info(5,"  snapshot opened\n");
  // 2 open parameter set
  if(!get_tag_ok(static_cast< ::stream >(INPUT.stream()),ParametersTag)) {
    get_tes(static_cast< ::stream >(INPUT.stream()),SnapShotTag);
    INPUT.SNAP_IN = 0;
    falcON_THROW("cannot read parameters from nemo input stream");
  }
  get_set(static_cast< ::stream >(INPUT.stream()),ParametersTag);
  debug_info(5,"  parameter set opened\n");
  // 3 read parameter set
  // 3.1 read total # bodies 
  if(!get_tag_ok(static_cast< ::stream >(INPUT.stream()),NobjTag)) {
    get_tes(static_cast< ::stream >(INPUT.stream()),ParametersTag);
    get_tes(static_cast< ::stream >(INPUT.stream()),SnapShotTag);
    INPUT.SNAP_IN = 0;
    falcON_THROW("cannot read # bodies from nemo input stream");
  }
  get_data(static_cast< ::stream >(INPUT.stream()),NobjTag,IntType,&NBOD,0);
  debug_info(5,"  read Nobj = %d\n",NBOD);
  // 3.2 try to read # SPH bodies
  if(get_tag_ok(static_cast< ::stream >(INPUT.stream()),NGasTag)) {
    get_data(static_cast< ::stream >(INPUT.stream()),NGasTag,IntType,&NSPH,0);
    HAS_NSPH = true;
    debug_info(5,"  read Nsph = %d\n",NSPH);
    if(NSPH > NBOD)
      falcON_THROW("read nemo data: more SPH bodies than total");
  }
  // 3.3 try to read simulation time
  if(get_tag_ok(static_cast< ::stream >(INPUT.stream()),TimeTag)) {
    HAS_TIME = true;
    char* time_type = get_type(static_cast< ::stream >(INPUT.stream()),TimeTag);
    if(0 == std::strcmp(time_type, DoubleType))
      get_data(static_cast< ::stream >(INPUT.stream()),
	       TimeTag,DoubleType,&TIME,0);
    else if(0 == std::strcmp(time_type, FloatType)) {
      float __TIME;
      get_data(static_cast< ::stream >(INPUT.stream()),
	       TimeTag,FloatType,&__TIME,0);
      TIME = __TIME;
    } else
      warning("nemo input: unknown type '%s' for time\n",time_type);
  }
  if(HAS_TIME)
    debug_info(5,"  read time = %f\n",TIME);
  // 4 close parameter set
  get_tes(static_cast< ::stream >(INPUT.stream()),ParametersTag);
  debug_info(5,"  parameter set read & closed\n");
  // 5 open particle set
  if(!get_tag_ok(static_cast< ::stream >(INPUT.stream()),ParticlesTag)) {
    get_tes(static_cast< ::stream >(INPUT.stream()),SnapShotTag);
    INPUT.SNAP_IN = 0;
    falcON_THROW("cannot read parameters from nemo input stream");
  }
  get_set(static_cast< ::stream >(INPUT.stream()),ParticlesTag);
  debug_info(5,"  particles set opened\n");
}
//------------------------------------------------------------------------------
snap_in::~snap_in() falcON_THROWING
{
  if(DATA_IN) {
    debug_info(4,"snap_in::~snap_in(): closing open data_in first\n");
    DATA_IN->~data_in();
  }
  HAS_NSPH = false;
  HAS_TIME = false;
  get_tes(static_cast< ::stream >(INPUT.stream()),ParticlesTag);
  get_tes(static_cast< ::stream >(INPUT.stream()),SnapShotTag);
  INPUT.SNAP_IN = 0;
  debug_info(4,"snap_in closed\n");
}
//------------------------------------------------------------------------------
bool snap_in::has(nemo_io::Field f) const
{
  return 
    !has_been_read(f) &&
    get_tag_ok(static_cast< ::stream >(INPUT.stream()),NemoTag(f));
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::data_in                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
data_in::data_in(snap_in const&snap, nemo_io::Field f) falcON_THROWING :
  FIELD(f), INPUT(snap), NREAD(0), NTOT(0), TYPE(nemo_io::Null), SUBN(0)
{
  debug_info(5,"data_in::data_in(%s) ...\n",NemoTag(FIELD));
  if( INPUT.DATA_IN )
    falcON_THROW("cannot read %s: nemo input still engaged",NemoTag(FIELD));
  if( !INPUT.has(f) )
    falcON_THROW("cannot read %s: not given with nemo input",NemoTag(FIELD));
  if( INPUT.has_been_read(f) )
    falcON_THROW("cannot read %s: already read from nemo input",
		    NemoTag(FIELD));
  // 1 get type of data on file and check for mismatch
  char*TypeTag=get_type(static_cast< ::stream >(INPUT.stream()),NemoTag(FIELD));
  TYPE = Type(TypeTag);
  if(nemo_io::mismatch(TYPE, nemo_io::type(FIELD)))
    falcON_THROW("cannot read %s: type mismatch (got %s, expect %s)",
		 NemoTag(FIELD), NemoType(TYPE), 
		 NemoType(nemo_io::type(FIELD)));
  debug_info(6,"  data type: %s\n",nemo_io::type_name(TYPE));
  // 2 get dimensiones of data on file
  int*dim=get_dims(static_cast< ::stream >(INPUT.stream()),NemoTag(FIELD));
  if(!dim)
    falcON_THROW("cannot read # %s data",NemoTag(FIELD));
  NTOT = dim[0];
  // 3 check for dimensions to match expectations AND open data set
  if(dim[0] != snap.N(FIELD) ) {
    // 3.1 dim[0] != snap_in::N: error out
    falcON_THROW("nemo input of %s: found %d data, expected %d",
		 NemoTag(FIELD), dim[0], snap.N(FIELD));
  } else if(dim[1] == 0) {
    // 3.2 array of scalars
    if(!is_scalar(FIELD) )
      falcON_THROW("nemo input of %s: found scalars",NemoTag(FIELD));
    debug_info(6,"  opening data set for %d scalars\n",NTOT);
    get_data_set(static_cast< ::stream >(INPUT.stream()),
		 NemoTag(FIELD),TypeTag,NTOT,0);
    SUBN = 1;
  } else if(dim[2] == 0) {
    // 3.2 array of vectors
    if(!is_vector(FIELD) )
      falcON_THROW("nemo input of %s: found vectors",NemoTag(FIELD));
    if(dim[1] != NDIM)
      falcON_THROW("nemo input of %s: Ndim mismatch",NemoTag(FIELD));
    debug_info(6,"  opening data set for %d vectors\n",NTOT);
    get_data_set(static_cast< ::stream >(INPUT.stream()),
		 NemoTag(FIELD),TypeTag,NTOT,NDIM,0);
    SUBN = NDIM;
  } else if(dim[3] == 0) {
    // 3.2 array of phases
    if(!is_phases(FIELD) )
      falcON_THROW("nemo input of %s: found phases",NemoTag(FIELD));
    if(dim[1] != 2 && dim[2] != NDIM)
      falcON_THROW("nemo input of %s: Ndim mismatch",NemoTag(FIELD));
    debug_info(6,"  opening data set for %d phases\n",NTOT);
    get_data_set(static_cast< ::stream >(INPUT.stream()),
		 NemoTag(FIELD),TypeTag,NTOT,2,NDIM,0);
    SUBN = 2*NDIM;
  } else {
    // 3.3 array of yet higher rank
    falcON_THROW("nemo input of %s: found high-rank data",NemoTag(FIELD));
  }
  INPUT.DATA_IN = this;
}
//------------------------------------------------------------------------------
data_in::~data_in()
{
  get_data_tes(static_cast< ::stream >(INPUT.stream()),NemoTag(FIELD));
  INPUT.DATA_IN      = 0;
  INPUT.FIELDS_READ |= FIELD;
  debug_info(5,"data_in(%s) closed\n",NemoTag(FIELD));
}
//------------------------------------------------------------------------------
void data_in::read(void*data, unsigned n)
{
  if(NREAD + n > NTOT) {
    warning("nemo input of %s: cannot read %d, only %d data left",
	    NemoTag(FIELD), n, NTOT-NREAD);
    n = NTOT - NREAD;
  }
  if(n) {
    get_data_blocked(static_cast< ::stream >(INPUT.stream()),
		     NemoTag(FIELD), data, n*SUBN);
    debug_info(5,"data_in::read(): %d %s read\n",n,NemoTag(FIELD));
    NREAD += n;
  }
}
//------------------------------------------------------------------------------
void data_in::read(void*data)
{
  if(NREAD < NTOT) {
    unsigned n = NTOT - NREAD;
    get_data_blocked(static_cast< ::stream >(INPUT.stream()),
		     NemoTag(FIELD), data, n*SUBN);
    debug_info(5,"data_in::read(): %d %s read\n",n,NemoTag(FIELD));
    NREAD += n;
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::nemo_out                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void nemo_out::close() falcON_THROWING
{
  if(STREAM == 0) return;
  if(SNAP_OUT) {
    debug_info(4,"nemo_out::close(): closing open snap_out first\n");
    SNAP_OUT->~snap_out();
  }
  if(IS_PIPE) {
    close_stdout();
    IS_PIPE = false;
  }
  IS_SINK = false;
  nemo_io::close();
  debug_info(4,"nemo_out: closed stream\n");
}
//------------------------------------------------------------------------------
namespace {
  // if the last character in name equals c then                                
  // - we copy name into copy, except for character c                           
  // - we return true                                                           
  // otherwise                                                                  
  // - we return false                                                          
  inline bool is_appended(const char*name, char c, char*copy)
  {
    char *end = strrchr(name,c);
    if(end && end[1] == 0) {
      strcpy(copy,name);
      copy[size_t(end-name)] = 0;
      return true;
    } else
      return false;
  }
}
//------------------------------------------------------------------------------
nemo_out& nemo_out::open(const char* file,
			 bool appending) falcON_THROWING
{
  close();
  if(file == 0 || file[0] == 0) return *this;
  IS_PIPE = !std::strcmp(const_cast<char*>(file),"-");
  IS_SINK = !std::strcmp(const_cast<char*>(file),".");
  if(IS_PIPE) open_stdout();
  SNAP_OUT = 0;
  char copy[1024];
  if(appending) {                         // appending anyway?
    if(is_appended(file,'!',copy)) {        // overwrite & append
      nemo_io::open(copy,"a!");
      debug_info(4,"nemo_out: opened file '%s' for appending with overwrite\n",
 		 copy);
    } else if(is_appended(file,'@',copy)) { // append 
      nemo_io::open(copy,"a");
      debug_info(4,"nemo_out: opened file '%s' for appending\n",copy);
    } else {                                // append
      nemo_io::open(file,"a");
      debug_info(4,"nemo_out: opened file '%s' for appending\n",file);
    }
  } else {                                // append only with trailing '@'
    if       (is_appended(file,'!',copy)) { // overwrite existing file 
      nemo_io::open(copy,"w!");
      debug_info(4,"nemo_out: opened file '%s' with overwrite\n",copy);
    } else if(is_appended(file,'@',copy)) { // append to existing file
      nemo_io::open(copy,"a");
      debug_info(4,"nemo_out: opened file '%s' for appending\n",copy);
    } else {                                // open new file
      nemo_io::open(file,"w");
      debug_info(4,"nemo_out: opened file '%s'\n",file);
    }
  }
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::snap_out                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
snap_out::snap_out(nemo_out const&out,
		   unsigned       nb,
		   unsigned       ns,
		   double         time) falcON_THROWING :
  OUTPUT(out), DATA_OUT(0), FIELDS_WRITTEN(0), NBOD(nb), NSPH(ns)
{
  debug_info(4,"snap_out::snap_out() ...\n");
  if(OUTPUT.SNAP_OUT)
    falcON_THROW("cannot open 2nd snapshot from nemo output stream");
  // 1 open snapshot set
  put_set(static_cast< ::stream >(OUTPUT.stream()),SnapShotTag);
  OUTPUT.SNAP_OUT = this;
  debug_info(5,"  snapshot opened\n");
  // 2 write parameter set
  put_set(static_cast< ::stream >(OUTPUT.stream()),ParametersTag);
  put_data(static_cast< ::stream >(OUTPUT.stream()),NobjTag,IntType,&NBOD,0);
  put_data(static_cast< ::stream >(OUTPUT.stream()),NGasTag,IntType,&NSPH,0);
  put_data(static_cast< ::stream >(OUTPUT.stream()),
	   TimeTag,DoubleType,&time,0);
  put_tes(static_cast< ::stream >(OUTPUT.stream()),ParametersTag);
  debug_info(5,"  parameter written: Nbod=%d, Nsph=%d, time=%f\n",
	     NBOD, NSPH, time);
  // 3 open particle set
  put_set(static_cast< ::stream >(OUTPUT.stream()),ParticlesTag);
  int CS = CSCode(Cartesian,Ndim,2);
  put_data(static_cast< ::stream>(OUTPUT.stream()),
	   CoordSystemTag,IntType,&CS,0);
}
//------------------------------------------------------------------------------
snap_out::~snap_out() falcON_THROWING
{
  if(DATA_OUT) {
    debug_info(4,"snap_out::~snap_out(): closing open data_out first\n");
    DATA_OUT->~data_out();
  }
  NBOD = 0;
  NSPH = 0;
  put_tes(static_cast< ::stream >(OUTPUT.stream()),ParticlesTag);
  put_tes(static_cast< ::stream >(OUTPUT.stream()),SnapShotTag);
  OUTPUT.SNAP_OUT = 0;
  debug_info(4,"snap_out closed\n");
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::data_out                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
data_out::data_out(snap_out const&snap, nemo_io::Field f) falcON_THROWING
  : FIELD(f),
    OUTPUT(snap),
    NWRITTEN(0),
    NTOT (nemo_io::is_sph(FIELD)? OUTPUT.NSPH : OUTPUT.NBOD),
    TYPE (nemo_io::type(FIELD)),
    SUBN (is_scalar(FIELD)? 1:
	  is_vector(FIELD)? NDIM : 2*NDIM)
{
  debug_info(5,"data_out::data_out(%s) ...\n",NemoTag(FIELD));
  if( OUTPUT.DATA_OUT )
    falcON_THROW("cannot write %s: nemo output still engaged",
		 NemoTag(FIELD));
  if( OUTPUT.has_been_written(FIELD) )
    falcON_THROW("cannot write %s: has already been written",
		 NemoTag(FIELD));
  if(is_scalar(FIELD)) {
    put_data_set(static_cast< ::stream >(OUTPUT.stream()),
		 NemoTag(FIELD),NemoType(TYPE),NTOT,0);
    debug_info(6,"  opening data set for %d scalars\n",NTOT);
  } else if(is_vector(FIELD)) {
    put_data_set(static_cast< ::stream >(OUTPUT.stream()),
		 NemoTag(FIELD),NemoType(TYPE),NTOT,NDIM,0);
    debug_info(6,"  opening data set for %d vectors\n",NTOT);
  } else {
    put_data_set(static_cast< ::stream >(OUTPUT.stream()),
		 NemoTag(FIELD),NemoType(TYPE),NTOT,2,NDIM,0);
    debug_info(6,"  opening data set for %d phases\n",NTOT);
  }
  OUTPUT.DATA_OUT = this;
}
//------------------------------------------------------------------------------
data_out::~data_out()
{
  if(NWRITTEN != NTOT)
    warning("nemo output of %s: assigned %d, written only %d bodies\n",
	    NemoTag(FIELD), NTOT, NWRITTEN);
  put_data_tes(static_cast< ::stream >(OUTPUT.stream()),NemoTag(FIELD));
  OUTPUT.DATA_OUT        = 0;
  OUTPUT.FIELDS_WRITTEN |= FIELD;
  debug_info(5,"data_out(%s) closed\n",NemoTag(FIELD));
}
//------------------------------------------------------------------------------
void data_out::write(const void*data, unsigned n)
{
  if(NWRITTEN + n > NTOT) {
    warning("nemo output of %s: cannot write %d, only %d free spaces left",
	    NemoTag(FIELD), n, NTOT-NWRITTEN);
    n = NTOT - NWRITTEN;
  }
  put_data_blocked(static_cast< ::stream >(OUTPUT.stream()),
		   NemoTag(FIELD), const_cast<void*>(data), n*SUBN);
  debug_info(6,"  %d %s written\n",n,NemoTag(FIELD));
  NWRITTEN += n;
}
//------------------------------------------------------------------------------
void data_out::write(const void*data)
{
  if(NWRITTEN < NTOT) {
    unsigned n = NTOT - NWRITTEN;
    put_data_blocked(static_cast< ::stream >(OUTPUT.stream()),
		     NemoTag(FIELD), const_cast<void*>(data), n*SUBN);
    debug_info(6,"  %d %s written\n",n,NemoTag(FIELD));
    NWRITTEN += n;
  }
}
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  static int within_count = 0;
  inline bool within(double const&val, char*range, real const&fuzz)
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
//                                                                            //
// function falcON::time_in_range                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
bool falcON::time_in_range(double t, const char*times)
{
  return  times == 0
    ||    std::strcmp(const_cast<char*>(times),"all") == 0
    ||    falcON::within(t,const_cast<char*>(times),0.0005);
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::FortranIRec                                                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
unsigned falcON::FortranIRec::read_size() throw(falcON::exception)
{
  if(HSZE == 4) {
    uint32 S;
    IN.read(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint32));
    return S;
  } else if(HSZE == 8) {
    uint64 S;
    IN.read(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint64));
    return S;
  } else 
    throw exception("FortranIRec: header size must be 4 or 8\n");
  return 0;
}
//------------------------------------------------------------------------------
falcON::FortranIRec::FortranIRec(input& in, unsigned rec)
  throw(falcON::exception) : IN(in), HSZE(rec), READ(0)
{
  debug_info(8,"FortranIRec: opening ... \n");
  if(!IN) throw exception("FortranIRec::FortranIRec(): input corrupted");
  if(IN.FREC)
    throw exception("trying to open 2nd FortranIRec to same input\n");
  IN.FREC = this;
  SIZE = read_size();
  debug_info(6,"FortranIRec: opened with %u bytes\n",SIZE);
}
//------------------------------------------------------------------------------
unsigned falcON::FortranIRec::read_bytes(char*buf, unsigned n)
  throw(falcON::exception)
{
  if(!IN) throw exception("FortranIRec::read_bytes(): input corrupted");
  if(READ + n > SIZE) {
    warning("FortranIRec::read(): cannot read %u, but only %u bytes\n",
	    n, SIZE-READ);
    n = SIZE - READ;
  }
  IN.read(buf,n);
  if(!IN) throw exception("FortranIRec: input corrupted");
  READ += n;
  debug_info(6,"FortranIRec: read %u bytes\n",n);
  return n;
}
//------------------------------------------------------------------------------
void falcON::FortranIRec::skip_bytes(unsigned n) {
  if(READ + n > SIZE) n = SIZE - READ;
  if(n && !IN) throw exception("FortranIRec::skip_bytes(): input corrupted");
  for(char C; n; --n,++READ) IN.read(&C,1);
}
//------------------------------------------------------------------------------
void falcON::FortranIRec::close() throw(falcON::exception)
{
  if(!IN) throw exception("FortranIRec::close(): input corrupted");
  if(READ != SIZE) {
    warning("FortranIRec: only %u of %u bytes read on closing record\n",
	    READ, SIZE);
    for(char C; READ!=SIZE; ++READ) IN.read(&C,1);
  }
  unsigned S = read_size();
  IN.FREC = 0;
  if(S != SIZE) throw exception("FortranIRec: record size mismatch");
  debug_info(6,"FortranIRec: closed with %u bytes\n",SIZE);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::FortranORec                                                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void falcON::FortranORec::write_size() throw(falcON::exception)
{
  if(HSZE == 4) {
    uint32 S = SIZE;
    OUT.write(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint32));
  } else if(HSZE == 8) {
    uint64 S;
    OUT.write(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint64));
  } else 
    throw exception("FortranORec: header size must be 4 or 8\n");
}
//------------------------------------------------------------------------------
falcON::FortranORec::FortranORec(output& out, unsigned size, unsigned rec)
  throw(falcON::exception) : OUT(out), HSZE(rec), SIZE(size), WRITTEN(0)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(OUT.FREC)
    throw exception("trying to open 2nd FortranORec to same output\n");
  OUT.FREC = this;
  write_size();
  debug_info(6,"FortranORec: opened for %u bytes\n",SIZE);
}
//------------------------------------------------------------------------------
unsigned falcON::FortranORec::write_bytes(const char*buf, unsigned n)
  throw(falcON::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN + n > SIZE) {
    warning("FortranORec::write(): cannot read %u, but only %u bytes\n",
	    n, SIZE-WRITTEN);
    n = SIZE - WRITTEN;
  }
  OUT.write(buf,n);
  if(!OUT) throw exception("FortranORec: ostream corrupted");
  WRITTEN += n;
  debug_info(6,"FortranORec: written %u bytes\n",n);
  return n;
}
//------------------------------------------------------------------------------
void falcON::FortranORec::fill_bytes(unsigned n, char C) {
  if(WRITTEN + n > SIZE) n = SIZE - WRITTEN;
  for(; n; --n,++WRITTEN) OUT.write(&C,1);
}
//------------------------------------------------------------------------------
void falcON::FortranORec::close() throw(falcON::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN!=SIZE) {
    warning("FortranORec: only %u of %u bytes written on closing record"
	    " ... padding with 0\n", WRITTEN, SIZE);
    for(char C=0; WRITTEN!=SIZE; ++WRITTEN) OUT.write(&C,1);
  }
  write_size();
  OUT.FREC = 0;
  debug_info(6,"FortranORec: closed with %u bytes\n",SIZE);
}
////////////////////////////////////////////////////////////////////////////////
