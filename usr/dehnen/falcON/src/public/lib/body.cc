// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.cc                                                                     |
//                                                                             |
// Copyright (C) 2000-2008 Walter Dehnen                                       |
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
#include <body.h>                                  // falcON::bodies etc        
#include <iostream>                                // C++ basic I/O             
#include <fstream>                                 // C++ file I/O              
#include <sstream>                                 // C++ string I/O            
#include <iomanip>                                 // C++ I/O formating         
#include <cstring>                                 // C++ strings               
#include <public/io.h>                             // utilities for NEMO I/O    
#include <numerics.h>

#ifdef falcON_NEMO                                 // compiler option           
  extern "C" {
#   include <stdinc.h>                             // NEMO basics               
  }
#endif

using namespace falcON;

falcON_TRAITS(falcON::bodies::block,"bodies::block");
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// struct falcON::bodies::block                                                 
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void bodies::block::clone(block*that)
{
  if(that == this) return;
  DebugInfo(3,"bodies::block::clone(): cloning block with %d [%d] %s\n",
	    that->NBOD,that->NALL,that->TYPE.name());
  if(this->TYPE != that->TYPE)
    falcON_THROW("bodies::block::clone(): bodytype mismatch ('%s' vs '%s')\n",
		 this->TYPE.name(), that->TYPE.name());
  for(fieldbit f; f; ++f) {
    this->del_field(f);
    this->set_data_void(f,that->data_void(f));
    that->set_data_void(f,0);
  }
  this->NALL  = that->NALL;
  this->NBOD  = that->NBOD;
  this->FIRST = that->FIRST;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::reset_flags() const
{
  if(0 != DATA[fieldbit::f]) {
    if(TYPE.is_sph()) 
      for(int n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sph;
    else
    if(TYPE.is_sink()) 
      for(int n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sink;
    else
      for(int n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::empty;
  }
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::flag_all_as_active() const falcON_THROWING
{
  if(0 != DATA[fieldbit::f])
    for(int n=0; n!=NALL; ++n)
      datum<fieldbit::f>(n).add(flags::active);
  else 
    falcON_THROW("in bodies::flag_all_as_active(): flags not supported");
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::reset_data(fieldset b) const falcON_THROWING {
#define RESETDATA(BIT,NAME)			\
  if(DATA[BIT] && b.contain(BIT) && NBOD)	\
    for(int n=0; n!=NBOD; ++n)			\
      field_traits<BIT>::set_zero(datum<BIT>(n));
 DEF_NAMED(RESETDATA)
#undef RESETDATA
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::add_field (fieldbit f) falcON_THROWING {
  if(TYPE.allows(f) && 0 == DATA[value(f)] ) {
    if(debug(4)) DebugInfo("bodies::block::add_field(): "
			   "allocating data for %s bodies: %u %c (%s)\n",
			   TYPE.name(),NALL,letter(f),fullname(f));
    set_data_void(f, falcON_NEW(char,NALL*falcON::size(f)));
    if(f == fieldbit::f) reset_flags();
  }
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::del_field (fieldbit f) falcON_THROWING {
  if(DATA[value(f)]) {
    if(debug(4)) DebugInfo("bodies::block::del_field(): "
			   "de-allocating data for %s bodies: %c (%s)\n",
			   TYPE.name(),letter(f),fullname(f));
    falcON_DEL_A(static_cast<char*>(DATA[value(f)]));
  }
  set_data_void(f,0);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::swap_bytes(fieldbit f) falcON_THROWING {
  if(DATA[value(f)]) {
    if(debug(4)) DebugInfo("bodies::block::swap_bytes(): "
			   "swapping bytes for %c (%s)\n",
			   letter(f),fullname(f));
    falcON::swap_bytes(DATA[value(f)], falcON::size(f), NALL);
  }
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::add_fields(fieldset b) falcON_THROWING {
  for(fieldbit f; f; ++f)
    if(b.contain(f)) add_field(f);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::del_fields(fieldset b) falcON_THROWING {
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) del_field(f);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::set_fields(fieldset b) falcON_THROWING {
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) add_field(f);
    else             del_field(f);
}
////////////////////////////////////////////////////////////////////////////////
bodies::block::~block() falcON_THROWING {
  for(fieldbit f; f; ++f)
    del_field(f);
}
////////////////////////////////////////////////////////////////////////////////
bodies::block::block(unsigned no,                  // I: our No                 
		     unsigned na,                  // I: data to allocate       
		     unsigned nb,                  // I: # bodies <= na_b       
		     unsigned fst,                 // I: first body index       
		     bodytype type,                // I: hold sph bodies?       
		     fieldset bits,                // I: data to allocate       
		     bodies  *bods)                // I: pointer to my bodies   
  falcON_THROWING
  : TYPE ( type ),
    NALL ( na ), 
    NBOD ( nb ), 
    NO   ( no ),
    FIRST( fst ),
    NEXT ( 0 ),
    BODS ( bods )
{
  if(na<nb)
    falcON_THROW("in bodies::block::block(): N_alloc < N_bodies");
  if(debug(6))
    DebugInfo("bodies::block: na=%d, bits=%s, type=%s allowed bits=%s\n",
	      na, word(bits), TYPE.name(), word(bits&TYPE.allows()));
  bits &= TYPE.allows();
  for(fieldbit f; f; ++f)
    set_data_void(f,0);
  try {
    add_fields(bits);
  } catch(exception E) {
    del_fields(fieldset::all);
    for(fieldbit f; f; ++f)
      del_field(f);
    falcON_RETHROW(E);
  }
}
///////////////////////////////////////////////////////////////////////////////
template<unsigned BIT=0, unsigned END=BodyData::NQUANT> struct CopyBody {
  static const fieldset::value_type BD = fieldset::one <<BIT;
  /// copies data field BIT for a single body within the same block
  static void copy(void    **data,
		   unsigned  from,
		   unsigned  to  ,
		   fieldset  b,
		   fieldset &c) {
    if(data[BIT] && b.contain(fieldbit(BIT)) ) {
      memcpy(static_cast<      char*>(data[BIT])+to  *BodyData::ZQUANT[BIT],
	     static_cast<const char*>(data[BIT])+from*BodyData::ZQUANT[BIT],
	     BodyData::ZQUANT[BIT]);
      c |= fieldset(fieldbit(BIT));
    }
    CopyBody<BIT+1, END>::copy(data,from,to,b,c);
  }
};
template<unsigned BIT> struct CopyBody<BIT,BIT> {
  static void copy(void**, unsigned, unsigned, fieldset, fieldset& ) {}
};
//------------------------------------------------------------------------------
fieldset bodies::block::copy_body(unsigned from, unsigned to, fieldset b)
{
  fieldset copied;
  if(from != to)
    CopyBody<0>::copy(DATA,from,to,b,copied);
  return copied;
}
////////////////////////////////////////////////////////////////////////////////
template<unsigned BIT=0, unsigned END=BodyData::NQUANT> struct CopyBodies {
  static const fieldset::value_type BD = fieldset::one <<BIT;
  static void copy(void*const*data_fr,
		   void*const*data_to,
		   unsigned  fr,
		   unsigned  to,
		   unsigned  num,
		   fieldset  b,
		   fieldset&c) {
    if(data_fr[BIT] && data_to[BIT] && b & fieldset(fieldbit(BIT)) ) {
      memcpy(static_cast<      char*>(data_to[BIT])+to*BodyData::ZQUANT[BIT],
	     static_cast<const char*>(data_fr[BIT])+fr*BodyData::ZQUANT[BIT],
	     num*BodyData::ZQUANT[BIT]);
      c |= fieldset(fieldbit(BIT));
    }
    CopyBodies<BIT+1, END>::copy(data_fr,data_to,fr,to,num,b,c);
  }
};
template<unsigned BIT> struct CopyBodies<BIT,BIT> {
  static void copy(void*const*, void*const*, unsigned, unsigned, unsigned,
		   fieldset, fieldset&) {}
};
//------------------------------------------------------------------------------
fieldset bodies::block::copy_bodies(const block*other,
				    unsigned    from,
				    unsigned    to,
				    unsigned    num,
				    fieldset    copy) falcON_THROWING
{
  fieldset copied;
  if(this == other)
    falcON_THROW("in bodies::block::copy_bodies(): this == other");
  else
    CopyBodies<0>::copy(other->DATA,DATA,from,to,num,copy,copied);
  return copied;
}
////////////////////////////////////////////////////////////////////////////////
inline void bodies::block::skip(unsigned&from,
				flags    copyflag) const falcON_THROWING
{
  if(copyflag)
    for(; from<NBOD && !(flag(from).are_set(copyflag)); ++from );
}
//------------------------------------------------------------------------------
// copy up to NALL bodies                                                 
// - we copy only bodies of the same type as hold here                    
// - we only copy bodies whose flag matches last argument                 
// - if the block copied is finished, we take its NEXT, starting at i=0   
// - on return the block pointer is either NULL or together with i they   
//   give the first body which was not copied because NALL was exceeded   
// - NBOD is set to the bodies copied                                     
fieldset bodies::block::copy(const block*&From,
			     unsigned    &from,
			     fieldset     copydata,
			     flags        copyflag) falcON_THROWING
{
  if( From == this )
    falcON_THROW("in bodies::block::copy(): cannot copy from self");
  NBOD = 0u;
  if( From == 0) return fieldset::empty;
  unsigned copy;
  unsigned free = NALL;
  fieldset copied;
  if(copyflag && ! has_field(fieldbit::f) )
    falcON_THROW("in bodies::block::copy(): "
		 "copyflag!=0 but flags not supported");
  // skip bodies not to be copied                                               
  From->skip(from,copyflag);
  while(free &&                              // WHILE  we have still space      
	From &&                              //   AND  the copied block is valid
	From->TYPE == TYPE &&                //   AND  its type is ours         
	from < From->NBOD ) {                //   AND  the index is valid too   
    // determine number of bodies to be copied to position from                 
    if(copyflag) {
      copy = 0u;
      for(unsigned to=from;
	  to < From->NBOD && From->flag(to).are_set(copyflag) && copy<free;
	  ++copy, ++to);
    } else
      copy = min(free, From->NBOD - from);
    // if any body to be copied, copy data, adjust free, NBOD, from, copied     
    if(copy) {
      fieldset c = copy_bodies(From, from, NBOD, copy, copydata);
      free -= copy;
      NBOD += copy;
      from += copy;
      copied = copied? c : copied & c;
    }
    // skip bodies not to be copied                                             
    From->skip(from,copyflag);
    // end of input block? then take next block                                 
    if(from == From->NBOD) {
      From = From->NEXT;
      if(From == this) 
	falcON_THROW("in bodies::block::copy(): cannot copy from self");
      from = 0;
      if(From) From->skip(from,copyflag);
    }
  }
  return copied;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::remove(unsigned &removed) falcON_THROWING
{
  if(NBOD == 0) return;
  if(0 == DATA[fieldbit::f] )
    falcON_THROW("in bodies::remove(): flags needed but not supported");
  unsigned lo=0u, hi=NBOD-1;
  while(lo < hi) {
    while(! to_remove(const_datum<fieldbit::f>(lo)) && lo < hi) ++lo;
    while(  to_remove(const_datum<fieldbit::f>(hi)) && lo < hi) --hi;
    if(lo < hi) copy_body(hi--,lo++);
  }
  if(lo == hi && ! to_remove(const_datum<fieldbit::f>(lo))) ++lo;
  removed += NBOD - lo;
  NBOD     = lo;
}
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
void bodies::block::read_data(data_in &inpt,
			      unsigned from,
			      unsigned N) falcON_THROWING
{
  fieldbit f= nemo_io::bit(inpt.field());
  if(!TYPE.allows(f))
    falcON_THROW("bodies::block::read_data(%c): not allowed by our type",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::read_data(%c): cannot read that many",
		 letter(f));
  add_field(f);
  if(inpt.must_coerce()) {
    if(debug(1))
      DebugInfo("bodies::block::read_data(%c): must convert from %s to %s",
		letter(f),
		nemo_io::type_name(nemo_io::NotReal),
		nemo_io::type_name(nemo_io::Real));
    unsigned ntot = N*inpt.sub_N();
    notreal* data = falcON_NEW(notreal, ntot);
    inpt.read(data, N);
    const notreal* d = data;
    real         * D = static_cast<real*>(DATA[value(f)])+from;
    for(unsigned i=0; i!=ntot; ++i,++d,++D) *D = *d;
    falcON_DEL_A(data);
  } else
    inpt.read(static_cast<char*>(DATA[value(f)])+from*falcON::size(f), N);
}
////////////////////////////////////////////////////////////////////////////////
namespace {
  typedef tupel<Ndim,notreal> Vect;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::read_posvel(data_in &inpt,
				unsigned from,
				unsigned N,
				fieldset want) falcON_THROWING
{
  if(inpt.field() != nemo_io::posvel)
    falcON_THROW("bodies::block::read_posvel(): input has not phases");
  if(from + N > NBOD)
    falcON_THROW("bodies::block::read_posvel(): cannot read that many");
  const bool coerce = inpt.must_coerce();
  if(coerce && debug(1))
    DebugInfo("bodies::read_posvel(): must convert from %s to %s",
	      nemo_io::type_name(nemo_io::NotReal),
	      nemo_io::type_name(nemo_io::Real));
  void* phases = coerce?
    static_cast<void*>(falcON_NEW(Vect,2*N)) : 
    static_cast<void*>(falcON_NEW(vect,2*N)) ;
  inpt.read(phases, N);
  if(want.contain(fieldbit::x)) {
    add_field(fieldbit::x);
    vect*to = static_cast<vect*>(DATA[fieldbit::x]) + from;
    if(coerce) {
      const Vect*ph = static_cast<const Vect*>(phases);
      for(unsigned n=0; n!=N; ++n,++to,ph+=2) *to = *ph;
    } else {
      const vect*ph = static_cast<const vect*>(phases);
      for(unsigned n=0; n!=N; ++n,++to,ph+=2) *to = *ph;
    }
  }
  if(want.contain(fieldbit::v)) {
    add_field(fieldbit::v);
    vect*to = static_cast<vect*>(DATA[fieldbit::v]) + from;
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
  if(debug(2)) DebugInfo("bodies::block::read_posvel(): read %s",
			 word(want&fieldset::phases));
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::write_data(data_out&outp,
			       unsigned from,
			       unsigned N) const falcON_THROWING
{
  fieldbit f= nemo_io::bit(outp.field());
  if(0 == DATA[value(f)])
    falcON_THROW("bodies::block::write_data(%c): data not supported",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::write_data(%c): cannot write that many",
		 letter(f));
  outp.write(static_cast<char*>(DATA[value(f)])+from*falcON::size(f), N);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::write_potpex(data_out&outp,
				 unsigned from,
				 unsigned N) const falcON_THROWING
{
  if(outp.field() != nemo_io::pot)
    falcON_THROW("bodies::block::write_potpex(): wrong field");
  if(0==DATA[fieldbit::p] || 0==DATA[fieldbit::q])
    falcON_THROW("bodies::block::write_potpex(): data not supported");
  if(from + N > NBOD)
    falcON_THROW("bodies::block::write_potpex(): cannot write that many");
  real *P = falcON_NEW(real,N);
  for(int n=0,m=from; n!=N; ++n,++m)
    P[n] = const_datum<fieldbit::p>(m) + const_datum<fieldbit::q>(m);
  outp.write(P,N);
  falcON_DEL_A(P);
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_REAL_IS_FLOAT
void bodies::block::read_Fortran(FortranIRec&I, fieldbit f, unsigned from,
				 unsigned N, bool swap) falcON_THROWING
{
  if(!TYPE.allows(f))
    falcON_THROW("bodies::block::read_Fortran(%c): not allowed by our type",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::read_Fortran(%c): cannot read that many",
		 letter(f));
  add_field(f);
  char    *C = static_cast<char*>(DATA[value(f)])+from*falcON::size(f);
  unsigned R = I.read_bytes(C, N*falcON::size(f));
  if(swap) 
    if(is_vector(f))
      falcON::swap_bytes(static_cast<void*>(C), sizeof(real), Ndim*N);
    else
      falcON::swap_bytes(static_cast<void*>(C), falcON::size(f), N);

  if(R != N*falcON::size(f))
    falcON_THROW("bodies::block::read_Fortran(%c): "
		 "could only read %u of %u bytes\n",R,N*falcON::size(f));
  if(debug(4))
    DebugInfo("bodies::block::read_Fortran(): read %u `%s'\n",N,fullname(f));
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::write_Fortran(FortranORec&O, fieldbit f, unsigned from,
				  unsigned N) const falcON_THROWING
{
  if(0 == DATA[value(f)])
    falcON_THROW("bodies::block::write_Fortran(%c): data not supported",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::write_Fortran(%c): cannot write that many",
		 letter(f));
  unsigned W = O.write_bytes(static_cast<const char*>(DATA[value(f)])
			     +from*falcON::size(f), N*falcON::size(f));
  if(W != N*falcON::size(f))
    falcON_THROW("bodies::block::write_Fortran(%c): "
		 "could only write %u of %u bytes\n",W,N*falcON::size(f));
  if(debug(4))
    DebugInfo("bodies::block::write_Fortran(): written %u `%s'\n",
	      N,fullname(f));
}
#endif
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::bodies::iterator                                               
//                                                                              
////////////////////////////////////////////////////////////////////////////////
bodies::iterator& bodies::iterator::read_data(data_in&D, unsigned R)
  falcON_THROWING
{
  if(R == 0 || R > D.N_unread()) R = D.N_unread();
  while(is_valid() && R) {
    unsigned r = min(N-K, R);
    const_cast<block*>(B)->read_data(D,K,r);
    R -= r;
    K += r;
    if(K >= N) next_block();
  }
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
bodies::iterator& bodies::iterator::write_data(data_out&D, unsigned W)
  falcON_THROWING
{
  if(W == 0 || W > D.N_free()) W = D.N_free();
  while(is_valid() && W) {
    unsigned w = min(N-K, W);
    B->write_data(D,K,w);
    W -= w;
    K += w;
    if(K >= N) next_block();
  }
}
////////////////////////////////////////////////////////////////////////////////
bodies::iterator& bodies::iterator::write_potpex(data_out&D, unsigned W)
  falcON_THROWING
{
  if(W == 0 || W > D.N_free()) W = D.N_free();
  while(is_valid() && W) {
    unsigned w = min(N-K, W);
    B->write_potpex(D,K,w);
    W -= w;
    K += w;
    if(K >= N) next_block();
  }
}
////////////////////////////////////////////////////////////////////////////////
bodies::iterator& bodies::iterator::read_posvel(data_in& D, fieldset get,
						unsigned R)
  falcON_THROWING
{
  if(R == 0 || R > D.N_unread()) R = D.N_unread();
  while(is_valid() && R) {
    unsigned r = min(N-K, D.N_unread());
    const_cast<block*>(B)->read_posvel(D,K,r,get);
    R -= r;
    K += r;
    if(K >= N) next_block();
  }
  return *this;
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_REAL_IS_FLOAT
bodies::iterator& bodies::iterator::read_Fortran(FortranIRec&I, fieldbit f, 
						 unsigned R, bool swap)
  falcON_THROWING
{
  if(R * falcON::size(f) > I.bytes_unread())
    falcON_THROW("body::read_Fortran: want %u `%s' (%u bytes) but "
		 "only %u bytes left on Fortran record\n",
		 R, fullname(f), R*falcON::size(f), I.bytes_unread());
  while(is_valid() && R) {
    unsigned r = min(N-K, R);
    const_cast<block*>(B)->read_Fortran(I,f,K,r,swap);
    R -= r;
    K += r;
    if(K >= N) next_block();
  }
  if(R) falcON_THROW("body::read_Fortran: %u data remain unread\n",R);
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
bodies::iterator& bodies::iterator::write_Fortran(FortranORec&O,
						  fieldbit f, unsigned W)
  falcON_THROWING
{
  if(W * falcON::size(f) > O.bytes_free())
    falcON_THROW("body::write_Fortran: want %u `%s' (%u bytes) but "
		 "only %u bytes left free on Fortran record\n",
		 W, fullname(f), W*falcON::size(f), O.bytes_free());
  while(is_valid() && W) {
    unsigned w = min(N-K, W);
    B->write_Fortran(O,f,K,w);
    W -= w;
    K += w;
    if(K >= N) next_block();
  }
  if(W) falcON_THROW("body::write_Fortran: %u data remain unwritten\n",W);
  return *this;
}
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::bodies                                                         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
// replace a given block
#if(0)
void bodies::replace_block(block*Bold, block*Bnew)
{
  // 1  get entry in BLOCK[]
  if(Bold == 0) falcON_THROW("bodies::replace_block(): Bold=0\n");
  bool found = false;
  for(int p=0; p!=index::max_blocks; ++p)
    if(BLOCK[p] == Bold) {
      BLOCK[p] = Bnew;
      found = true;
      break;
    }
  if(!found) falcON_THROW("bodies::replace_block(): Bold (%p) not found\n",
			  Bold);
  // 2  get block::NEXT correct
  Bnew->NEXT = Bold->NEXT;
  for(int p=0; p!=index::max_blocks; ++p)
    if(BLOCK[p] && BLOCK[p]->NEXT == Bold)
      BLOCK[p]->NEXT = Bnew;
  // 3  check for TYPES[] and FIRST
  for(bodytype t; t; ++t)
    if(TYPES[t] == Bold) TYPES[t] = Bnew;
  if(FIRST == Bold) FIRST = Bnew;
}
#endif
////////////////////////////////////////////////////////////////////////////////
// reset blocks' FIRST entries
void bodies::reset_firsts(int first[BT_NUM])
{
  for(bodytype t; t; ++t) {
    int F = first[t];
    for(block*B=TYPES[t]; B; B=B->next_of_same_type()) {
      B->set_first(F);
      F += B->N_bodies();
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
// # bodies not flagged to be ignored
unsigned bodies::N_subset() const
{
  if(!have(fieldbit::f)) return N_bodies();
  unsigned n = 0;
  LoopAllBodies(this,b) if(in_subset(b)) ++n;
  return n;
}
// delete all blocks and reset related data
void bodies::del_data() falcON_THROWING
{
  for(unsigned i=0; i!=index::max_blocks; ++i) if(BLOCK[i]) {
    falcON_DEL_O(BLOCK[i]);
    BLOCK[i] = 0;
  }
  NBLK = 0u;
  for(bodytype t; t; ++t) {
    NALL [t] = 0u;
    NBOD [t] = 0u;
    TYPES[t] = 0;
  }
  FIRST = 0;
}
////////////////////////////////////////////////////////////////////////////////
// destruction: delete all data
bodies::~bodies() falcON_THROWING
{
  if(debug(6)) DebugInfo("bodies::~bodies(): destructing bodies");
  BITS = fieldset::empty;
  if(C_FORTRAN)
    for(fieldbit f; f; ++f)
      const_cast<block*>(FIRST)->set_data_void(f,0);
  del_data();
}
////////////////////////////////////////////////////////////////////////////////
// set blocks' FIRST entries
void bodies::set_firsts()
{
  unsigned n = 0;
  for(const block* p=FIRST; p; p=p->next()) {
    const_cast<block*>(p)->set_first(n);
    n += p->N_bodies();
  }
}

////////////////////////////////////////////////////////////////////////////////
// set up blocks to hold N[t] bodies of type t
void bodies::set_data(const unsigned N[BT_NUM]) falcON_THROWING
{
  if(debug(5)) DebugInfo("bodies::set_data(): N=[%d,%d,%d], BITS=%s\n",
			 N[0],N[1],N[2],word(BITS));
  NBLK = 0u;
  NTOT = 0u;
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  try {
    block   *last = 0;
    unsigned i    = 0;
    for(bodytype t; t; ++t) {
      NBOD[t]  = NALL[t] = N[t];
      NTOT    += NBOD[t];
      NDEL[t]  = 0u;
      NNEW[t]  = 0u;
      TYPES[t] = 0;
      for(unsigned a,n=0u; n < NALL[t]; n+=a) {
	if(NBLK == index::max_blocks)
	  falcON_THROW("bodies: # blocks exceeds limit");
	a = min(NALL[t]-n, unsigned(index::max_bodies));
	block *b = new block(NBLK,a,a,i,t,BITS,this);
	if(debug(10)) DebugInfo("allocated %s @ %p\n",
				nameof(block),static_cast<void*>(b));
	i += a;
	if(last) last->link(b);
	last = b;
	if(n==0u) TYPES[t] = b;
	BLOCK[NBLK++] = b;
      }
    }
  } catch(falcON::exception E) {
    del_data();
    falcON_RETHROW(E);
  }
  FIRST = BLOCK[0];
}
////////////////////////////////////////////////////////////////////////////////
// construction 0: construction with N=0, but data fields
bodies::bodies(fieldset bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 )
{
  unsigned n[BT_NUM]={0u};
  if(debug(3))
    DebugInfo("bodies::bodies(): constructing bodies: n=%u,%u,%u, bits=%s",
	      n[0],n[1],n[2],word(BITS));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
// /////////////////////////////////////////////////////////////////////////////
// construction 1, new version
bodies::bodies(const unsigned n[BT_NUM],
	       fieldset       bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 )
{
  if(debug(3))
    DebugInfo("bodies::bodies(): constructing bodies: n=%u,%u,%u, bits=%s",
	      n[0],n[1],n[2],word(bits));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
// /////////////////////////////////////////////////////////////////////////////
// resets N, data; same as destruction followed by constructor 1            
void bodies::reset(const unsigned n[BT_NUM],
		   fieldset       bits) falcON_THROWING
{
  bool keepN = true;
  for(bodytype t; t; ++t) keepN = keepN && NALL[t] == n[t];
  if(keepN) {
    NTOT = 0u;
    for(bodytype t; t; ++t) {
      NBOD[t] = NALL[t];
      NDEL[t] = 0u;
      NNEW[t] = 0u;
      NTOT   += NBOD[t];
    }
    for(unsigned i=0; i!=index::max_blocks; ++i) 
      if(BLOCK[i]) BLOCK[i]->NBOD = BLOCK[i]->NALL;
    del_fields(BITS - bits);
    add_fields(bits - BITS);
  } else {
    del_data();
    BITS = bits;
    set_data(n);
  }
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
// construction 2:                                                          
// just make a copy of existing bodies:                                     
// - only copy data specified by 2nd argument                               
// - only copy bodies whose flags matches 3rd argument                      
bodies::bodies(bodies const&Other,
	       fieldset     copydata,
	       flags        copyflag) falcON_THROWING :
  BITS      ( copydata & Other.BITS ),
  C_FORTRAN ( 0 )
{
  if(copyflag && !Other.have_flag() ) 
    falcON_THROW("in bodies::bodies(): "
		 "copyflag !=0, but other bodies not supporting flag");
  unsigned n[BT_NUM];
  for(bodytype t; t; ++t) {
    if(copyflag) {
      LoopTypedBodies(&Other,i,t)
	if( flag(i).are_set(copyflag) ) ++(n[t]);
    } else 
      n[t] = Other.NBOD[t];
  }
  set_data(n);
  for(bodytype t; t; ++t) if(TYPES[t]) {
    block      *p =TYPES[t];
    const block*op=Other.TYPES[t];
    unsigned    oi=0;
    while(p && op && oi < op->N_bodies()) {
      p->copy(op,oi,copydata,copyflag);
      p = p->next();
    }
  }
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
// construction for C & FORTRAN support                                     
bodies::bodies(char, const unsigned n[BT_NUM]) falcON_THROWING
: BITS      ( fieldset::empty ),
  C_FORTRAN ( 1 )
{
  if(debug(3))
    DebugInfo("bodies::bodies(): constructing bodies for C & FORTRAN: n=%u,%u",
	      n[0],n[1]);
  for(bodytype t; t; ++t)
    if(n[t] > index::max_bodies)
      falcON_THROW("too many bodies\n");
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
void bodies::reset(char, fieldbit f, void*D) falcON_THROWING
{
  if(!C_FORTRAN || !FIRST || BLOCK[0] != FIRST)
    falcON_THROW("bodies::reset() called from wrongly initialized bodies");
  if(D) {
    char* DATA = static_cast<char*>(D);
    BITS |= fieldset(f);
    for(bodytype t; t; ++t)
      if( TYPES[t] && t.allows(f) ) {
	TYPES[t]->set_data_void(f,DATA);
	DATA += falcON::size(f) * TYPES[t]->NALL;
      }
  }
}
////////////////////////////////////////////////////////////////////////////////
void bodies::swap_bytes(fieldbit f) falcON_THROWING
{
  if(!BITS.contain(f))
    for(const block*p=FIRST; p; p=p->next())
      const_cast<block*>(p)->swap_bytes(f);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::add_field(fieldbit f) falcON_THROWING
{
  if(!BITS.contain(f)) {
    for(const block*p=FIRST; p; p=p->next())
      const_cast<block*>(p)->add_field(f);
    BITS |= fieldset(f);
    if(f == fieldbit::k) reset_keys();
  }
}
////////////////////////////////////////////////////////////////////////////////
void bodies::add_fields(fieldset b) falcON_THROWING
{
  if(!BITS.contain(b)) {
    for(const block *p=FIRST; p; p=p->next())
      const_cast<block*>(p)->add_fields(b);
    if(!BITS.contain(fieldbit::k) && b.contain(fieldbit::k)) reset_keys();
    BITS |= b;
  }
}
////////////////////////////////////////////////////////////////////////////////
void bodies::del_field(fieldbit f) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->del_field(f);
  BITS &= ~(fieldset(f));
}
////////////////////////////////////////////////////////////////////////////////
void bodies::del_fields(fieldset b) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->del_fields(b);
  BITS &= ~b;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::remove() falcON_THROWING {
  for(bodytype t; t; ++t)
    NBOD[t] = 0u;
  NTOT = 0u;
  for(block *p=FIRST; p; p=p->next()) {
    p->remove(NDEL[p->type()]);
    p->set_first(NTOT);
    NBOD[p->type()] += p->N_bodies();
    NTOT            += p->N_bodies();
  }
}
////////////////////////////////////////////////////////////////////////////////
// update FIRST and link TYPES[] together                                     
void bodies::link_blocks() {
  block **L = &FIRST, *P;
  for(bodytype t; t; ++t) {
    P = TYPES[t];
    if(P) {
      *L = P;
      while(P->NEXT && P->NEXT->TYPE == t) P = P->NEXT;
      L = &(P->NEXT);
    }
  }
  *L = 0;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::merge(bodies&Other) falcON_THROWING {
  if(NBLK + Other.NBLK > index::max_blocks)
    falcON_THROW("bodies::merge(): too many blocks\n");
  // loop Other.BLOCK[] and add them at the head of our TYPES[]                 
  for(unsigned n=0; n!=Other.NBLK; ++n) {
    block*B = Other.BLOCK[n];
    B->set_fields(BITS);
    BLOCK[NBLK] = B;
    B->NO   = NBLK++;
    B->NEXT = TYPES[B->TYPE];
    TYPES[B->TYPE]  = B;
    NALL [B->TYPE] += B->NALL;
    NBOD [B->TYPE] += B->NBOD;
    NTOT           += B->NBOD;
  }
  // link blocks together and reset FIRST                                       
  link_blocks();
  // set block::FIRST                                                           
  set_firsts();
  // finally reset all entries of OTHER                                         
  Other.FIRST = 0;
  for(bodytype t; t; ++t) {
    Other.TYPES[t] = 0;
    Other.NALL [t] = 0;
    Other.NBOD [t] = 0;
  }
  Other.NTOT = 0;
  Other.NBLK = 0;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::create(unsigned N, bodytype t) falcON_THROWING
{
  if(N > index::max_bodies) 
    falcON_THROW("bodies::create(): asked for %u > %u bodies\n",
		 N,  index::max_bodies);
  if(NBLK >= index::max_blocks)
    falcON_THROW("bodies::create(): number of blocks exceeded\n");
  // allocate new block and add it to the head of the list TYPES[t]             
  block* p = new block(NBLK,N,0,0,t,BITS,this);
  BLOCK[NBLK++]  = p;
  p->link(TYPES[t]);
  TYPES[t]  = p;
  NALL [t] += N;
  // link blocks together and reset FIRST                                       
  link_blocks();
  // set block::FIRST                                                           
  set_firsts();
  if(debug(2))
    DebugInfo("bodies::create(): created %u new bodies of type %s\n",
	      N,t.name());
}
////////////////////////////////////////////////////////////////////////////////
bodies::iterator bodies::new_body(bodytype t) falcON_THROWING
{
  if(0 == N_free(t)) {
    falcON_Warning("bodies::new_body(): no body available\n");
    return iterator(0);
  }
  for(const block* b=TYPES[t]; b; b=b->next_of_same_type())
    if(b->NALL > b->NBOD) {
      iterator i(b,const_cast<block*>(b)->NBOD++);
      set_firsts();
      NBOD[t]++;
      NNEW[t]++;
      NTOT++;
      if(have(fieldbit::f)) i.flag().add(flags::newbody);
      return i;
    }
  falcON_THROW("bodies::new_body(): cannot find free block\n");
  return iterator(0);
}
////////////////////////////////////////////////////////////////////////////////
falcON::real bodies::TotalMass(bodytype t) const
{
  if(!t || TYPES[t]==0 || !(TYPES[t]->has_field(fieldbit::m)) )
    return zero;
  real M(zero);
  for(const block* b=TYPES[t]; b; b=b->next_of_same_type())
    for(unsigned i=0; i!=b->N_bodies(); ++i)
      M += b->const_datum<fieldbit::m>(i);
  return M;
}
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// data I/O                                                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////
// read a single nemo snapshot:                                             
// - start at body position given                                           
// - we must have enough space to accomodate input                          
// - add_fields if required                                                 
// - type must match (best only read bodies of one type)                    
fieldset bodies::read_snapshot(snap_in  const&snap,
			       fieldset       want,
			       iterator const&start,
			       unsigned       Nread,
			       bool           warn) falcON_THROWING
{
  fieldset get = want & fieldset::nemoin;
  if(Nread == 0 || Nread > snap.Ntot()) Nread = snap.Ntot();
  if(start.my_index() + Nread > N_bodies())
    falcON_THROW("bodies::read_snapshot(): not enough space for data");
  fieldset read;
  // if phases given & x or v wanted, read them from phases
  if(get&fieldset::phases && snap.has(nemo_io::posvel)) {
    data_in inpt(snap,nemo_io::posvel);
    body b(start);
    b.read_posvel(inpt,get,Nread);
    if(inpt.N_read() != Nread)
      falcON_THROW("bodies::read_snapshot(): couldn't read all phase data");
    if(debug(2)) DebugInfo("bodies::read_snapshot(): phases read");
    read |= get & fieldset::phases;
    BITS |= get & fieldset::phases;
  }
  // now read data field by field
  for(fieldbit f; f; ++f) if(get.contain(f)) {
    if(debug(6))
      DebugInfo("bodies::read_snapshot(): f=%c: %s\n",letter(f),
		read.contain(f)? "already read" :
		!snap.has(nemo_io::field(f))? "not present" : "to be read");
    if( !read.contain(f) && snap.has(nemo_io::field(f)) ) {
      data_in inpt(snap,nemo_io::field(f));
      body b(start);
      b.read_data(inpt,Nread);
      if(min(Nread,inpt.N()) > inpt.N_read())
	falcON_THROW("bodies::read_snapshot(): "
		     "could only read %u of %u %c data",
		     inpt.N_read(), inpt.N(), letter(f));
      if(debug(2)) DebugInfo("bodies::read_snapshot(): %u %c read",
			     inpt.N_read(), letter(f));
      BITS |= fieldset(f);
      read |= fieldset(f);
    }
  }
  if(debug(1)) DebugInfo("bodies::read_snapshot(): read=%s\n",word(read));
  if(read & fieldset::source) mark_srce_data_changed();
  if(read & fieldset::sphmax) mark_sph_data_changed();
  if(warn && want != read)
    falcON_Warning("bodies::read_snapshot: couldn't read %s",
		   word(read.missing(want)));
  return read;
}
////////////////////////////////////////////////////////////////////////////////
void bodies::write_snapshot(snap_out const&snap,
			    fieldset       put,
			    iterator const&start,
			    unsigned       Nwrite) const falcON_THROWING
{
  if(this != start.my_bodies())
    falcON_THROW("bodies::write_snapshot(): start body is not ours");
  if(Nwrite == 0 || Nwrite > snap.Ntot()) Nwrite = snap.Ntot();
  if(start.my_index() + Nwrite > N_bodies())
    falcON_THROW("bodies::write_snapshot(): not enough data to write");
  put &= BITS;
  put &= fieldset::nemo;
  fieldset written;
  if(put&fieldset::p && put&fieldset::q) {
    data_out outp(snap,nemo_io::pot);
    body b(start);
    b.write_potpex(outp,Nwrite);
    if(outp.N_written() != Nwrite)
      falcON_THROW("bodies::write_snapshot(): couldn't write all pq data");
    if(debug(2)) DebugInfo("bodies::write_snapshot(): written pq");
    written |= fieldset::potent;
  }
  for(fieldbit f; f; ++f)
    if(put.contain(f) && !written.contain(f) &&
       (!falcON::is_sph(f) || N_sph() ) ) {
      data_out outp(snap,nemo_io::field(f));
      body b(start);
      b.write_data(outp,Nwrite);
      if(outp.N() != outp.N_written())
	falcON_THROW("bodies::write_snapshot(): "
		     "could only write %u of %u %c data",
		     outp.N_written(), outp.N(), letter(f));
      if(debug(2)) DebugInfo("bodies::write_snapshot(): written %u %c",
		outp.N_written(), letter(f));
      written |= fieldset(f);
    }
  if(debug(1)) DebugInfo("bodies::write_snapshot(): "
			 "written=%s for %u SPH & %u STD bodies\n",
			 word(written), N_sph(), N_std());
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
namespace {
  template<int BIT> 
  void Read(std::istream&in, body &B) {
    in >> B. template datum<BIT>();
  }
  template<int BIT> 
  void Write(std::ostream&to, const body &B) {
    to << ' ' << field_traits<BIT>::word()
       << '=' << B. template const_dat<BIT>();
  }
  typedef void (*p_reader)(std::istream&, body&);
  typedef void (*p_writer)(std::ostream&, const body&);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::read_simple_ascii(std::istream  &in,
			       const fieldbit*item,
			       unsigned       Ni,
			       const unsigned N[BT_NUM])
{
  // 1. create table of readers                                                 
  fieldset get;
  p_reader read [BT_NUM][100] = {0};
  p_writer write[BT_NUM][100] = {0};
  if(Ni > 100) {
    Ni = 100;
    falcON_Warning(" can only read the first 100 data entries\n");
  }
  for(int i=0; i!=Ni; ++i) {
    if(debug(6)) DebugInfo("bodies::read_simple_ascii(): item[%d]=%c\n",
			   i,letter(item[i]));
    if(get.contain(item[i]))
      falcON_Warning("bodies::read_simple_ascii(): "
		     "reading item '%c' more than once",
	      letter(item[i]));
    get |= fieldset(item[i]);
    switch(value(item[i])) {
#define SET_READ(BIT,NAME)				\
    case BIT:						\
      for(bodytype t; t; ++t)				\
        if(t.allows(BIT)) {				\
          read[t][i] = &Read<BIT>;			\
          if(debug(20)) write[t][i]= &Write<BIT>;	\
	}						\
      break;
      DEF_NAMED(SET_READ);
#undef SET_READ
    default: 
      for(bodytype t; t; ++t)
        read[t][i] = 0;
    }
  }
  if(Ni==100)
    falcON_Warning("bodies::read_simple_ascii(): cannot read >100 items\n");
  // 2. reset N & add fields                                                    
  reset(N, BITS|get);
  // 3. loop bodies & read data whereby ignoring lines starting with '#'        
  for(bodytype t; t; ++t)
    if(N[t]) {
      if(debug(4))
	DebugInfo("bodies::read_simple_ascii(): now reading %d %s bodies...\n",
		  N[t],t.name());
      LoopTypedBodies(this,Bi,t) {
	while( in && eat_line(in,'#') );
	if(!in) falcON_THROW("bodies::read_simple_ascii(): "
			     "end of input before data have been read");
	for(int i=0; i!=Ni; ++i)
	  if(read[t][i]) {
	    read [t][i](in,Bi);
	    if(write[t][i]) write[t][i](std::cerr,Bi);
	  }
	if(debug(20)) std::cerr<<'\n';
	SwallowRestofLine(in);
      }
    }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// sorted index table                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void bodies::sorted(Array<index>&table, 
		    real       (*func)(iterator const&)) const falcON_THROWING
{
  const int n = N_subset();
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  if(have(fieldbit::f)) {
    int i = 0;
    LoopSubsetBodies(this,b) {
      I[i] = static_cast<index>(b);
      Q[i] = func(b);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = static_cast<index>(b);
      Q[i] = func(b);
      ++i;
    }
  }
  int*R = falcON_NEW(int,n);
  HeapIndex(Q,n,R);
  table.reset(n);
  for(int i=0; i!=n; ++i)
    table[i] = I[R[i]];
  falcON_DEL_A(Q);
  falcON_DEL_A(I);
  falcON_DEL_A(R);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::sorted(Array<index>&table, 
		    Array<real> &quant, 
		    real       (*func)(iterator const&)) const falcON_THROWING
{
  const int n = N_subset();
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  if(have(fieldbit::f)) {
    int i = 0;
    LoopSubsetBodies(this,b) {
      I[i] = static_cast<index>(b);
      Q[i] = func(b);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = static_cast<index>(b);
      Q[i] = func(b);
      ++i;
    }
  }
  int*R = falcON_NEW(int,n);
  HeapIndex(Q,n,R);
  table.reset(n);
  quant.reset(n);
  for(int i=0; i!=n; ++i) {
    table[i] = I[R[i]];
    quant[i] = Q[R[i]];
  }
  falcON_DEL_A(Q);
  falcON_DEL_A(I);
  falcON_DEL_A(R);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::snapshot                                                       
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace {
  class PointerBank {
    //--------------------------------------------------------------------------
    struct PterWithKey {
      friend class falcON::traits<PterWithKey>;
      const void  *pter;
      char        *key,*name;
      size_t       size;
      PterWithKey *next;
      //........................................................................
      PterWithKey(const void* p, const char*k, size_t s, const char*n,
		  PterWithKey*x)
	: pter(p),
	  size(s),
	  next(x),
	  key (falcON_NEW(char, strlen(k)+strlen(n)+2)),
	  name(key + strlen(k) + 1) {
	strcpy(key ,k);
	strcpy(name,n);
      }
      //........................................................................
      ~PterWithKey() {
	falcON_DEL_A(key);
      }
    } *HEAD;
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    /// default constructor
    PointerBank() : HEAD(0) {}
    //--------------------------------------------------------------------------
    /// copy constructor
    PointerBank(PointerBank const&PB) : HEAD(0) {
      for(PterWithKey*P=PB.HEAD; P; P=P->next)
	HEAD = new PterWithKey(P->pter, P->key, P->size, P->name, HEAD);
    }
    //--------------------------------------------------------------------------
    /// destructor
    ~PointerBank() {
      PterWithKey*P=HEAD,*N;
      while(P) {
	N=P->next;
	falcON_DEL_O(P);
	P=N;
      }
    }
    //--------------------------------------------------------------------------
    /// add a pointer: key must not yet be known in bank
    void add(const void*p, const char* k, size_t s, const char* n) {
      for(PterWithKey*P=HEAD; P; P=P->next)
	if(0==strcmp(P->key, k))
	  falcON_THROW("snapshot::add_pointer(): "
		       "key '%s' is already in bank\n",k);
      HEAD = new PterWithKey(p,k,s,n,HEAD);
    }
    //--------------------------------------------------------------------------
    /// set a pointer: add if new key, else replace (type & size must match)
    void set(const void*p, const char* k, size_t s, const char* n) {
      for(PterWithKey*P=HEAD; P; P=P->next)
	if(0==strcmp(P->key, k)) {
	  if(strcmp(P->name, n))
	    falcON_THROW("snapshot::set_pointer(): "
			 "name mismatch ('%s' : '%s')",P->name,n);
	  if(P->size != s)
	    falcON_THROW("snapshot::set_pointer(): "
			 "size mismatch (%u : %u)",P->size,s);
	  P->pter = p;
	  return;
	}
      HEAD = new PterWithKey(p,k,s,n,HEAD);
    }
    //--------------------------------------------------------------------------
    /// delete an entry from the bank
    void del(const char* k, bool warn = 0) {
      PterWithKey **PP=&HEAD, *P=HEAD;
      for(; P; PP=&(P->next), P=P->next)
	if(0==strcmp(P->key, k)) {
	  (*PP) = P->next;
	  falcON_DEL_O(P);
	}
      if(warn)
	falcON_Warning("snapshot::del_pointer(): "
		       "key '%s' not found in bank\n",k);
    }
    //--------------------------------------------------------------------------
    /// return a pointer referred to by a given key
    const void*get(const char*k, size_t s, const char*n, const char*func) const
      falcON_THROWING {
      for(PterWithKey*P=HEAD; P; P=P->next)
	if(0==strcmp(P->key, k)) {
	  if(s != P->size)
	    falcON_THROW("snapshot::%s(): "
			 "size (%u) does not match value in bank (%u)\n",
			 func,s,P->size);
	  if(strcmp(n,P->name))
	    falcON_THROW("snapshot::%s(): "
			 "name (%s) does not match value in bank (%s)\n",
			 func,n,P->name);
	  return P->pter;
	}
      return 0;
    }
  };// class PointerBank
} // namespace {
falcON_TRAITS(::PointerBank,"{body.cc}::PointerBank");
#if(0) // gcc 3.3.5 doesn't like this
falcON_TRAITS(::PointerBank::PterWithKey,
	      "{body.cc}::PointerBank::PterWithKey");
#endif
////////////////////////////////////////////////////////////////////////////////
void snapshot::__add_pointer(const void*p,
			     const char*k,
			     size_t     s,
			     const char*n) const falcON_THROWING
{
  if(debug(4))
    DebugInfo("snapshot::add_pointer() %p to '%s' under \"%s\"\n",p,n,k);
  if(p) {
    if(PBNK == 0) const_cast<snapshot*>(this)->PBNK = new PointerBank();
    static_cast<PointerBank*>(PBNK)->add(p,k,s,n);
  } else if(PBNK)
    // NULL pointer: just check match of name & type and non-existence
    if(static_cast<PointerBank*>(PBNK)->get(k,s,n,"add_pointer"))
      falcON_THROW("snapshot::add_pointer(): key '%s' is already in bank\n",k);
}
////////////////////////////////////////////////////////////////////////////////
void snapshot::__set_pointer(const void*p,
			     const char*k,
			     size_t     s,
			     const char*n) const falcON_THROWING
{
  if(debug(4))
    DebugInfo("snapshot::set_pointer() %p to '%s' under \"%s\"\n",p,n,k);
  if(p) {
    // non-NULL pointer: add or replace (if type & size match)
    if(PBNK == 0) const_cast<snapshot*>(this)->PBNK = new PointerBank();
    static_cast<PointerBank*>(PBNK)->set(p,k,s,n);
  } else if(PBNK)
    // NULL pointer: delete from bank
    static_cast<PointerBank*>(PBNK)->del(k);
}
////////////////////////////////////////////////////////////////////////////////
const void* snapshot::__get_pointer(const char*k,
				    size_t     s,
				    const char*n) const falcON_THROWING
{
  const void*p = PBNK? static_cast<PointerBank*>(PBNK)->get(k,s,n,
							    "get_pointer") : 0;
  if(debug(4))
    DebugInfo("snapshot::get_pointer() %p to '%s' under \"%s\"\n",p,n,k);
  return p;
}
////////////////////////////////////////////////////////////////////////////////
void snapshot::del_pointer(const char*k) const
{
  if(debug(4))
    DebugInfo("snapshot::del_pointer() under \"%s\"\n",k);
  if(PBNK) static_cast<PointerBank*>(PBNK)->del(k);
}
////////////////////////////////////////////////////////////////////////////////
snapshot::snapshot(snapshot const&S,
		   fieldset       Bd,
		   flags          F) falcON_THROWING
: bodies ( S,Bd,F ),
  INIT   ( S.INIT ),
  TINI   ( S.TINI ),
  TIME   ( S.TIME ),
  PBNK   ( S.PBNK? new PointerBank(*(static_cast<PointerBank*>(S.PBNK))) : 0 ),
  PARA   ( 0 )
{}
////////////////////////////////////////////////////////////////////////////////
#if(0) // not yet implemented due to bodies::copy() missing
void snapshot::copy(snapshot const&S,
		    fieldset       Bd,
		    flags          F) falcON_THROWING
{
  bodies::copy(S,Bd,F);
  TIME = S.TIME;
  if(PBNK) falcON_DEL_O(static_cast<PointerBank*>(PBNK));
  PBNK = S.PBNK? new PointerBank(*(static_cast<PointerBank*>(S.PBNK))) : 0;
}
#endif
////////////////////////////////////////////////////////////////////////////////
snapshot::~snapshot()
{
  if(PBNK) { falcON_DEL_O(static_cast<PointerBank*>(PBNK)); PBNK = 0; }
}
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_NEMO
bool snapshot::read_nemo(                          // R: was time in range?     
			 nemo_in const&i,          // I: nemo input             
			 fieldset     &r,          // O: what has been read     
			 fieldset      g,          //[I: what to read]          
			 const char   *t,          //[I: time range]            
			 bool          w)          //[I: warn: missing data]    
  falcON_THROWING
{
  if(!i.has_snapshot())
    falcON_THROW("snapshot::read_nemo(): no snapshot to read");
  snap_in s(i);

  if(s.has_time()) {
    if(t && !time_in_range(s.time(),t)) {
      r = fieldset::empty;
      return false;
    }
    TIME = s.time();
    if(!INIT) { 
      TINI = TIME;
      INIT = true;
    }
  } else
    TIME = 0.;
  bool need_reset = false;
  for(bodytype t; t; ++t)
    if(s.Nbod(t) != N_bodies(t)) need_reset = true;
  if(need_reset) reset(s.Nbod(), fieldset::o);
  r = read_snapshot(s,g,begin_all_bodies(),N_bodies(),w);
  return true;
}
////////////////////////////////////////////////////////////////////////////////
fieldset snapshot::read_part(                      // R: what has been read     
			     snap_in  const&s,     // I: nemo input             
			     fieldset       g,     // I: what to read           
			     iterator const&b,     // I: start position         
			     bool           w,     //[I: warn: missing data]    
			     unsigned       n)     //[I: #, def: all in input]  
  falcON_THROWING
{
  if(s.has_time()) {
    TIME = s.time();
    if(!INIT) { 
      TINI = TIME;
      INIT = true;
    }
  } else
    TIME = 0.;
  return read_snapshot(s,g,b,n,w);
}
////////////////////////////////////////////////////////////////////////////////
void snapshot::write_nemo(nemo_out const&o,        // I: nemo output            
			  fieldset       w,        // I: what to write          
			  iterator const&b,        // I: starting here          
			  unsigned       n) const  //[I: #, default: all]       
  falcON_THROWING
{
  unsigned i = bodyindex(b);
  if(this != b.my_bodies())
    falcON_THROW("snapshot::write_nemo() start body is not ours\n");
  if(n == 0) n = N_bodies()-i;
  else if(i + n > N_bodies()) {
    falcON_Warning("snapshot::write_nemo() cannot write %u bodies, "
		   "will only write %u\n",n,N_bodies()-i);
    n = N_bodies()-i;
  }
  unsigned nb[BT_NUM]={0}, nt=n, nc(0u);
  for(bodytype t; t; ++t)
    if(i < (nc+=N_bodies(t))) {
      nb[t] = min(nc-i, nt);
      i  += nb[t];
      nt -= nb[t];
    }
  snap_out s(o,nb,TIME);
  write_snapshot(s,w,b,n);
}
////////////////////////////////////////////////////////////////////////////////
void snapshot::write_nemo(nemo_out const&o,        // I: nemo output            
			  fieldset       w) const  // I: what to write          
  falcON_THROWING
{
  snap_out s(o,N_bodies_per_type(),TIME);
  write_snapshot(s,w,begin_all_bodies(),N_bodies());
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_REAL_IS_FLOAT
//------------------------------------------------------------------------------
#define READ(BIT)							\
  if(!is_sph(BIT) && nd || ns) {					\
    unsigned nr = is_sph(BIT)? ns : ns+nd;				\
    unsigned fr = rec+rec + nr*field_traits<BIT>::size;			\
    if(fsze+fr > fsize) {						\
      fieldset miss = fgot.missing(read);				\
      falcON_Warning("bodies::read_gadget(): "				\
		     "reached end of file \"%s\"; "			\
		     "cannot read %s\n", file, word(miss));		\
      read &= fgot;							\
      goto NextFile;							\
    }									\
    FortranIRec F(in, rec, swap);					\
    if(read.contain(BIT)) {						\
      if(F.size() != nr*field_traits<BIT>::size)			\
	falcON_THROW("bodies::read_gadget(): mismatch reading %u %c: "	\
		     "expected %u bytes, found %u\n",			\
		     nr,field_traits<BIT>::word(),			\
		     nr*field_traits<BIT>::size,F.size());		\
      add_field(BIT);							\
      if(ns) {								\
	body sph(SPH);							\
	sph.read_Fortran(F, BIT, ns, swap);				\
      }									\
      if(nd && !is_sph(BIT)) {						\
	body std(STD);							\
        std.read_Fortran(F, BIT, nd, swap);				\
      }									\
      if(debug(2)) DebugInfo("bodies::read_gadget(): read %u %c\n",	\
	                     nr, field_traits<BIT>::word());		\
      fgot |= fieldset(BIT);						\
    } else {								\
      F.skip_bytes(F.size());						\
      if(debug(3)) DebugInfo("bodies::read_gadget(): skip %u %c\n",	\
			     nr, field_traits<BIT>::word());		\
    }									\
    fsze += fr;								\
  }
//------------------------------------------------------------------------------
double bodies::read_gadget(const char*fname,
			   fieldset   read,
			   fieldset  &got,
			   unsigned   rec) falcON_THROWING
{
  read &= fieldset("mxvkURHpa");
  got   = fieldset::empty;
  // 1 open first data file, read header, and determine number of data files
  GadgetHeader header0,headeri,*header=&header0;
  input        in;
  int          nfile=1;
  char         filename[256];
  const char  *file=0;
  bool         swap;
  size_t       fsize = FileSize(fname);
  if(fsize) {
    // 1.1 try single file "fname"
    file = fname;
    try {
      in.open(file);
      if(! header->Read(in, rec, swap))
      falcON_THROW("bodies::read_gadget(): "
		   "file \"%s\" not a gadget data file\n",file);
    }
    catch(exception E) { falcON_RETHROW(E); }
    if(swap && debug(1))
      DebugInfo("bodies::read_gadget(): need to swap bytes ...\n");
    nfile = header->num_files;
    if(nfile==0) nfile=1;
    if(debug(1)) {
      DebugInfo("bodies::read_gadget(): header read from file \"%s\":\n",
		 fname);
      header->dump(std::clog);
    }
    if(nfile!=1)
      falcON_THROW("bodies::read_gadget(): num_files=%u, expected 0 or 1\n",
		   nfile);
    header->check_simple_npart_error();
  } else {
    // 1.2 try file "fname.0" (with possibly more to follow)
    SNprintf(filename,256,"%s.%u",fname,0);
    fsize = FileSize(filename);
    in.open(filename);
    if(!in) falcON_THROW("bodies::read_gadget(): cannot open file \"%s\" "
			 "nor file \"%s\"\n", fname, filename);
    file = filename;
    try {
      if(! header->Read(in, rec, swap) )
	falcON_THROW("bodies::read_gadget(): "
		     "file \"%s\" not a gadget data file\n",file);
    }
    catch(exception E) { falcON_RETHROW(E); }
    if(swap && debug(1))
      DebugInfo("bodies::read_gadget(): need to swap bytes ...\n");
    nfile = header->num_files;
    if(debug(1)) {
      DebugInfo("bodies::read_gadget(): header read from file \"%s\":\n",file);
      header->dump(std::clog);
    }
  }
  // 2 establish number of SPH and non-SPH particles and allocate memory.
  unsigned NB[BT_NUM] = {0}, NP[6] ={0};
  NB[bodytype::gas] = header->npartTotal[0];
  for(int k=1; k!=6; ++k) NB[bodytype::std] += header->npartTotal[k];
  if(NB[bodytype::gas] == 0) read &= fieldset(fieldset::nonSPH);
  reset(NB,read);
  // 3 loop data files
  body SPH = begin_typed_bodies(bodytype::gas);
  body STD = begin_typed_bodies(bodytype::std);
  for(int ifile=0; ifile!=nfile; ) {
    fieldset fgot;
    size_t   fsze = rec+rec + sizeof(GadgetHeader);
    for(int k=0; k!=6; ++k) {
      NP[k] += header->npart[k];
      if(NP[k] > header->npartTotal[k])
	falcON_THROW("bodies::read_gadget(): corrupted data file(s): "
		     "Sum npart[%u]=%u > npartTotal[%u]=%u\n",
		     k,NP[k],k,header->npartTotal[k]);
    }
    // 3.1 determine number of sph and std bodies in this file
    unsigned ns = header->npart[0], nd=0, nm=0;
    for(int k=0; k!=6; ++k) {
      if(k) nd += header->npart[k];
      if(header->masstab[k] == 0) nm += header->npart[k];
    }
    // 3.2 read positions
    READ(fieldbit::x);
    if(read == fgot) goto NextFile;
    // 3.3 read velocties
    READ(fieldbit::v);
    if(read == fgot) goto NextFile;
    // 3.4 read keys
    READ(fieldbit::k);
    if(read == fgot) goto NextFile;
    // 3.5 read masses --- OR assign them ...
    if(nm) {
      unsigned fr = rec+rec + nm*sizeof(real);
      if(fsze+fr > fsize) {
	fieldset miss = fgot.missing(read);
	falcON_Warning("bodies::read_gadget(): reached end of file \"%s\"; "
		       "cannot read %s\n", file, word(miss));
	read &= fgot;
	goto NextFile;
      }
      FortranIRec F(in, rec, swap);
      if(read.contain(fieldbit::m)) {
	if(F.size() != nm*sizeof(real))
	  falcON_THROW("bodies::read_gadget(): mismatch reading %u m: "
		       "expected %u bytes, found %u\n",
		       nm,nm*sizeof(real), F.size());
	body sph(SPH), std(STD);
	add_field(fieldbit::m);
	if(header->npart[0]) {
	  if(header->masstab[0])
	    for(int b=0; b!=header->npart[0]; ++b,++sph)
	      sph.mass() = header->masstab[0];
	  else
	    sph.read_Fortran(F, fieldbit::m, header->npart[0], swap);
	}
	for(int k=1; k!=6; ++k) if(header->npart[k]) {
	  if(header->masstab[k])
	    for(int b=0; b!=header->npart[k]; ++b,++std)
	      std.mass() = header->masstab[k];
	  else
	    std.read_Fortran(F, fieldbit::m, header->npart[k], swap);
	}
	if(debug(2)) DebugInfo("bodies::read_gadget(): read %u m\n", nm);
	fgot |= fieldset(fieldbit::m);
      } else {
	F.skip_bytes(F.size());
	if(debug(3)) DebugInfo("bodies::read_gadget(): skip %u m\n", nm);
      }
      fsze += fr;
    } else if(read.contain(fieldbit::m)) {
      body sph(SPH), std(STD);
      add_field(fieldbit::m);
      if(header->npart[0])
	for(int b=0; b!=header->npart[0]; ++b,++sph)
	  sph.mass() = header->masstab[0];
      for(int k=1; k!=6; ++k) if(header->npart[k]) {
	for(int b=0; b!=header->npart[k]; ++b,++std)
	  std.mass() = header->masstab[k];
      }
      fgot |= fieldset(fieldbit::m);
    }
    if(read == fgot) goto NextFile;
    // 3.6 read gas internal energies
    READ(fieldbit::U);
    if(read == fgot) goto NextFile;
    // 3.7 read gas densities
    READ(fieldbit::R);
    if(read == fgot) goto NextFile;
    // 3.8 read SPH smoothing lengths
    READ(fieldbit::H);
    if(read == fgot) goto NextFile;
    // 3.9 read potentials
    READ(fieldbit::p);
    if(read == fgot) goto NextFile;
    // 3.10 read accelerations
    READ(fieldbit::a);
    if(read == fgot) goto NextFile;
  NextFile:
    got |= fgot;
    // 3.11 open next data file (unless all have been read) and move SPH & STD
    if(++ifile != nfile) {
      SPH += ns;
      STD += nd;
      SNprintf(filename,256,"%s.%u",fname,ifile);
      in.open(file);
      if(!in) falcON_THROW("bodies::read_gadget(): cannot open file \"%s\"\n",
			   file);
      bool swap_old = swap;
      try { 
	if(! headeri.Read(in, rec, swap))
	  falcON_THROW("bodies::read_gadget(): "
		       "file \"%s\" not a gadget data file\n",file);
      }
      catch(exception E) { falcON_RETHROW(E); }
      if(swap_old != swap)
	falcON_THROW("bodies::read_gadget(): different endianess "
		     "amongst data files for the same snapshot\n");
      if(header0.mismatch(headeri))
	falcON_THROW("bodies::read_gadget(): header mismatch\n");
      header = &headeri;
      if(debug(2)) {
	DebugInfo("bodies::read_gadget(): header read from file \"%s\":\n",
		   file);
	header->dump(std::clog);
      }
    }
  }
  if(debug(2))
    DebugInfo("bodies::read_gadget(): read %s for %u SPH & %u STD bodies\n",
	      word(got), NB[bodytype::gas], NB[bodytype::std]);
  return header->time;
}
////////////////////////////////////////////////////////////////////////////////
#define WRITE(BIT)							  \
  if(!is_sph(BIT) || N_sph()) {						  \
    unsigned nw = is_sph(BIT)? N_sph() : N_bodies();			  \
    FortranORec F(out, nw * field_traits<BIT>::size, rec);		  \
    if(have(BIT)) {							  \
      if(N_sph())							  \
	begin_typed_bodies(bodytype::gas).write_Fortran(F, BIT, N_sph()); \
      if(!is_sph(BIT) && N_std())					  \
        begin_typed_bodies(bodytype::std).write_Fortran(F, BIT, N_std()); \
      if(debug(2)) DebugInfo("bodies::write_gadget(): written %u %c\n",	  \
			     nw, field_traits<BIT>::word());		  \
    } else {								  \
      if(warn)								  \
	falcON_Warning("bodies::write_gadget(): "			  \
                       "don't have %c, write out zeros\n",		  \
                       field_traits<BIT>::word());			  \
      F.fill_bytes(nw);							  \
      if(debug(2))							  \
        DebugInfo(2,"bodies::write_gadget(): written %u 0 for %c\n",	  \
		  nw, field_traits<BIT>::word());			  \
    }									  \
    written |= fieldset(BIT);						  \
  }
//------------------------------------------------------------------------------
void bodies::write_gadget(output&out, double time, fieldset write,
			  bool warn, unsigned rec) const falcON_THROWING
{
  // ensure we have keys ("ids" in gadget)
  write |= fieldset("mxvkU");
  const bool had_keys = have(fieldbit::k);
  if(!had_keys) const_cast<bodies*>(this)->add_field(fieldbit::k);
  // initialize and write header
  GadgetHeader header0, *header=&header0;
  header->num_files = 1;
  header->npart[0] = header->npartTotal[0] = N_sph();
  header->npart[1] = header->npartTotal[1] = N_std();
  FortranORec::Write(out,header,1,rec);
  // write out data
  fieldset written;
  WRITE(fieldbit::x);               // positions:   compulsary in gadget
  WRITE(fieldbit::v);               // velocities:  compulsary in gadget
  WRITE(fieldbit::k);               // keys:        compulsary in gadget
  WRITE(fieldbit::m);               // masses:      compulsary in gadget
  WRITE(fieldbit::U);               // U_internal:  compulsary in gadget
  if(write & fieldset("RHpa")) {
  WRITE(fieldbit::R);               // gas density:       optional
  if(write & fieldset("Hpa")) {
  WRITE(fieldbit::H);               // smoothing lengths: optional
  if(write & fieldset("pa")) {
  WRITE(fieldbit::p);               // potentials:        optional
  if(write & fieldset::a)
  WRITE(fieldbit::a);               // accelerations:     optional
  } } }
  if(debug(1))
    DebugInfo("bodies::write_gadget(): written %s for "
	      "%u SPH & %u STD bodies\n", word(written), N_sph(), N_std());
  // make sure we have no keys if we didn't have them before ...
  if(!had_keys) const_cast<bodies*>(this)->del_field(fieldbit::k);
}
//------------------------------------------------------------------------------
#endif // falcON_REAL_IS_FLOAT
////////////////////////////////////////////////////////////////////////////////
