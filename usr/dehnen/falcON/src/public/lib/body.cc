// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.cc                                                                     |
//                                                                             |
// Copyright (C) 2000-2007 Walter Dehnen                                       |
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
void bodies::block::reset_flags() const
{
  if(0 != DATA[fieldbit::f]) {
    if(TYPE.is_sph()) 
      for(int n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sph;
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
    falcON_ExceptF("flags not supported","bodies::flag_all_as_active()");
}
////////////////////////////////////////////////////////////////////////////////
void bodies::block::flag_as_sph() const falcON_THROWING
{
  if(0 != DATA[fieldbit::f]) {
    for(int n=0; n!=NALL; ++n)
      datum<fieldbit::f>(n).add(flags::sph);
  } else
    falcON_ExceptF("flags not supported","bodies::flag_as_sph()");
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
inline
void bodies::block::add_field (fieldbit f) falcON_THROWING {
  if(TYPE.allows(f) && 0 == DATA[value(f)] ) {
    debug_info(4,"bodies::block::add_field(): allocating data for %u %c (%s)\n",
	       NALL,letter(f),name(f));
    set_data_void(f, falcON_NEW(char,NALL*falcON::size(f)));
    if(f == fieldbit::f) reset_flags();
  }
}
////////////////////////////////////////////////////////////////////////////////
inline
void bodies::block::del_field (fieldbit f) falcON_THROWING {
  if(DATA[value(f)]) {
    debug_info(4,"bodies::block::del_field(): "
	       "de-allocating data for %c (%s)\n",letter(f),name(f));
    falcON_DEL_A(static_cast<char*>(DATA[value(f)]));
  }
  set_data_void(f,0);
}
////////////////////////////////////////////////////////////////////////////////
inline
void bodies::block::add_fields(fieldset b) falcON_THROWING {
  for(fieldbit f; f; ++f)
    if(b.contain(f)) add_field(f);
}
////////////////////////////////////////////////////////////////////////////////
inline
void bodies::block::del_fields(fieldset b) falcON_THROWING {
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) del_field(f);
}
////////////////////////////////////////////////////////////////////////////////
inline
void bodies::block::set_fields(fieldset b) falcON_THROWING {
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) add_field(f);
    else             del_field(f);
}
////////////////////////////////////////////////////////////////////////////////
inline
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
  : NALL ( na ), 
    NBOD ( nb ), 
    NO   ( no ),
    TYPE ( type ),
    FIRST( fst ),
    NEXT ( 0 ),
    BODS ( bods )
{
  if(na<nb)
    falcON_ExceptF("N_alloc < N_bodies","bodies::block::block()");
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
  static const unsigned BD = 1<<BIT;
  static void copy(void    **data,
		   unsigned  from,
		   unsigned  to  ,
		   fieldset  b,
		   fieldset &c) {
    if(data[BIT] && b.contain(fieldbit(BIT)) ) {
      memcpy(static_cast<      char*>(data[BIT])+to  *BodyData::ZQUANT[BIT],
	     static_cast<const char*>(data[BIT])+from*BodyData::ZQUANT[BIT],
	     BodyData::ZQUANT[BIT]);
      c |= fieldset(1<<BIT);
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
  fieldset copied(0);
  if(from != to)
    CopyBody<0>::copy(DATA,from,to,b,copied);
  return copied;
}
////////////////////////////////////////////////////////////////////////////////
template<unsigned BIT=0, unsigned END=BodyData::NQUANT> struct CopyBodies {
  static const unsigned BD = 1<<BIT;
  static void copy(void*const*data_fr,
		   void*const*data_to,
		   unsigned  fr,
		   unsigned  to,
		   unsigned  num,
		   fieldset  b,
		   fieldset&c) {
    if(data_fr[BIT] && data_to[BIT] && b & fieldset(1<<BIT) ) {
      memcpy(static_cast<      char*>(data_fr[BIT])+to*BodyData::ZQUANT[BIT],
	     static_cast<const char*>(data_to[BIT])+fr*BodyData::ZQUANT[BIT],
	     num*BodyData::ZQUANT[BIT]);
      c |= fieldset(1<<BIT);
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
  fieldset copied(0);
  if(this == other)
    falcON_ExceptF("this == other","bodies::block::copy_bodies()");
  else
    CopyBodies<0>::copy(other->DATA,DATA,from,to,num,copy,copied);
  return copied;
}
////////////////////////////////////////////////////////////////////////////////
inline void bodies::block::skip(unsigned&from,
				flags    copyflag) const falcON_THROWING
{
  if(copyflag) {
    if(! has_field(fieldbit::f) )
      falcON_ExceptF("copyflag!=0 but flags not supported",
		     "bodies::block::copy()");
    for(; from<NBOD && !(flag(from).are_set(copyflag)); ++from );
  }
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
    falcON_ExceptF("cannot copy from self","bodies::block::copy()");
  NBOD = 0u;
  if( From == 0)
    return fieldset(0);
  unsigned copy;
  unsigned free = NALL;
  fieldset copied (0);
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
	falcON_ExceptF("cannot copy from self","bodies::block::copy()");
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
    falcON_ExceptF("flags needed but not supported","bodies::remove()");
  unsigned lo=0u, hi=NBOD-1;
  while(lo < hi) {
    while(! to_remove(const_datum<fieldbit::f>(lo)) && lo < hi)
      ++lo;
    while(  to_remove(const_datum<fieldbit::f>(hi)) && lo < hi) {
      ++removed;
      --hi;
    }
    if(lo < hi) copy_body(hi--,lo++);
  }
  if(lo == hi && ! to_remove(const_datum<fieldbit::f>(lo))) ++lo;
  NBOD = lo;
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
    debug_info(1,"bodies::block::read_data(%c): must convert from %s to %s",
	       letter(f),
	       nemo_io::type_name(nemo_io::NotReal),
	       nemo_io::type_name(nemo_io::Real));
    unsigned ntot = N*inpt.sub_N();
    notreal* data = falcON_NEW(notreal, ntot);
    inpt.read(data, N);
    const notreal* d = data;
    real         * D = static_cast<real*>(DATA[value(f)]);
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
  if(coerce)
    debug_info(1,"bodies::read_posvel(): must convert from %s to %s",
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
  debug_info(2,"bodies::block::read_posvel(): read %s",
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
void bodies::block::read_Fortran(FortranIRec&I, fieldbit f, unsigned from,
				 unsigned N) falcON_THROWING
{
  if(!TYPE.allows(f))
    falcON_THROW("bodies::block::read_Fortran(%c): not allowed by our type",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::read_Fortran(%c): cannot read that many",
		 letter(f));
  add_field(f);
  unsigned R = I.read_bytes(static_cast<char*>(DATA[value(f)])
			    +from*falcON::size(f), N*falcON::size(f));
  if(R != N*falcON::size(f))
    falcON_THROW("bodies::block::read_Fortran(%c): "
		 "could only read %u of %u bytes\n",R,N*falcON::size(f));
  debug_info(4,"bodies::block::read_Fortran(): read %u %s\n",N,name(f));
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
  debug_info(4,"bodies::block::write_Fortran(): written %u %s\n",N,name(f));
}
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
bodies::iterator& bodies::iterator::read_Fortran(FortranIRec&I,
						 fieldbit f, unsigned R)
  falcON_THROWING
{
  if(R * falcON::size(f) > I.bytes_unread())
    falcON_THROW("body::read_Fortran: want %u %s (%u bytes) but "
		 "only %u bytes left on Fortran record\n",
		 R, name(f), R*falcON::size(f), I.bytes_unread());
  while(is_valid() && R) {
    unsigned r = min(N-K, R);
    const_cast<block*>(B)->read_Fortran(I,f,K,r);
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
    falcON_THROW("body::write_Fortran: want %u %s (%u bytes) but "
		 "only %u bytes left free on Fortran record\n",
		 W, name(f), W*falcON::size(f), O.bytes_free());
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
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::bodies                                                         
//                                                                              
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
  debug_info(6,"bodies::~bodies(): destructing bodies");
  BITS = fieldset(0);
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
void bodies::set_data(const unsigned *N) falcON_THROWING
{
  NBLK = 0u;
  NTOT = 0u;
  try {
    block   *last = 0;
    unsigned i    = 0;
    for(bodytype t; t; ++t) {
      NBOD[t] = NALL[t] = N? N[t] : 0u;
      NTOT   += NBOD[t];
      NDEL[t] = 0u;
      NNEW[t] = 0u;
      TYPES[t] = 0;
      for(unsigned a,n=0u; n < NALL[t]; n+=a) {
	if(NBLK == index::max_blocks)
	  falcON_ExceptF("# blocks exceeds limit","bodies");
	a = min(NALL[t]-n, unsigned(index::max_bodies));
	block *b = new block(NBLK,a,a,i,t,BITS,this);
	i+= a;
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
// construction 1:                                                          
// allocate Nbod (1st arg) bodies with properties given by 2nd argument     
// 3rd argument gives # SPH bodies within Nbod                              
bodies::bodies(unsigned nb,
	       fieldset bits,
	       unsigned ns) falcON_THROWING :
  BITS      ( bits ),
  C_FORTRAN ( 0 )
{
  debug_info(3,"bodies::bodies(): constructing bodies: nb=%u, ns=%u, bits=%s",
	     nb,ns,word(bits));
  unsigned n[BT_NUM] = {ns, nb>ns? nb-ns:0};
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
bodies::bodies(const unsigned*nall,
	       fieldset       bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 )
{
  unsigned none[BT_NUM]={0,0};
  const unsigned*n = nall? nall : none;
  debug_info(3,"bodies::bodies(): constructing bodies: n=%u,%u, bits=%s",
	     n[0],n[1],word(bits));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
// resets N, data; same as destruction followed by constructor 1            
void bodies::reset(const unsigned*nall,
		   fieldset       bits) falcON_THROWING
{
  unsigned none[BT_NUM]={0,0};
  const unsigned*n= nall? nall : none;
  bool keepN = true;
  for(bodytype t; t; ++t) keepN = keepN && NALL[t] == n[t];
  if(keepN) {
    NTOT = 0u;
    for(bodytype t; t; ++t) NTOT += (NBOD[t] = NALL[t]);
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
    falcON_ExceptF("copyflag !=0, but other bodies not supporting flag",
		   "bodies::bodies()");
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
: BITS      ( 0 ),
  C_FORTRAN ( 1 )
{
  debug_info(3,"bodies::bodies(): constructing bodies for C & FORTRAN: n=%u,%u",
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
  debug_info(2,"bodies::create(): created %u new bodies of type %s\n",
	     N,t.name());
}
////////////////////////////////////////////////////////////////////////////////
bodies::iterator bodies::new_body(bodytype t) falcON_THROWING
{
  if(0 == N_free(t)) {
    warning("bodies::new_body(): no body available\n");
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
  if(Nread == 0 || Nread > snap.Nbod()) Nread = snap.Nbod();
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
    debug_info(2,"bodies::read_snapshot(): phases read");
    read |= get & fieldset::phases;
    BITS |= get & fieldset::phases;
  }
  // now read data field by field
  for(fieldbit f; f; ++f) if(get.contain(f)) {
    debug_info(6,"bodies::read_snapshot(): f=%c: %s\n",letter(f),
	       read.contain(f)? "already read" :
	       !snap.has(nemo_io::field(f))? "not present" : "to be read");
    if( !read.contain(f) && snap.has(nemo_io::field(f)) ) {
      data_in inpt(snap,nemo_io::field(f));
      body b(start);
      b.read_data(inpt,Nread);
      if(inpt.N() != inpt.N_read())
	falcON_THROW("bodies::read_snapshot(): "
		     "could only read %u of %u %c data",
		     inpt.N_read(), inpt.N(), letter(f));
      debug_info(2,"bodies::read_snapshot(): %u %c read",
		 inpt.N_read(), letter(f));
      BITS |= fieldset(f);
      read |= fieldset(f);
    }
  }
  debug_info(1,"bodies::read_snapshot(): read=%s\n",word(read));
  if(read & fieldset::source) mark_srce_data_changed();
  if(read & fieldset::sphmax) mark_sph_data_changed();
  if(warn && want != read) warning("bodies::read_snapshot: couldn't read %s",
				  word(want.missing(read)));
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
  if(Nwrite == 0 || Nwrite > snap.Nbod()) Nwrite = snap.Nbod();
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
    debug_info(2,"bodies::write_snapshot(): written pq");
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
      debug_info(2,"bodies::write_snapshot(): written %u %c",
		 outp.N_written(), letter(f));
      written |= fieldset(f);
    }
  debug_info(1,"bodies::write_snapshot(): "
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
  void Ignore(std::istream&, body &) {}
  typedef void (*p_reader)(std::istream&, body &);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::read_simple_ascii(std::istream  &in,
			       const fieldbit*item,
			       unsigned       Ni,
			       unsigned       Nb,
			       unsigned       Ns)
{
  // 1. create table of readers                                                 
  fieldset get;
  p_reader readsph[100]={0}, readstd[100]={0};
  if(Ni > 100) {
    Ni = 100;
    warning(" can only read the first 100 data entries\n");
  }
  for(int i=0; i!=Ni; ++i) {
    if(get.contain(item[i]))
      warning("bodies::read_simple_ascii: reading item more than once");
    get |= fieldset(item[i]);
    switch(value(item[i])) {
#define SET_READ_STD(BIT,NAME)				\
    case BIT:						\
      readsph[i]= &Read<BIT>;				\
      readstd[i]= &Read<BIT>;				\
      break;
      DEF_NAMED_NONSPH(SET_READ_STD);
#undef SET_READ_STD
#define SET_READ_SPH(BIT,NAME)				\
    case BIT:						\
      readsph[i]= &Read<BIT>;				\
      readstd[i]= 0;					\
      break;
      DEF_NAMED_SPH(SET_READ_SPH);
#undef SET_READ_SPH
    default: 
      readsph[i]= 0;
      readstd[i]= 0;
      break;
    }
  }
  if(Ni==100)
    warning("bodies::read_simple_ascii(): cannot read >100 items\n");
  // 2. reset N & add fields                                                    
  resetN(Nb,Ns);
  add_fields(get);
  // 3. loop bodies & read data whereby ignoring lines starting with '#'        
  LoopSPHBodies(this,Bi) {
    while( in && eat_line(in,'#') );
    if(!in) falcON_ExceptF("end of input before data have been read",
			   "bodies::read_simple_ascii()");
    for(int i=0; i!=Ni; ++i)
      if(readsph[i]) readsph[i](in,Bi);
    SwallowRestofLine(in);
  }
  LoopSTDBodies(this,Bi) {
    while( in && eat_line(in,'#') );
    if(!in) falcON_ExceptF("end of input before data have been read",
			   "bodies::read_simple_ascii()");
    for(int i=0; i!=Ni; ++i)
      if(readstd[i]) readstd[i](in,Bi);
    SwallowRestofLine(in);
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
	warning("snapshot::del_pointer()"
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
  debug_info(4,"snapshot::add_pointer() %p to '%s' under \"%s\"\n",p,n,k);
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
  debug_info(4,"snapshot::set_pointer() %p to '%s' under \"%s\"\n",p,n,k);
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
  debug_info(4,"snapshot::get_pointer() %p to '%s' under \"%s\"\n",p,n,k);
  return p;
}
////////////////////////////////////////////////////////////////////////////////
void snapshot::del_pointer(const char*k) const
{
  debug_info(4,"snapshot::del_pointer() under \"%s\"\n",k);
  if(PBNK) static_cast<PointerBank*>(PBNK)->del(k);
}
////////////////////////////////////////////////////////////////////////////////
snapshot::snapshot(fieldset Bd) falcON_THROWING
: bodies ( static_cast<const unsigned*>(0), Bd ),
  INIT   ( false ),
  TINI   ( 0. ),
  TIME   ( 0. ),
  PBNK   ( 0 ) {}
////////////////////////////////////////////////////////////////////////////////
snapshot::snapshot(double   t,
		   unsigned Nb,
		   fieldset Bd,
		   unsigned Ns) falcON_THROWING
: bodies ( Nb,Bd,Ns ),
  INIT   ( true ),
  TINI   ( t ),
  TIME   ( t ),
  PBNK   ( 0 ) {}
////////////////////////////////////////////////////////////////////////////////
snapshot::snapshot(double         t,
		   const unsigned*N,
		   fieldset       Bd) falcON_THROWING
: bodies ( N,Bd ),
  INIT   ( true ),
  TINI   ( t ),
  TIME   ( t ),
  PBNK   ( 0 ) {}
////////////////////////////////////////////////////////////////////////////////
snapshot::snapshot(double       t,
		   bodies const&B,
		   fieldset     Bd,
		   flags        F) falcON_THROWING
: bodies ( B,Bd,F ),
  INIT   ( true ),
  TINI   ( t ),
  TIME   ( t ),
  PBNK   ( 0 ) {}
////////////////////////////////////////////////////////////////////////////////
snapshot::snapshot(snapshot const&S,
		   fieldset       Bd,
		   flags          F) falcON_THROWING
: bodies ( S,Bd,F ),
  INIT   ( S.INIT ),
  TINI   ( S.TINI ),
  TIME   ( S.TIME ),
  PBNK   ( S.PBNK? new PointerBank(*(static_cast<PointerBank*>(S.PBNK))) : 0 )
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
      r = fieldset(0);
      return false;
    }
    TIME = s.time();
    if(!INIT) { 
      TINI = TIME;
      INIT = true;
    }
  } else
    TIME = 0.;
  if(s.Nbod() != N_bodies() ||
     s.Nsph() != N_bodies(bodytype::gas) )
    reset(s.Nbod(), fieldset::o, s.Nsph());
  r = read_snapshot(s,g,begin_all_bodies(),N_bodies(),w);
  return true;
}
////////////////////////////////////////////////////////////////////////////////
fieldset snapshot::read_nemo(                      // R: what has been read     
			     snap_in  const&s,     // I: nemo input             
			     fieldset       g,     // I: what to read           
			     iterator const&b,     // I: start position         
			     unsigned       n,     //[I: #, def: all in input]  
			     bool           w)     //[I: warn: missing data]    
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
    warning("snapshot::write_nemo() cannot write %u bodies, "
	    "will only write %u\n",n,N_bodies()-i);
    n = N_bodies()-i;
  }
  unsigned ns = N_bodies(bodytype::gas)>i? N_bodies(bodytype::gas)-i : 0;
  snap_out s(o,n,ns,TIME);
  write_snapshot(s,w,b,n);
}
////////////////////////////////////////////////////////////////////////////////
void snapshot::write_nemo(nemo_out const&o,        // I: nemo output            
			  fieldset       w) const  // I: what to write          
  falcON_THROWING
{
  snap_out s(o,N_bodies(),N_bodies(bodytype::gas),TIME);
  write_snapshot(s,w,begin_all_bodies(),N_bodies());
}
#endif // falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
namespace {
  struct GadgetHeader {  // structure taken from gadget/allvars.h
    int          npart[6];
    double       masstab[6];
    double       time;
    double       redshift;
    int          flag_sfr;
    int          flag_feedback;
    unsigned int npartTotal[6];
    int          flag_cooling;
    int          num_files;
    double       BoxSize;
    double       Omega0;
    double       OmegaLambda;
    double       HubbleParam;
    int          flag_stellarage;
    int          flag_metals;
    unsigned int npartTotalHighWord[6];
    int          flag_entropy_instead_u;
    char         fill[60];
    //--------------------------------------------------------------------------
    GadgetHeader() :
      time(0.), redshift(0.), flag_sfr(0), flag_feedback(0), flag_cooling(0),
      num_files(0), BoxSize(0.), Omega0(0.), OmegaLambda(0.), HubbleParam(0.),
      flag_stellarage(0), flag_metals(0), flag_entropy_instead_u(0)
    {
      for(int k=0; k!=6; ++k) {
	npart[k] = 0;
	npartTotal[k] = 0;
	npartTotalHighWord[k] = 0;
	masstab[k] = 0.;
      }
    }
    //--------------------------------------------------------------------------
    bool mismatch(GadgetHeader const&H) const {
      bool okay = true;
#define CHECK_I(FIELD,FIELDNAME)				\
      if(FIELD != H.FIELD) {					\
	okay = false;						\
	warning("GadgetHeader \"%s\" mismatch (%u vs %u)\n",	\
		FIELDNAME, FIELD, H.FIELD);			\
      }
#define CHECK_D(FIELD,FIELDNAME)				\
      if(FIELD != H.FIELD) {					\
	okay = false;						\
	warning("GadgetHeader\"%s\" mismatch (%f vs %f)\n",	\
		FIELDNAME, FIELD, H.FIELD);			\
      }
      CHECK_D(time,"time");
      CHECK_D(redshift,"redshift");
      CHECK_I(flag_sfr,"flag_sfr");
      CHECK_I(flag_feedback,"flag_feedback");
      CHECK_I(npartTotal[0],"npartTotal[0]");
      CHECK_I(npartTotal[1],"npartTotal[1]");
      CHECK_I(npartTotal[2],"npartTotal[2]");
      CHECK_I(npartTotal[3],"npartTotal[3]");
      CHECK_I(npartTotal[4],"npartTotal[4]");
      CHECK_I(npartTotal[5],"npartTotal[5]");
      CHECK_I(flag_cooling,"flag_cooling");
      CHECK_I(num_files,"num_files");
      CHECK_D(BoxSize,"BoxSize");
      CHECK_D(Omega0,"Omega0");
      CHECK_D(OmegaLambda,"OmegaLambda");
      CHECK_D(HubbleParam,"HubbleParam");
      CHECK_I(flag_stellarage,"flag_stellarage");
      CHECK_I(flag_metals,"flag_metals");
      CHECK_I(flag_entropy_instead_u,"flag_entropy_instead_u");
      CHECK_I(npartTotalHighWord[0],"npartTotalHighWord[0]");
      CHECK_I(npartTotalHighWord[1],"npartTotalHighWord[1]");
      CHECK_I(npartTotalHighWord[2],"npartTotalHighWord[2]");
      CHECK_I(npartTotalHighWord[3],"npartTotalHighWord[3]");
      CHECK_I(npartTotalHighWord[4],"npartTotalHighWord[4]");
      CHECK_I(npartTotalHighWord[5],"npartTotalHighWord[5]");
      return !okay;
#undef CHECK_I
#undef CHECK_D
    }
    //--------------------------------------------------------------------------
    void check_simple_npart_error() {
      for(int k=0; k!=6; ++k)
	if(npart[k] > npartTotal[k]) {
	  warning("GadgetHeader: npart[%u]=%u > npartTotal[%u]=%u: "
		  "we will try to fix by setting npartTotal[%u]=%u\n",
		  k,npart[k],k,npartTotal[k],k,npart[k]);
	  npartTotal[k] = npart[k];
	}
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream&out) const {
      out<<" gadget header dump:";
      for(int k=0; k!=6; ++k)
	out<<"\n type "<<k
	   <<": npart="<<std::setw(8)<<npart[k]
	   <<" npartTotal="<<std::setw(8)<<npartTotal[k]
	   <<" masstab="<<masstab[k];
      out<<"\n redshift               = "<<redshift
	 <<"\n flag_sfr               = "<<flag_sfr
	 <<"\n flag_feedback          = "<<flag_feedback
	 <<"\n flag_cooling           = "<<flag_cooling
	 <<"\n num_files              = "<<num_files
	 <<"\n BoxSize                = "<<BoxSize
	 <<"\n Omega0                 = "<<Omega0
	 <<"\n OmegaLambda            = "<<OmegaLambda
	 <<"\n HubbleParam            = "<<HubbleParam
	 <<"\n flag_stellarage        = "<<flag_stellarage
	 <<"\n flag_metals            = "<<flag_metals
	 <<"\n flag_entropy_instead_u = "<<flag_entropy_instead_u
	 <<std::endl;
    }
  };
} // namespace {
//------------------------------------------------------------------------------
falcON_TRAITS(GadgetHeader,"GadgetHeader");
////////////////////////////////////////////////////////////////////////////////
#define READ(BIT) if(!is_sph(BIT) && nd || ns) {			\
    FortranIRec F(in, rec);						\
    unsigned nr = is_sph(BIT)? ns : ns+nd;				\
    if(read.contain(BIT)) {						\
      if(F.size() != nr*field_traits<BIT>::size)			\
	falcON_THROW("bodies::read_gadget(): mismatch reading %u %c: "	\
		     "expected %u bytes, found %u\n",			\
		     nr,field_traits<BIT>::word(),			\
		     nr*field_traits<BIT>::size,F.size());		\
      add_field(BIT);							\
      if(ns) {								\
	body sph(SPH);							\
	sph.read_Fortran(F, BIT, ns);					\
      }									\
      if(!is_sph(BIT) && nd) {						\
	body std(STD);							\
        std.read_Fortran(F, BIT, nd);					\
      }									\
      debug_info(2,"bodies::read_gadget(): read %u %c\n",		\
	         nr, field_traits<BIT>::word());			\
      fgot |= fieldset(BIT);						\
    } else {								\
      F.skip_bytes(F.size());						\
      debug_info(3,"bodies::read_gadget(): skip %u %c\n",		\
                 nr, field_traits<BIT>::word());			\
    }									\
  }
//------------------------------------------------------------------------------
double bodies::read_gadget(const char*fname,
			   fieldset   read,
			   unsigned   rec) falcON_THROWING
{
  read &= fieldset("mxvkURHpa");
  fieldset got;
  // 1 open first data file, read header, and determine number of data files
  GadgetHeader header0,headeri,*header=&header0;
  input        in(fname);
  int          nfile=1;
  char         filename[256];
  const char  *file=0;
  if(in) {
    // 1.1 try single file "fname"
    file = fname;
    try { FortranIRec::Read(in, header, 1, rec); }
    catch(exception E) { falcON_RETHROW(E); }
    nfile = header->num_files;
    if(nfile==0) nfile=1;
    if(debug(1)) {
      debug_info("bodies::read_gadget(): header read from file \"%s\":\n",
		 fname);
      header->dump(std::clog);
    }
    if(nfile!=1)
      falcON_THROW("bodies::read_gadget(): num_files=%u, expected 0 or 1\n",
		   nfile);
    header->check_simple_npart_error();
  } else {
    // 1.2 try file "fname.0" (with possibly more to follow)
    sprintf(filename,"%s.%u",fname,0);
    in.open(filename);
    if(!in) falcON_THROW("bodies::read_gadget(): cannot open file \"%s\" "
			 "nor file \"%s\"\n", fname, filename);
    file = filename;
    try { FortranIRec::Read(in, header, 1, rec); }
    catch(exception E) { falcON_RETHROW(E); }
    nfile = header->num_files;
    if(debug(1)) {
      debug_info("bodies::read_gadget(): header read from file \"%s\":\n",file);
      header->dump(std::clog);
    }
  }
  // 2 establish number of SPH and non-SPH particles and allocate memory.
  unsigned NB[BT_NUM] = {0}, NP[6] ={0};
  NB[0] = header->npartTotal[0];
  for(int k=1; k!=6; ++k) NB[1] += header->npartTotal[k];
  if(NB[0] == 0) read &= fieldset(fieldset::nonSPH);
  reset(NB,read);
  // 3 loop data files
  body SPH = begin_sph_bodies();
  body STD = begin_std_bodies();
  for(int ifile=0; ifile!=nfile; ) {
    fieldset fgot;
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
    // read positions
    READ(fieldbit::x);
    if(read == fgot) goto NextFile;
    // read velocties
    READ(fieldbit::v);
    if(read == fgot) goto NextFile;
    // read keys
    READ(fieldbit::k);
    if(read == fgot) goto NextFile;
    // read masses --- OR assign them ...
    if(nm) {
      FortranIRec F(in, rec);
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
	    sph.read_Fortran(F, fieldbit::m, header->npart[0]);
	}
	for(int k=1; k!=6; ++k) if(header->npart[k]) {
	  if(header->masstab[k])
	    for(int b=0; b!=header->npart[k]; ++b,++std)
	      std.mass() = header->masstab[k];
	  else
	    std.read_Fortran(F, fieldbit::m, header->npart[k]);
	}
	debug_info(2,"bodies::read_gadget(): read %u m\n", nm);
	fgot |= fieldset(fieldbit::m);
      } else {
	F.skip_bytes(F.size());
	debug_info(3,"bodies::read_gadget(): skip %u m\n", nm);
      }
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
    // read gas internal energies
    READ(fieldbit::U);
    if(read == fgot) goto NextFile;
    // read gas densities
    READ(fieldbit::R);
    if(read == fgot) goto NextFile;
    // read SPH smoothing lengths
    READ(fieldbit::H);
    if(read == fgot) goto NextFile;
    // read potentials
    READ(fieldbit::p);
    if(read == fgot) goto NextFile;
    // read accelerations
    READ(fieldbit::a);
    if(read == fgot) goto NextFile;
  NextFile:
    got |= fgot;
    // open next data file (unless all have been read) and move SPH & STD
    if(++ifile != nfile) {
      SPH += ns;
      STD += nd;
      sprintf(filename,"%s.%u",fname,ifile);
      in.open(file);
      if(!in) falcON_THROW("bodies::read_gadget(): cannot open file \"%s\"\n",
			   file);
      try { FortranIRec::Read(in, &headeri, 1, rec); }
      catch(exception E) { falcON_RETHROW(E); }
      if(header0.mismatch(headeri))
	falcON_THROW("bodies::read_gadget(): header mismatch\n");
      header = &headeri;
      if(debug(2)) {
	debug_info("bodies::read_gadget(): header read from file \"%s\":\n",
		   file);
	header->dump(std::clog);
      }
    }
  }
  debug_info(1,"bodies::read_gadget(): read %s for %u SPH & %u STD bodies\n",
	     word(got), NB[0], NB[1]);
  return header->time;
}
////////////////////////////////////////////////////////////////////////////////
#define WRITE(BIT)							\
  if(!is_sph(BIT) || N_sph()) {						\
    unsigned nw = is_sph(BIT)? N_sph() : N_bodies();			\
    FortranORec F(out, nw * field_traits<BIT>::size, rec);		\
    if(have(BIT)) {							\
      if(N_sph())							\
	begin_sph_bodies().write_Fortran(F, BIT, N_sph());		\
      if(!is_sph(BIT) && N_std())					\
        begin_std_bodies().write_Fortran(F, BIT, N_std());		\
      debug_info(2,"bodies::write_gadget(): written %u %c\n",		\
	         nw, field_traits<BIT>::word());			\
    } else {								\
      if(warn)								\
	warning("bodies::write_gadget(): don't have %c, write out zeros\n", \
		field_traits<BIT>::word());				\
      F.fill_bytes(nw);							\
      debug_info(2,"bodies::write_gadget(): written %u 0 for %c\n",	\
	         nw, field_traits<BIT>::word());			\
    }									\
    written |= fieldset(BIT);						\
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
  WRITE(fieldbit::x);
  WRITE(fieldbit::v);
  WRITE(fieldbit::k);
  WRITE(fieldbit::m);
  WRITE(fieldbit::U);
  if(write & fieldset("RHpa")) {
  WRITE(fieldbit::R);
  if(write & fieldset("Hpa")) {
  WRITE(fieldbit::H);
  if(write & fieldset("pa")) {
  WRITE(fieldbit::p);
  if(write & fieldset::a)
  WRITE(fieldbit::a);
  } } }
  debug_info(1,"bodies::write_gadget(): written %s for "
	     "%u SPH & %u STD bodies\n", word(written), N_sph(), N_std());
  // make sure we have no keys if we didn't have them before ...
  if(!had_keys) const_cast<bodies*>(this)->del_field(fieldbit::k);
}
////////////////////////////////////////////////////////////////////////////////
