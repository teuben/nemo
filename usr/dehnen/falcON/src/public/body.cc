// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.cc                                                                     |
//                                                                             |
// Copyright (C) 2000-2005 Walter Dehnen                                       |
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
#include <public/inline_io.h>                      // utilities for C++ I/O     
#include <public/numerics.h>                       // falcON numeric utilities  

#ifdef falcON_NEMO                                 // compiler option           
  extern "C" {
#   include <stdinc.h>                             // NEMO basics               
  }
#endif

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  falcON_TRAITS(bodies::block,"bodies::block","bodies::blocks");
}
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
    debug_info(4,"allocating data for %d %c (%s)\n",NALL,letter(f),name(f));
    set_data_void(f, falcON_NEW(char,NALL*falcON::size(f)));
    if(f == fieldbit::f) reset_flags();
  }
}
////////////////////////////////////////////////////////////////////////////////
inline
void bodies::block::del_field (fieldbit f) falcON_THROWING {
  if(DATA[value(f)]) {
    debug_info(4,"de-allocating data for %c (%s)\n",letter(f),name(f));
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
template<unsigned BIT=0, unsigned END=BD_NQUANT> struct CopyBody {
  static const unsigned BD = 1<<BIT;
  static void copy(void    **data,
		   unsigned  from,
		   unsigned  to  ,
		   fieldset  b,
		   fieldset &c) {
    if(data[BIT] && b.contain(fieldbit(BIT)) ) {
      memcpy(static_cast<      char*>(data[BIT])+to  *BD_ZQUANT[BIT],
	     static_cast<const char*>(data[BIT])+from*BD_ZQUANT[BIT],
	     BD_ZQUANT[BIT]);
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
template<unsigned BIT=0, unsigned END=BD_NQUANT> struct CopyBodies {
  static const unsigned BD = 1<<BIT;
  static void copy(void*const*data_from,
		   void*const*data_to,
		   unsigned  from,
		   unsigned  to,
		   unsigned  num,
		   fieldset  b,
		   fieldset&c) {
    if(data_from[BIT] && data_to[BIT] && b & fieldset(1<<BIT) ) {
      memcpy(static_cast<      char*>(data_from[BIT])+to  *BD_ZQUANT[BIT],
	     static_cast<const char*>(data_to  [BIT])+from*BD_ZQUANT[BIT],
	     num*BD_ZQUANT[BIT]);
      c |= fieldset(1<<BIT);
    }
    CopyBodies<BIT+1, END>::copy(data_from,data_to,from,to,num,b,c);
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
//                                                                              
// class falcON::bodies                                                         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
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
  debug_info(6,"destructing bodies");
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
void bodies::set_data(unsigned *N) falcON_THROWING
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
  debug_info(3,"constructing bodies: nb=%d, ns=%d, bits=%s",nb,ns,word(bits));
  unsigned n[BT_NUM] = {ns, nb>ns? nb-ns:0};
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
bodies::bodies(unsigned*nall,
	       fieldset bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 )
{
  unsigned none[BT_NUM]={0,0}, *n = nall? nall : none;
  debug_info(3,"constructing bodies: n=%d,%d, bits=%s",n[0],n[1],word(bits));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
////////////////////////////////////////////////////////////////////////////////
// resets N, data; same as destruction followed by constructor 1            
void bodies::reset(unsigned*nall,
		   fieldset bits) falcON_THROWING
{
  unsigned none[BT_NUM]={0,0}, *n= nall? nall : none;
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
  C_FORTRAN ( 1 ),
  NBLK      ( 0 ),
  NTOT      ( 0 )
{
  block **LAST=&FIRST;
  for(bodytype t; t; ++t)
    if(n[t] > index::max_bodies)
      falcON_THROW("too many bodies\n");
    else if(n[t] > 0) {
      TYPES[t]      = new block(NBLK,n[t],n[t],NTOT,t,fieldset::o,this);
      NALL [t]      = n[t];
      NBOD [t]      = n[t];
      BLOCK[NBLK++] = TYPES[t];
      *LAST         = TYPES[t];
      LAST          = &(TYPES[t]->NEXT);
      NTOT         += n[t];
    } else {
      NALL [t]      = 0;
      NBOD [t]      = 0;
      TYPES[t]      = 0;
    }
  *LAST = 0;
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
  for(const block*p=FIRST; p; p=p->next())
    const_cast<block*>(p)->add_field(f);
  BITS |= fieldset(f);
}
////////////////////////////////////////////////////////////////////////////////
void bodies::add_fields(fieldset b) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->add_fields(b);
  BITS |= b;
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
    falcON_THROW("bodies::create(): asked for %d > %d bodies\n",
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
  debug_info(2,"bodies::create(): created %d new bodies of type %s\n",
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
		     "could only read %d of %d %c data",
		     inpt.N_read(), inpt.N(), letter(f));
      debug_info(2,"bodies::read_snapshot(): %d %c read",
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
    if(put.contain(f) && !written.contain(f) ) {
      data_out outp(snap,nemo_io::field(f));
      body b(start);
      b.write_data(outp,Nwrite);
      if(outp.N() != outp.N_written())
	falcON_THROW("bodies::write_snapshot(): "
		     "could only write %d of %d %c data",
		     outp.N_written(), outp.N(), letter(f));
      debug_info(2,"bodies::write_snapshot(): written %d %c",
		 outp.N_written(), letter(f));
      written |= fieldset(f);
    }
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
		    real       (*func)(iterator const&),
		    iterator     b,
		    unsigned     n) const falcON_THROWING
{
  if(this != b.my_bodies())
    falcON_THROW("bodies::sorted(): body not from this set of bodies");
  if(n == 0) n = N_bodies() - bodyindex(b);
  else if(bodyindex(b) + n > N_bodies()) {
    warning("bodies::sorted(): can only sort %d instead of %d bodies",
	    N_bodies()-bodyindex(b),n);
    n = N_bodies()-bodyindex(b);
  }
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  for(int i=0; i!=n; ++b,++i) {
    I[i] = static_cast<index>(b);
    Q[i] = func(b);
  }
  int  *R = falcON_NEW(int,n);
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
		    real       (*func)(iterator const&),
		    iterator     b,
		    unsigned     n) const falcON_THROWING
{
  if(this != b.my_bodies())
    falcON_THROW("bodies::sorted(): body not from this set of bodies");
  if(n == 0) n = N_bodies() - bodyindex(b);
  else if(bodyindex(b) + n > N_bodies()) {
    warning("bodies::sorted(): can only sort %d instead of %d bodies",
	    N_bodies()-bodyindex(b),n);
    n = N_bodies()-bodyindex(b);
  }
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  for(int i=0; i!=n; ++b,++i) {
    I[i] = static_cast<index>(b);
    Q[i] = func(b);
  }
  int  *R = falcON_NEW(int,n);
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
    warning("snapshot::write_nemo() cannot write %d bodies, "
	    "will only write %d\n",n,N_bodies()-i);
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
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
