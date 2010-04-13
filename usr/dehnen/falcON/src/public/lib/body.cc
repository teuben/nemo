// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/lib/body.cc
///
/// \author  Walter Dehnen
///
/// \date    2000-2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2010  Walter Dehnen
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
#include <body.h>                                  // falcON::bodies etc        
#include <iostream>                                // C++ basic I/O             
#include <fstream>                                 // C++ file I/O              
#include <sstream>                                 // C++ string I/O            
#include <iomanip>                                 // C++ I/O formating         
#include <cstring>                                 // C++ strings               
#include <public/nemo++.h>                         // utilities for NEMO I/O    
#include <public/bodyfunc.h>
#include <utils/numerics.h>
#include <utils/heap.h>

namespace falcON {
  using namespace WDutils;
}

using namespace falcON;
typedef long unsigned lu; // will convert size_t to this type in printf

falcON_TRAITS(falcON::bodies::block,"bodies::block");
//                                                                              
// struct falcON::bodies::block                                                 
//                                                                              
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
  this->NALL       = that->NALL;
  this->NBOD       = that->NBOD;
  this->FIRST      = that->FIRST;
  this->LOCALFIRST = that->LOCALFIRST;
}
//
void bodies::block::reset_flags() const
{
  if(0 != DATA[fieldbit::f]) {
    if(TYPE.is_sph()) 
      for(unsigned n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sph;
    else
    if(TYPE.is_sink()) 
      for(unsigned n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sink;
    else
      for(unsigned n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::empty;
  }
}
//
void bodies::block::flag_all_as_active() const falcON_THROWING
{
  if(0 != DATA[fieldbit::f])
    for(unsigned n=0; n!=NALL; ++n)
      datum<fieldbit::f>(n).add(flags::active);
  else 
    falcON_THROW("in bodies::flag_all_as_active(): flags not supported");
}
//
void bodies::block::reset_data(fieldset b) const falcON_THROWING
{
#define RESETDATA(BIT,NAME)				\
  if(DATA[BIT] && b.contain(BIT) && NBOD)		\
    for(unsigned n=0; n!=NBOD; ++n)			\
      field_traits<BIT>::set_zero(datum<BIT>(n));
 DEF_NAMED(RESETDATA)
#undef RESETDATA
}
//
void bodies::block::add_field (fieldbit f) falcON_THROWING
{
  if(TYPE.allows(f) && 0 == DATA[value(f)] ) {
    DebugInfo(4,"bodies::block::add_field(): "
	      "allocating data for %s bodies: %u %c (%s)\n",
	      TYPE.name(),NALL,letter(f),fullname(f));
    set_data_void(f, falcON_NEW(char,NALL*falcON::size(f)));
    if(f == fieldbit::f) reset_flags();
  }
}
//
void bodies::block::del_field (fieldbit f) falcON_THROWING
{
  if(DATA[value(f)]) {
    DebugInfo(4,"bodies::block::del_field(): "
	      "de-allocating data for %s bodies: %c (%s)\n",
	      TYPE.name(),letter(f),fullname(f));
    falcON_DEL_A(static_cast<char*>(DATA[value(f)]));
  }
  set_data_void(f,0);
}
//
void bodies::block::swap_bytes(fieldbit f) falcON_THROWING
{
  if(DATA[value(f)]) {
    DebugInfo(4,"bodies::block::swap_bytes(): swapping bytes for %c (%s)\n",
	      letter(f),fullname(f));
    falcON::swap_bytes(DATA[value(f)], falcON::size(f), NALL);
  }
}
//
void bodies::block::add_fields(fieldset b) falcON_THROWING
{
  for(fieldbit f; f; ++f)
    if(b.contain(f)) add_field(f);
}
//
void bodies::block::del_fields(fieldset b) falcON_THROWING
{
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) del_field(f);
}
//
void bodies::block::set_fields(fieldset b) falcON_THROWING
{
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) add_field(f);
    else             del_field(f);
}
//
bodies::block::~block() falcON_THROWING
{
  for(fieldbit f; f; ++f)
    del_field(f);
}
//
bodies::block::block(unsigned no,                  // I: our No                 
		     unsigned na,                  // I: data to allocate       
		     unsigned nb,                  // I: # bodies <= na_b       
		     unsigned fst,                 // I: first body index       
		     bodytype typ,                 // I: hold sph bodies?       
		     fieldset bits,                // I: data to allocate       
		     bodies  *bods)                // I: pointer to my bodies   
  falcON_THROWING
: TYPE       ( typ ),
  NALL       ( na ), 
  NBOD       ( nb ), 
  NO         ( no ),
  FIRST      ( fst ),
  LOCALFIRST ( fst ),
  NEXT       ( 0 ),
  BODS       ( bods )
{
  if(na<nb)
    falcON_THROW("in bodies::block::block(): N_alloc < N_bodies");
  DebugInfo(6,"bodies::block: na=%d, bits=%s, type=%s allowed bits=%s\n",
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
//
fieldset bodies::block::copy_body(unsigned fr, unsigned to, fieldset _copy)
{
  if(fr>=NALL)
    falcON_THROW("in bodies::block::copy_body(): "
		 "from=%d > NALL=%d\n", fr,NALL);
  if(to>=NALL)
    falcON_THROW("in bodies::block::copy_body(): "
		 "to=%d > NALL=%d\n", to,NALL);
  fieldset copied(fieldset::empty);
  if(fr!=to) {
    for(fieldbit f; f; ++f)
      if(_copy.contain(f) && data_void(f)) {
	memcpy(static_cast<      char*>(data_void(f))+to*falcON::size(f),
	       static_cast<const char*>(data_void(f))+fr*falcON::size(f),
	       falcON::size(f));
	copied |= fieldset(f);
      }
    DebugInfo(8,"bodies::block::copy_body(): copied %s from %d to %d\n",
	      word(copied),fr,to);
  }
  return copied;
}
//
fieldset bodies::block::copy_bodies(const block*that,
				    unsigned    fr,
				    unsigned    to,
				    unsigned    n,
				    fieldset   _copy) falcON_THROWING
{
  if(this == that)
    falcON_THROW("in bodies::block::copy_bodies() from same block");
  if(to+n > NALL)
    falcON_THROW("in bodies::block::copy_bodies(): "
		 "to+n=%d > NALL=%d\n", to+n,NALL);
  if(fr+n > that->NALL)
    falcON_THROW("in bodies::block::copy_bodies(): "
		 "from+n=%d > that->NALL=%d\n", fr+n,that->NALL);
  fieldset copied(fieldset::empty);
  for(fieldbit f; f; ++f)
    if(_copy.contain(f) && data_void(f) && that->data_void(f)) {
      memcpy(static_cast<      char*>(this->data_void(f))+to*falcON::size(f),
	     static_cast<const char*>(that->data_void(f))+fr*falcON::size(f),
	     n*falcON::size(f));
      copied |= fieldset(f);
    }
  return copied;
}
//
inline void bodies::block::skip(unsigned&from,
				flags    copyflag) const falcON_THROWING
{
  if(copyflag)
    for(; from<NBOD && !(flag(from).are_set(copyflag)); ++from ) {}
}
//
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
  if(From == 0) return fieldset::empty;
  unsigned _copy;
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
      _copy = 0u;
      for(unsigned to=from;
	  to < From->NBOD && From->flag(to).are_set(copyflag) && _copy<free;
	  ++_copy, ++to) {}
    } else
      _copy = min(free, From->NBOD - from);
    // if any body to be copied, copy data, adjust free, NBOD, from, copied     
    if(_copy) {
      fieldset c = copy_bodies(From, from, NBOD, _copy, copydata);
      free -= _copy;
      NBOD += _copy;
      from += _copy;
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
//
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
  DebugInfo(6,"bodies::block::remove(): removed %d: NBOD=%d\n",removed,NBOD);
}
#ifdef falcON_NEMO
//
void bodies::block::read_data(data_in &inpt,
			      unsigned from,
			      unsigned N) falcON_THROWING
{
  fieldbit f= nemo_io::bit(inpt.field());
  if(!TYPE.allows(f))
    falcON_THROW("bodies::block::read_data(%c): not allowed by our type",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::read_data(%c): cannot read %d from %d "
		 "(NBOD=%d)\n",letter(f),N,from,NBOD);
  add_field(f);
  inpt.read(static_cast<char*>(DATA[value(f)])+from*falcON::size(f), N);
  DebugInfo(2,"bodies::block::read_data(): read %d %c",N,letter(f));
}
//
void bodies::block::read_posvel(data_in &inpt,
				unsigned from,
				unsigned N,
				fieldset want) falcON_THROWING
{
  if(inpt.field() != nemo_io::posvel)
    falcON_THROW("bodies::block::read_posvel(): input has not phases");
  if(from + N > NBOD)
    falcON_THROW("bodies::block::read_posvel(): cannot read %d from %d "
		 "(NBOD=%d)\n",N,from,NBOD);
  if(want.contain(fieldbit::x)) add_field(fieldbit::x);
  if(want.contain(fieldbit::v)) add_field(fieldbit::v);
  inpt.read_phases(want.contain(fieldbit::x)?
		   static_cast<vect*>(DATA[fieldbit::x]) + from : 0,
		   want.contain(fieldbit::v)?
		   static_cast<vect*>(DATA[fieldbit::v]) + from : 0, N);
  DebugInfo(2,"bodies::block::read_posvel(): read %d, %s",N,
	    word(want&fieldset::phases));
}
//
void bodies::block::write_data(data_out&outp,
			       unsigned from,
			       unsigned N) const falcON_THROWING
{
  fieldbit f= nemo_io::bit(outp.field());
  if(0 == DATA[value(f)])
    falcON_THROW("bodies::block::write_data(%c): data not supported",
		 letter(f));
  if(from + N > NBOD)
    falcON_THROW("bodies::block::write_data(%c): "
		 "cannot write %d from %d (NBOD=%d)",letter(f),N,from,NBOD);
  outp.write(static_cast<char*>(DATA[value(f)])+from*falcON::size(f), N);
}
//
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
  for(unsigned n=0,m=from; n!=N; ++n,++m)
    P[n] = const_datum<fieldbit::p>(m) + const_datum<fieldbit::q>(m);
  outp.write(P,N);
  falcON_DEL_A(P);
}
#endif // falcON_NEMO
//
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
  if(swap) {
    if(is_vector(f))
      falcON::swap_bytes(static_cast<void*>(C), sizeof(real), Ndim*N);
    else
      falcON::swap_bytes(static_cast<void*>(C), falcON::size(f), N);
  }
  if(R != N*falcON::size(f))
    falcON_THROW("bodies::block::read_Fortran(%c): "
		 "could only read %u of %lu bytes\n",
		 letter(f),R,lu(N*falcON::size(f)));
  DebugInfo(4,"bodies::block::read_Fortran(): read %u `%s'\n",N,fullname(f));
}
//
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
		 "could only write %u of %lu bytes\n",
		 letter(f),W,lu(N*falcON::size(f)));
  DebugInfo(4,"bodies::block::write_Fortran(): written %u `%s'\n",
	    N,fullname(f));
}
#endif
//
#ifdef falcON_NEMO
//                                                                              
// class falcON::bodies::iterator                                               
//                                                                              
bodies::iterator& bodies::iterator::read_data(data_in&D, unsigned R)
  falcON_THROWING
{
  if(R == 0 || R > D.N_unread()) R = D.N_unread();
  while(is_valid() && R) {
    unsigned r = min(B->N_bodies()-K, R);
    const_cast<block*>(B)->read_data(D,K,r);
    R -= r;
    K += r;
    if(K >= B->N_bodies()-K) next_block();
  }
  return *this;
}
//
bodies::iterator& bodies::iterator::write_data(data_out&D, unsigned W)
  falcON_THROWING
{
  if(W == 0 || W > D.N_free()) W = D.N_free();
  while(is_valid() && W) {
    unsigned w = min(B->N_bodies()-K, W);
    B->write_data(D,K,w);
    W -= w;
    K += w;
    if(K >= B->N_bodies()) next_block();
  }
  return *this;
}
//
bodies::iterator& bodies::iterator::write_potpex(data_out&D, unsigned W)
  falcON_THROWING
{
  if(W == 0 || W > D.N_free()) W = D.N_free();
  while(is_valid() && W) {
    unsigned w = min(B->N_bodies()-K, W);
    B->write_potpex(D,K,w);
    W -= w;
    K += w;
    if(K >= B->N_bodies()) next_block();
  }
  return *this;
}
//
bodies::iterator& bodies::iterator::read_posvel(data_in& D, fieldset get,
						unsigned R)
  falcON_THROWING
{
  if(R == 0 || R > D.N_unread()) R = D.N_unread();
  while(is_valid() && R) {
    unsigned r = min(B->N_bodies()-K, D.N_unread());
    const_cast<block*>(B)->read_posvel(D,K,r,get);
    R -= r;
    K += r;
    if(K >= B->N_bodies()) next_block();
  }
  return *this;
}
#endif // falcON_NEMO
//
#ifdef falcON_REAL_IS_FLOAT
bodies::iterator& bodies::iterator::read_Fortran(FortranIRec&I, fieldbit f, 
						 unsigned R, bool swap)
  falcON_THROWING
{
  if(R * falcON::size(f) > I.bytes_unread())
    falcON_THROW("body::read_Fortran(%c): want %u `%s' (%lu bytes) but "
		 "only %u bytes left on Fortran record\n",letter(f),
		 R, fullname(f), lu(R*falcON::size(f)), I.bytes_unread());
  while(is_valid() && R) {
    unsigned r = min(B->N_bodies()-K, R);
    const_cast<block*>(B)->read_Fortran(I,f,K,r,swap);
    R -= r;
    K += r;
    if(K >= B->N_bodies()) next_block();
  }
  if(R) falcON_THROW("body::read_Fortran: %u data remain unread\n",R);
  return *this;
}
//
bodies::iterator& bodies::iterator::write_Fortran(FortranORec&O,
						  fieldbit f, unsigned W)
  falcON_THROWING
{
  if(W * falcON::size(f) > O.bytes_free())
    falcON_THROW("body::write_Fortran(%c): want %u `%s' (%lu bytes) but "
		 "only %u bytes left free on Fortran record\n",letter(f),
		 W, fullname(f), lu(W*falcON::size(f)), O.bytes_free());
  while(is_valid() && W) {
    unsigned w = min(B->N_bodies()-K, W);
    B->write_Fortran(O,f,K,w);
    W -= w;
    K += w;
    if(K >= B->N_bodies()) next_block();
  }
  if(W) falcON_THROW("body::write_Fortran: %u data remain unwritten\n",W);
  return *this;
}
#endif
//                                                                              
// class falcON::bodies                                                         
//                                                                              

// link a new block in
void bodies::add_block(block*B)
{
  // link to last block of same or earlier type, if any, otherwise make it first
  block**P=&FIRST;
  while(*P && (*P)->TYPE <= B->TYPE) P = &((*P)->NEXT);
  B->link(*P);
  *P = B;
  // update TYPES[]
  if(0==TYPES[B->type()])
    TYPES[B->type()]=B;
  // update BLOCK[] and block::NO
  for(int I=0; I!=index::max_blocks; ++I)
    if(BLOCK[I] == 0) {
      BLOCK[I] = B;
      B->NO    = I;
      break;
    }
  // update NBLK and block::BODS
  B->BODS  = this;
  NBLK ++;
  // update block::FIRST and NALL[], NBOD[], NTOT
  set_firsts();
}
// erase a block from our linkage
void bodies::erase_block(block*B)
{
  if(B==0) return;
  // remove from FIRST
  if(FIRST == B)
    FIRST = B->next();
  // remove from TYPES[]
  if(TYPES[B->type()] == B)
    TYPES[B->type()] = B->next_of_same_type();
  // remove from block::NEXT
  for(int i=0; i!=index::max_blocks; ++i)
    if(BLOCK[i] && BLOCK[i]->next() == B) {
      BLOCK[i]->link(B->next());
      break;
    }
  // remove from BLOCK[]
  bool found = false;
  for(int i=0; i!=index::max_blocks; ++i)
    if(BLOCK[i] == B) {
      BLOCK[i] = 0;
      found = true;
      break;
    }
  // reset # bodies info and block::FIRST
  if(found) {
    --NBLK;
    B->BODS = 0;
    set_firsts();
  } else
    falcON_Warning("bodies::erase_block(): block not found in table\n");
}
// remove empty blocks
void bodies::remove_empty_blocks(bool all) falcON_THROWING
{
  for(;;) {
    block*B=0;
    // search for empty block
    for(int i=0; i!=index::max_blocks; ++i)
      if(BLOCK[i] && (all? BLOCK[i]->N_alloc():BLOCK[i]->N_bodies())==0) {
	B=BLOCK[i];
	break;
      }
    // found one: remove it
    if(B) {
      erase_block(B);
      falcON_DEL_O(B);
    } else
      break;
  }
}
// create a new block and link it in
bodies::block* bodies::new_block(bodytype t, unsigned Na, unsigned Nb,
				 fieldset f) falcON_THROWING
{
  if(Nb > Na)
    falcON_THROW("bodies::new_block(): Nb=%u > Na=%u\n",Nb,Na);
  if(Na > index::max_bodies) 
    falcON_THROW("bodies::new_block(): asked for %u > %u bodies\n",
		 Na,  index::max_bodies);
  if(NBLK >= index::max_blocks)
    falcON_THROW("bodies::new_block(): number of blocks exceeded\n");
  block*B=new block(0,Na,Nb,0,t,f,this);
  NNEW[t] += Nb;
  add_block(B);
  DebugInfo(2,"bodies::new_block(): created block for up to %u bodies "
	    "(%u active) of type %s\n", Na,Nb,t.name());
  return B;
}
// reset blocks' FIRST entries (used in parallel code)
void bodies::reset_firsts(unsigned fst[bodytype::NUM])
{
  for(bodytype t; t; ++t) {
    unsigned L=0;
    for(block*B=TYPES[t]; B; B=B->next_of_same_type()) {
      B->set_first(L+fst[t], L);
      L+=B->N_bodies();
    }
  }
}
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
  for(unsigned i=0; i!=index::max_blocks; ++i) {
    if(BLOCK[i]) falcON_DEL_O(BLOCK[i]);
    BLOCK[i] = 0;
  }
  NBLK = 0u;
  for(bodytype t; t; ++t) {
    NALL [t] = 0u;
    NBOD [t] = 0u;
    TYPES[t] = 0;
  }
  NTOT  = 0;
  FIRST = 0;
}
// destruction: delete all data
bodies::~bodies() falcON_THROWING
{
  DebugInfo(6,"bodies::~bodies(): destructing bodies");
  BITS = fieldset::empty;
  if(C_FORTRAN)
    for(fieldbit f; f; ++f)
      const_cast<block*>(FIRST)->set_data_void(f,0);
  del_data();
}
// set block::FIRST, LOCALFIRST, NALL, NBOD & NTOT
void bodies::set_firsts()
{
  for(bodytype t; t; ++t) {
    NALL[t] = 0u;
    NBOD[t] = 0u;
  }
  NTOT = 0u;
  for(block*P=FIRST; P; P=P->next()) {
    P->set_first(NTOT);
    NALL[P->type()] += P->N_alloc ();
    NBOD[P->type()] += P->N_bodies();
    NTOT            += P->N_bodies();
  }
}
// set up blocks to hold N[t] bodies of type t
void bodies::set_data(const unsigned N[bodytype::NUM]) falcON_THROWING
{
  DebugInfo(5,"bodies::set_data(): N=[%d,%d,%d], BITS=%s\n",
	    N[0],N[1],N[2],word(BITS));
  del_data();
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
	DebugInfo(10,"allocated %s @ %p\n",
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
  DebugInfo(6,"bodies::set_data(): done\n");
}
// construction 0: construction with N=0, but data fields
bodies::bodies(fieldset bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 ),
  FORCES    ( 0 )
{
  unsigned n[bodytype::NUM]={0u};
  DebugInfo(2,"bodies::bodies(): constructing bodies @%p: n=%u,%u,%u, bits=%s",
	    this,n[0],n[1],n[2],word(BITS));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
  DebugInfo(2,"bodies::bodies(): constructed\n");
}
// construction 1, new version
bodies::bodies(const unsigned n[bodytype::NUM],
	       fieldset       bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 ),
  FORCES    ( 0 )
{
  DebugInfo(2,"bodies::bodies(): constructing bodies @%p: n=%u,%u,%u, bits=%s",
	    this,n[0],n[1],n[2],word(bits));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
// resets N, data; same as destruction followed by constructor 1            
void bodies::reset(const unsigned n[bodytype::NUM],
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
// construction 2:
// just make a copy of existing bodies:
// - only copy data specified by 2nd argument
// - only copy bodies whose flags matches 3rd argument
// - only copy bodies whose type is contained in 4th argument
bodies::bodies(bodies const&Other,
	       fieldset     copydata,
	       flags        copyflag,
	       bodytypes    copytypes) falcON_THROWING :
  BITS      ( copydata & Other.BITS ),
  C_FORTRAN ( 0 ),
  FORCES    ( 0 )
{
  if(copyflag && !Other.have_flag() ) 
    falcON_THROW("in bodies::bodies(): "
		 "copyflag !=0, but other bodies not supporting flag");
  unsigned n[bodytype::NUM]={0};
  for(bodytype t; t; ++t) if(copytypes.contain(t)) {
    if(copyflag) {
      LoopTypedBodies(&Other,i,t)
	  if( falcON::flag(i).are_set(copyflag) ) ++(n[t]);
    } else 
      n[t] = Other.NBOD[t];
  }
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
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
// construction for C & FORTRAN support                                     
bodies::bodies(char, const unsigned n[bodytype::NUM]) falcON_THROWING
: BITS      ( fieldset::empty ),
  C_FORTRAN ( 1 ),
  FORCES    ( 0 )
{
  DebugInfo(3,"bodies::bodies(): constructing bodies for C & FORTRAN: n=%u,%u",
	    n[0],n[1]);
  for(bodytype t; t; ++t)
    if(n[t] > index::max_bodies)
      falcON_THROW("too many bodies\n");
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
//
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
//
void bodies::swap_bytes(fieldbit f) falcON_THROWING
{
  if(!BITS.contain(f))
    for(const block*p=FIRST; p; p=p->next())
      const_cast<block*>(p)->swap_bytes(f);
}
//
void bodies::add_field(fieldbit f) falcON_THROWING
{
  if(!BITS.contain(f)) {
    for(const block*p=FIRST; p; p=p->next())
      const_cast<block*>(p)->add_field(f);
    BITS |= fieldset(f);
    if(f == fieldbit::k) reset_keys();
  }
}
//
void bodies::add_fields(fieldset b) falcON_THROWING
{
  if(!BITS.contain(b)) {
    for(const block *p=FIRST; p; p=p->next())
      const_cast<block*>(p)->add_fields(b);
    if(!BITS.contain(fieldbit::k) && b.contain(fieldbit::k)) reset_keys();
    BITS |= b;
  }
}
//
void bodies::del_field(fieldbit f) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->del_field(f);
  BITS &= ~(fieldset(f));
}
//
void bodies::del_fields(fieldset b) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->del_fields(b);
  BITS &= ~b;
}
//
void bodies::remove(bodytype t)
{
  for(block*P=TYPES[t]; P && P->TYPE == t; P=P->NEXT)
    P->remove(NDEL[t]);
  set_firsts();
  DebugInfo(5,"bodies::remove(%s): removed %d bodies\n", t.name(), NDEL[t]);
}
//
void bodies::remove() {
  for(block*P=FIRST; P; P=P->NEXT)
    P->remove(NDEL[P->TYPE]);
  set_firsts();
  DebugInfo(5,"bodies::remove(): removed %d,%d,%d bodies\n",
	    NDEL[0],NDEL[1],NDEL[2]);
}
//
void bodies::apply_filter(BodyFilter const&F, bool zm, bool warn)
  falcON_THROWING
{
  if(F.is_empty()) return;
  // remember original data fields
  fieldset ORIG = BITS;
  // ensure we have flags (needed for removal)
  if(!have(fieldbit::f)) {
    add_field(fieldbit::f);
    reset_flags();
  }
  // ensure we have all data needed for filter
  if(!BITS.contain(F.need())) {
    fieldset miss = BITS.missing(F.need());
    if(zm) {
      if(warn) 
	falcON_Warning("snapshot::apply_filter(): data '%s' required for filter"
		       " are not supported; will assume zero values instead\n",
		       word(miss));
      add_fields(miss);
      reset_data(miss);
    } else {
      if(!ORIG.contain(fieldbit::f)) del_field(fieldbit::f);
      falcON_THROW("snapshot::apply_filter(): data '%s' required for filter"
		   " are not supported\n", word(miss));
    }
  }
  // apply filter: remove bodies not matching the filter
  LoopAllBodies(this,b)
    if(!F(b)) b.flag_for_removal();
  remove();
  // delete any fields added
  del_fields(BITS-ORIG);
}
//
void bodies::merge(bodies&Other) falcON_THROWING
{
  if(NBLK + Other.NBLK > index::max_blocks)
    falcON_THROW("bodies::merge(): too many blocks\n");
  // add blocks
  for(block*P=Other.FIRST; P; P=P->next())
    add_block(P);
  // reset Other (so that its destructor will not harm us)
  Other.FIRST = 0;
  for(bodytype t; t; ++t) {
    Other.TYPES[t] = 0;
    Other.NALL [t] = 0;
    Other.NBOD [t] = 0;
    Other.NNEW [t] = 0;
    Other.NDEL [t] = 0;
  }
  Other.BITS = fieldset::empty;
  Other.NTOT = 0;
  Other.NBLK = 0;
  for(int i=0; i!=index::max_blocks; ++i)
    Other.BLOCK[i] = 0;
}
//
namespace {
  bodies *CopyFrom, *CopyTo;
  Array<bodies::index> IndexTable;
  template<int BIT> struct CopyInOrder {
    static void act(bodytype t)
    {
      unsigned i=0;
      LoopTypedBodies(CopyTo,b,t)
	b.datum<BIT>() = CopyFrom->const_datum<BIT>(IndexTable[i++]);
    }
  };
}
//
void bodies::apply_sort(BodyFunc<real> const&F, double time, fieldset copy,
			bool zm, bool warn) falcON_THROWING
{
  if(F.is_empty()) return;
  copy &= BITS;
  // ensure we have all data needed for function
  if(!BITS.contain(F.need())) {
    fieldset miss = BITS.missing(F.need());
    if(zm) {
      if(warn) 
	falcON_Warning("snapshot::apply_sort(): data '%s' required for function"
		       " are not supported; will assume zero values instead\n",
		       word(miss));
      add_fields(miss);
      reset_data(miss);
    } else
      falcON_THROW("snapshot::apply_sort(): data '%s' required for function"
		   " are not supported\n", word(miss));
  }
  // prepare other set of bodies to copy into
  bodies Other(NBOD,copy);
  CopyFrom = this;
  CopyTo   =&Other;
  // loop bodytypes
  for(bodytype t; t; ++t)
    if(NBOD[t]) {
      // obtain index table in ascending order of F(b)
      {
	Array<real>  Qf(NBOD[t]);
	Array<index> If(NBOD[t]);
	unsigned i=0;
	LoopTypedBodies(this,b,t) {
	  Qf[i] = F(b,time);
	  If[i] = index(b);
	  ++i;
	}
	Array<int> In(NBOD[t]);
	HeapIndex(Qf,In);
	IndexTable.reset(NBOD[t]);
	for(i=0; i!=NBOD[t]; ++i)
	  IndexTable[i] = If[In[i]];
      }
      // loop fields and copy data
      LoopFields<CopyInOrder>::const_some(t, copy & t.allows());
    }
  // replace us with Other
  del_data();
  merge(Other);
}
//
bodies::block* bodies::ensure_contiguous(unsigned N, bodytype t, unsigned Na)
{
  // do we have a contiguous set of at least N free bodies?
  block*P=TYPES[t];
  while(P && P->N_free()==0) P=P->next_of_same_type();
  unsigned Nf=0;
  for(block*B=P; B; B=B->next_of_same_type()) {
    if     (P==B)       Nf = B->N_free ();       // first part
    else if(B->NBOD==0) Nf+= B->N_alloc();       // further part
    else {                                       // oops: non-contiguous
      while(B && B->N_free()==0) B=B->next_of_same_type();
      Nf= B? B->N_free() : 0;
      P = B;
    }
    if(Nf >= N) break;
  }
  // YES: return block with first part
  if(Nf >= N) {
    DebugInfo(5,"bodies::ensure_contiguous(): found contiguous chunk\n");
    return P;
  }
  // NO:  have to make it
  DebugInfo(5,"bodies::ensure_contiguous(): making new block ...\n");
  return new_block(t,max(Na,N),0,BITS);
}
//
bodies::iterator bodies::new_bodies(unsigned N, bodytype t, unsigned Na)
  falcON_THROWING
{
  // ensure we have enough bodies to activate
  block*P=ensure_contiguous(N,t,Na);
  if(0==P || 0==P->N_free())
    falcON_THROW("bodies::new_bodies(): error in ensure_contiguous()\n");
  unsigned n=N;
  iterator I0(P,P->NBOD);
  // activate bodies
  for(block*B=P; n && B; B=B->next_of_same_type()) {
    unsigned s = min(B->N_free(),n);
    B->NBOD += s;
    n       -= s;
  }
  if(n) falcON_THROW("bodies::new_bodies(): cannot find enough free bodies\n");
  set_firsts();
  // flag as new
  if(have(fieldbit::f)) {
    iterator IN(I0,N);
    for(iterator I(I0); I!=IN; ++I)
      I.flag().add(flags::newbody);
  }
  return I0;
}
//
bodies::iterator bodies::new_body(bodytype t, unsigned Na) falcON_THROWING
{
  // ensure we have a body to activate
  block*P=ensure_contiguous(1,t,Na);
  if(0==P || 0==P->N_free())
    falcON_THROW("bodies::new_body(): error in ensure_contiguous()\n");
  // activate body
  iterator I0(P,P->NBOD++);
  set_firsts();
  // flag as new
  if(have(fieldbit::f)) I0.flag().add(flags::newbody);
  return I0;
}
//
void bodies::joinup(bodytype t) falcON_THROWING
{
  block*To=TYPES[t];
  bool copy=false;
  for(;;) {
    // find block of type t with free space
    while(To && To->N_free()==0) To=To->next_of_same_type();
    if(0 == To) break;
    // find later block of type t with data
    block*Fr=To->next_of_same_type();
    while(Fr && Fr->NBOD==0) Fr=Fr->next_of_same_type();
    if(0 == Fr) break;
    // copy data from Fr to To
    int nc=min(To->N_free(), Fr->NBOD);
    To->copy_bodies(Fr,Fr->NBOD-nc,To->NBOD,nc);
    To->NBOD += nc;
    Fr->NBOD -= nc;
    copy = true;
  }
  if(copy) set_firsts();
}
//
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
//                                                                              
// data I/O                                                                     
//                                                                              
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
    DebugInfo(2,"bodies::read_snapshot(): phases read");
    read |= get & fieldset::phases;
    BITS |= get & fieldset::phases;
  }
  // now read data field by field
  for(fieldbit f; f; ++f) if(get.contain(f)) {
    DebugInfo(6,"bodies::read_snapshot(): f=%c: %s\n",letter(f),
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
      DebugInfo(2,"bodies::read_snapshot(): %u %c read",
		inpt.N_read(), letter(f));
      BITS |= fieldset(f);
      read |= fieldset(f);
    }
  }
  DebugInfo(1,"bodies::read_snapshot(): read=%s\n",word(read));
  if(read & fieldset::source) mark_srce_data_changed();
  if(read & fieldset::sphmax) mark_sph_data_changed();
  if(warn && want != read)
    falcON_Warning("bodies::read_snapshot: couldn't read %s",
		   word(read.missing(want)));
  return read;
}
//
void bodies::write_snapshot(snap_out const&snap,
			    fieldset       put,
			    iterator const&start,
			    unsigned       Nwrite) const falcON_THROWING
{
  if(this != start.my_bodies())
    falcON_THROW("bodies::write_snapshot(): start body is not ours");
  if(Nwrite == 0 || Nwrite > snap.Ntot()) Nwrite = snap.Ntot();
  if(start.my_index() + Nwrite > N_bodies())
    falcON_THROW("bodies::write_snapshot(): not enough data to write: "
		 "start=%d, Nwrite=%d, Nbodies=%d\n",
		 start.my_index(), Nwrite, N_bodies());
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
//
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
//
void bodies::read_simple_ascii(std::istream  &in,
			       const fieldbit*item,
			       unsigned       Ni,
			       const unsigned N[bodytype::NUM])
{
  // 1. create table of readers                                                 
  fieldset get;
  p_reader read [bodytype::NUM][100] = {{0}};
  p_writer write[bodytype::NUM][100] = {{0}};
  if(Ni > 100) {
    Ni = 100;
    falcON_Warning(" can only read the first 100 data entries\n");
  }
  for(unsigned i=0; i!=Ni; ++i) {
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
  // 2. reset N & add fields                                                    
  reset(N, BITS|get);
  // 3. loop bodies & read data whereby ignoring lines starting with '#'        
  for(bodytype t; t; ++t)
    if(N[t]) {
      if(debug(4))
	DebugInfo("bodies::read_simple_ascii(): now reading %d %s bodies...\n",
		  N[t],t.name());
      LoopTypedBodies(this,Bi,t) {
	while( in && eat_line(in,'#') ) {}
	if(!in) falcON_THROW("bodies::read_simple_ascii(): "
			     "end of input before data have been read");
	for(unsigned i=0; i!=Ni; ++i)
	  if(read[t][i]) {
	    read [t][i](in,Bi);
	    if(write[t][i]) write[t][i](std::cerr,Bi);
	  }
	if(debug(20)) std::cerr<<'\n';
	SwallowRestofLine(in);
      }
    }
}
// sorted index table
void bodies::sorted(Array<index>&table, 
		    real       (*func)(iterator const&)) const falcON_THROWING
{
  const int n = N_subset();
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  if(have(fieldbit::f)) {
    int i = 0;
    LoopSubsetBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = index(b);
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
//
void bodies::sorted(Array<index>&table, double time,
		    BodyFunc<real> const&func) const falcON_THROWING
{
  const int n = N_subset();
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  if(have(fieldbit::f)) {
    int i = 0;
    LoopSubsetBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b,time);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b,time);
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
//
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
      I[i] = index(b);
      Q[i] = func(b);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = index(b);
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
//
namespace {
  struct Nbour { real Q; bodies::index I; };
  inline bool operator<(Nbour const&a, Nbour const&b) {return a.Q<b.Q;}
//   inline bool operator>(Nbour const&a, Nbour const&b) {return a.Q>b.Q;}
//   inline bool operator<(real q, Nbour const&b) {return q<b.Q;}
//   inline bool operator>(real q, Nbour const&b) {return q>b.Q;}
//   inline bool operator<(Nbour const&a, real q) {return a.Q<q;}
//   inline bool operator>(Nbour const&a, real q) {return a.Q>q;}
  real Huge = 1.e30;
}
falcON_TRAITS(Nbour,"<anonymous>::Nbour");
//
unsigned bodies::findNeighbours(const body&B, unsigned K, Array<index>&I) const
  falcON_THROWING
{ 
  if(!have_pos())
    falcON_THROW("bodies::findNeighbours(): have no positions\n");
  Nbour*List = falcON_NEW(Nbour,K);
  for(unsigned k=0; k!=K; ++k)
    List[k].Q = Huge;
  unsigned Niac=0;
  LoopSubsetBodies(this,b) {
    real q = dist_sq(falcON::pos(B),falcON::pos(b));
    if(List->Q > q) {
      List->Q = q;
      List->I = bodies::index(b);
      MaxHeap::after_top_replace(List,K);
      ++Niac;
    }
  }
  MaxHeap::sort(List,K);
  I.reset(K);
  if(Niac < K) K = Niac;
  for(unsigned k=0; k!=K; ++k)
    I[k] = List[k].I;
  falcON_DEL_A(List);
  return K;
}
//                                                                              
// class falcON::snapshot                                                       
//                                                                              
namespace {
  class PointerBank {
    struct PterWithKey {
      friend class falcON::traits<PterWithKey>;
      const void  *pter;
      char        *key,*name;
      size_t       size;
      PterWithKey *next;
      PterWithKey(const void* p, const char*k, size_t s, const char*n,
		  PterWithKey*x)
	: pter(p),
	  key (falcON_NEW(char, strlen(k)+strlen(n)+2)),
	  name(key + strlen(k) + 1),
	  size(s),
	  next(x) 
      {
	strcpy(key ,k);
	strcpy(name,n);
      }
      ~PterWithKey() {
	falcON_DEL_A(key);
      }
    } *HEAD;
    //
  public:
    /// default constructor
    PointerBank()
      : HEAD(0) {}
    /// copy constructor
    PointerBank(PointerBank const&PB) : HEAD(0)
    {
      for(PterWithKey*P=PB.HEAD; P; P=P->next)
	HEAD = new PterWithKey(P->pter, P->key, P->size, P->name, HEAD);
    }
    /// destructor
    ~PointerBank() {
      PterWithKey*P=HEAD,*N;
      while(P) {
	N=P->next;
	falcON_DEL_O(P);
	P=N;
      }
    }
    /// add a pointer: key must not yet be known in bank
    void add(const void*p, const char* k, size_t s, const char* n)
    {
      for(PterWithKey*P=HEAD; P; P=P->next)
	if(0==strcmp(P->key, k))
	  falcON_THROW("snapshot::add_pointer(): "
		       "key '%s' is already in bank\n",k);
      HEAD = new PterWithKey(p,k,s,n,HEAD);
    }
    /// set a pointer: add if new key, else replace (type & size must match)
    void set(const void*p, const char* k, size_t s, const char* n)
    {
      for(PterWithKey*P=HEAD; P; P=P->next)
	if(0==strcmp(P->key, k)) {
	  if(strcmp(P->name, n))
	    falcON_THROW("snapshot::set_pointer(): "
			 "name mismatch ('%s' : '%s')",P->name,n);
	  if(P->size != s)
	    falcON_THROW("snapshot::set_pointer(): "
			 "size mismatch (%lu : %lu)",lu(P->size),lu(s));
	  P->pter = p;
	  return;
	}
      HEAD = new PterWithKey(p,k,s,n,HEAD);
    }
    /// delete an entry from the bank
    void del(const char* k, bool warn = 0)
    {
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
    /// return a pointer referred to by a given key
    const void*get(const char*k, size_t s, const char*n, const char*func) const
      falcON_THROWING
    {
      for(PterWithKey*P=HEAD; P; P=P->next)
	if(0==strcmp(P->key, k)) {
	  if(s != P->size)
	    falcON_THROW("snapshot::%s(): "
			 "size (%lu) does not match value in bank (%lu)\n",
			 func,lu(s),lu(P->size));
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
//
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
//
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
//
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
//
void snapshot::del_pointer(const char*k) const
{
  if(debug(4))
    DebugInfo("snapshot::del_pointer() under \"%s\"\n",k);
  if(PBNK) static_cast<PointerBank*>(PBNK)->del(k);
}
//
snapshot::snapshot(snapshot const&S,
		   fieldset       Bd,
		   flags          F,
		   bodytypes      T) falcON_THROWING
: bodies ( S,Bd,F,T ),
  TIME   ( S.TIME ),
  PBNK   ( S.PBNK? new PointerBank(*(static_cast<PointerBank*>(S.PBNK))) : 0 ),
  PARA   ( 0 )
{}
//
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
//
snapshot::~snapshot()
{
  if(PBNK) { falcON_DEL_O(static_cast<PointerBank*>(PBNK)); PBNK = 0; }
}
//
void snapshot::apply_filter(BodyFilter&F, bool zm, bool warn) falcON_THROWING
{
  F.set_time(TIME);
  bodies::apply_filter(F,zm,warn);
}
//
#ifdef falcON_NEMO
bool snapshot::read_nemo(nemo_in const&i, fieldset&r, fieldset g,
			 const char*range, bool w) falcON_THROWING
{
  if(!i.has_snapshot())
    falcON_THROW("snapshot::read_nemo(): no snapshot to read");
  snap_in s(i);

  if(s.has_time()) {
    if(range && !time_in_range(s.time(),range)) {
      r = fieldset::empty;
      return false;
    }
    TIME = s.time();
  } else
    TIME = 0.;
  bool need_reset = false;
  for(bodytype t; t; ++t)
    if(s.Nbod(t) != N_bodies(t)) need_reset = true;
  if(need_reset) reset(s.Nbod(), fieldset::empty);
  r = read_snapshot(s,g,begin_all_bodies(),N_bodies(),w);
  return true;
}
//
fieldset snapshot::read_part(snap_in  const&s, fieldset g, iterator const&b,
			     bool w, unsigned n) falcON_THROWING
{
  TIME = s.has_time()? s.time() : 0.0;
  return read_snapshot(s,g,b,n,w);
}
//
void snapshot::write_nemo(nemo_out const&o,        // I: nemo output            
			  fieldset       w,        // I: what to write          
			  iterator const&b,        // I: starting here          
			  unsigned       n) const  //[I: #, default: all]       
  falcON_THROWING
{
  unsigned i = falcON::bodyindex(b);
  if(this != b.my_bodies())
    falcON_THROW("snapshot::write_nemo() start body is not ours\n");
  if(n == 0) n = N_bodies()-i;
  else if(i + n > N_bodies()) {
    falcON_Warning("snapshot::write_nemo() cannot write %u bodies, "
		   "will only write %u\n",n,N_bodies()-i);
    n = N_bodies()-i;
  }
  unsigned nb[bodytype::NUM]={0}, nt=n, nc(0u);
  for(bodytype t; t; ++t)
    if(i < (nc+=N_bodies(t))) {
      nb[t] = min(nc-i, nt);
      i  += nb[t];
      nt -= nb[t];
    }
  snap_out s(o,nb,TIME);
  write_snapshot(s,w,b,n);
}
//
void snapshot::write_nemo(nemo_out const&o,        // I: nemo output            
			  fieldset       w) const  // I: what to write          
  falcON_THROWING
{
  snap_out s(o,N_bodies_per_type(),TIME);
  write_snapshot(s,w,begin_all_bodies(),N_bodies());
}
#endif // falcON_NEMO
//
#ifdef falcON_REAL_IS_FLOAT

namespace {
  //                                                                            
  /// structure modelled after gadget/allvars.h                                 
  //                                                                            
  struct GadgetHeader {
    unsigned int npart[6];          ///< # particles per type in this file
    double       masstab[6];        ///< if non-zero: mass of particle of type
    double       time;              ///< simulation time of snapshot
    double       redshift;          ///< redshift of snapshot
    int          flag_sfr;
    int          flag_feedback;
    unsigned int npartTotal[6];     ///< # particles per type in whole snapshot
    int          flag_cooling;
    int          num_files;         ///< # file for this snapshot
    double       BoxSize;
    double       Omega0;
    double       OmegaLambda;
    double       HubbleParam;
    int          flag_stellarage;
    int          flag_metals;
    unsigned int npartTotalHighWord[6];
    int          flag_entropy_instead_u;
    char         fill[60];          ///< to get sizeof(GadgetHeader)=256
    /// default constructor: set all data to 0
    GadgetHeader();
    /// try to read a GadgetHeader from an input file
    ///
    /// If the size of the Fortran record == 256 == sizeof(GadgetHeader), we
    /// read the header and return true.\n
    /// If the size of the Fortran record == byte_swapped(256), then we assume
    /// the file is of different endianess. We read the header, byte-swap it
    /// and return true.\n
    /// Otherwise, the data are not consistent with a GadgetHeader, so we return
    /// false.
    /// \return have read successfully
    /// \param  in   input stream to read from
    /// \param  rec  size of Fortran record header (must be 4 or 8)
    /// \param  swap (output) need byte-swap?
    bool Read(input& in, unsigned rec, bool& swap)
      throw(falcON::exception);
    /// check whether two GadgetHeaders could possibly come from different data
    /// files for the same snapshot
    bool mismatch(GadgetHeader const&H) const;
    /// on some ICs, npartTotal[] = 0. Here we remedy for this error
    void check_simple_npart_error();
    /// dump all the header data
    void dump(std::ostream&out) const;
  };
  //
  GadgetHeader::GadgetHeader() :
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
  //
  bool GadgetHeader::Read(input& in, unsigned rec, bool& swap)
    throw(falcON::exception)
  {
    swap = 0;
    // read record header and determine swapping necessity
    if(rec == 4) {
      uint32 S;
      in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint32));
      if(S != sizeof(GadgetHeader)) {
	swap_bytes(S);
	if(S == sizeof(GadgetHeader)) swap = 1;
	else return false;
      }
    } else if(rec == 8) {
      uint64 S;
      in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint64));
      if(S != sizeof(GadgetHeader)) {
	swap_bytes(S);
	if(S == sizeof(GadgetHeader)) swap = 1;
	else return false;
      }
    } else
      throw falcON::exception("Fortran header size must be 4 or 8\n");
    // read full GadgetHeader
    in.read(static_cast<char*>
	    (static_cast<void*>(this)), sizeof(GadgetHeader));
    // if required swap bytes
    if(swap) {
      swap_bytes(npart,6);
      swap_bytes(masstab,6);
      swap_bytes(time);
      swap_bytes(redshift);
      swap_bytes(flag_sfr);
      swap_bytes(flag_feedback);
      swap_bytes(npartTotal,6);
      swap_bytes(flag_cooling);
      swap_bytes(num_files);
      swap_bytes(BoxSize);
      swap_bytes(Omega0);
      swap_bytes(OmegaLambda);
      swap_bytes(HubbleParam);
      swap_bytes(flag_stellarage);
      swap_bytes(flag_metals);
      swap_bytes(npartTotalHighWord,6);
      swap_bytes(flag_entropy_instead_u);
    }
    // read record trailer and check for consistency
    if(rec == 4) {
      uint32 S;
      in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint32));
      if(swap) swap_bytes(S);
      if(S != sizeof(GadgetHeader)) {
	falcON_Warning("GadgetHeader::Read(): record size mismatch\n");
	return false;
      }
    } else if(rec == 8) {
      uint64 S;
      in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint64));
      if(swap) swap_bytes(S);
      if(S != sizeof(GadgetHeader)) {
	falcON_Warning("GadgetHeader::Read(): record size mismatch\n");
	return false;
      }
    }
    return true;
  }
  //
  bool GadgetHeader::mismatch(GadgetHeader const&H) const {
    bool okay = true;
#define CHECK_I(FIELD,FIELDNAME)					\
    if(FIELD != H.FIELD) {						\
      okay = false;							\
      falcON_Warning("GadgetHeader \"%s\" mismatch (%u vs %u)\n",	\
		     FIELDNAME, FIELD, H.FIELD);			\
    }
#define CHECK_D(FIELD,FIELDNAME)					\
    if(FIELD != H.FIELD) {						\
      okay = false;							\
      falcON_Warning("GadgetHeader\"%s\" mismatch (%f vs %f)\n",	\
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
  //
  void GadgetHeader::check_simple_npart_error() {
    for(int k=0; k!=6; ++k)
      if(npart[k] > npartTotal[k]) {
	falcON_Warning("GadgetHeader: npart[%u]=%u > npartTotal[%u]=%u: "
		       "we will try to fix by setting npartTotal[%u]=%u\n",
		       k,npart[k],k,npartTotal[k],k,npart[k]);
	npartTotal[k] = npart[k];
      }
  }
  //
  void GadgetHeader::dump(std::ostream&out) const {
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
} // namespace {
falcON_TRAITS(::GadgetHeader,"GadgetHeader");
//
#define READ(BIT)							\
  if((!is_sph(BIT) && nd) || ns) {					\
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
		     "expected %lu bytes, found %u\n",			\
		     nr,field_traits<BIT>::word(),			\
		     lu(nr*field_traits<BIT>::size),F.size());		\
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
//
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
  unsigned NB[bodytype::NUM] = {0}, NP[6] ={0};
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
		       "expected %lu bytes, found %u\n",
		       nm,lu(nm*sizeof(real)), F.size());
	body sph(SPH), std(STD);
	add_field(fieldbit::m);
	if(header->npart[0]) {
	  if(header->masstab[0])
	    for(unsigned b=0; b!=header->npart[0]; ++b,++sph)
	      sph.mass() = header->masstab[0];
	  else
	    sph.read_Fortran(F, fieldbit::m, header->npart[0], swap);
	}
	for(int k=1; k!=6; ++k) if(header->npart[k]) {
	  if(header->masstab[k])
	    for(unsigned b=0; b!=header->npart[k]; ++b,++std)
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
	for(unsigned b=0; b!=header->npart[0]; ++b,++sph)
	  sph.mass() = header->masstab[0];
      for(int k=1; k!=6; ++k) if(header->npart[k]) {
	for(unsigned b=0; b!=header->npart[k]; ++b,++std)
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
//
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
  header->time = time;
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
//
#endif // falcON_REAL_IS_FLOAT
////////////////////////////////////////////////////////////////////////////////
