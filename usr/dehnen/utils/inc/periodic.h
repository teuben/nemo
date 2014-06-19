// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/periodic.h
///
/// \brief  periodic boundary conditions
///
/// \author Walter Dehnen
///
/// \date   2010-2012
///
/// \version Jun-2012 WD  created based on falcON.2's period.h
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010-2012 Walter Dehnen
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
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_periodic_h
#define WDutils_included_periodic_h

#ifndef WDutils_included_geometry_h
#  include "geometry.h"
#endif
#ifndef WDutils_included_inline_h
#  include "inline.h"
#endif

namespace WDutils {
#if __cplusplus < 201103L
#  define noexcept
#endif
  ///
  /// type used to hold offset bitset for positional offsets
  ///
  typedef uint8_t positional_offset_bits;
  ///
  /// trivial parts of PeriodicBox
  ///
  template<int _tD, typename _tX> class PeriodicBoxBase
  {
  public:
    typedef _tX real;
#if __cplusplus < 201103L
    typedef tupel<_tD,bool> bvec;
    typedef tupel<_tD,real> xvec;
#else
    typedef vector<_tD,bool> bvec;
    typedef vector<_tD,real> xvec;
#endif
  protected:
    /// is x within the periodic range?
    static bool inside(real x, real hp) noexcept
    { return -hp<=x && x<=hp; }
    /// shift x into the periodic range
    static void shift(real&x, real hp, real fp) noexcept
    {
      while(x<-hp) x+= fp;
      while(x> hp) x-= fp;
    }
    /// shift x into the periodic range, but not by more than a full period
    static void shift_short(real&x, real hp, real fp) noexcept
    {
      if     (x<-hp) x+=fp;
      else if(x> hp) x-=fp;
    }
    /// \name data
    //@{
    xvec HalfP;                 ///< half period per dimension
    xvec FullP;                 ///< full period per dimension
    xvec InvFP;                 ///< inverse of full periods
    real MinHP;                 ///< minimum non-zero HalfP
    bvec PiDim;                 ///< HalfP[d]>0?
    bool PiAny;                 ///< periodic in any d?
    bool PiAll;                 ///< periodic in all d?
    //@}
    //  deprecated
    typedef positional_offset_bits offsetbits;
    ///
    /// default ctor: open boundary in all dimensions
    ///
    PeriodicBoxBase() noexcept
      : HalfP(real(0))
      , FullP(real(0))
      , InvFP(real(0))
      , MinHP(0)
      , PiDim(false)
      , PiAny(false)
      , PiAll(false) {}
    ///
    /// ctor from given half range
    /// \param[in]  hp   half the period in each dimension
    /// \note if hp[d]=0, no periodic boundary is imposed in dimension d
    ///
    explicit PeriodicBoxBase(xvec const&hp) noexcept
    { SetPeriod(hp); }
  public:
    ///
    /// set data (like re-construction)
    /// \param[in]  hp   half the period in each dimension
    /// \note if hp[d]=0, no periodic boundary is imposed in dimension d
    ///
    void SetPeriod(xvec const&hp) noexcept
    {
      MinHP = 0;
      PiDim = 0;
      PiAny = 0;      
      PiAll = 1;      
      for(int d=0; d!=_tD; ++d) {
	HalfP[d] = abs(hp[d]);
	FullP[d] = HalfP[d]+HalfP[d];
	InvFP[d] = FullP[d]>0? real(1)/FullP[d] : real(0);
	if(HalfP[d]>0) {
	  PiDim[d] = true;
	  PiAny    = true;
	  if(MinHP<=0 || HalfP[d]<MinHP) MinHP=HalfP[d];
	} else
	  PiAll    = false;
      }
    }
    ///
    /// are we implementing the same periodic box as other?
    ///
    bool operator==(PeriodicBoxBase const&other) const noexcept
    { return HalfP == other.HalfP; }
    ///
    /// are we implementing another periodic box as other?
    ///
    bool operator!=(PeriodicBoxBase const&other) const noexcept
    { return HalfP != other.HalfP; }
    ///
    /// half period in dimension d
    ///
    real const&half_period(int d) const noexcept
    { return HalfP[d]; }
    ///
    /// full period in dimension d
    ///
    real const&full_period(int d) const noexcept
    { return FullP[d]; }
    ///
    /// is there any periodicity in any dimension?
    ///
    bool is_periodic() const noexcept
    { return PiAny; }
    ///
    /// is there periodicity in all dimensions?
    ///
    bool is_fully_periodic() const noexcept
    { return PiAll; }
    ///
    /// is there periodicity in given dimension?
    ///
    bool is_periodic(int d) const noexcept
    { return PiDim[d]; }
    ///
    /// minimum non-zero half-period
    ///
    real min_half_period() const noexcept
    { return MinHP; }
    ///
    /// is a given (smoothing) length okay for find_ghosts()?
    ///
    bool length_okay(real h) const noexcept
    { return !PiAny || h<MinHP; }
    ///
    /// is a given (smoothing) length okay for find_ghosts()?
    ///
    template<typename scalar>
    bool length_okay(scalar h) const noexcept
    { return !PiAny || h<MinHP; }
    ///
    /// is a given (smoothing) length-squared okay for find_ghosts()?
    ///
    bool length_sq_okay(real q) const noexcept
    { return !PiAny || q<MinHP*MinHP; }
    ///
    /// given an offset bitset, return the positional offset in given dimension
    ///
    template<int D>
    real offset(offsetbits off) const noexcept
    { return off&(1<<(D+D))? -FullP[D] : off&(2<<(D+D))? FullP[D] : real(0); }
  };// class PeriodicBoxBase<,> 
  
  ///
  /// support for periodic boundary conditions
  ///
  /// In each spatial dimension d, the periodic box extends from -HalfP[d] to
  /// +HalfP[d] unless HalfP[d]=0, in which case we assume an open boundary.
  ///
  template<int _tD, typename _tX> class PeriodicBox;

  ///
  /// specialisation for 2D
  ///
  template<typename _tX> class PeriodicBox<2,_tX>
  : public PeriodicBoxBase<2,_tX> 
  {
    typedef PeriodicBoxBase<2,_tX> base;
    using base::HalfP;
    using base::FullP;
    using base::InvFP;
    using base::MinHP;
    using base::PiDim;
    using base::PiAny;
    using base::inside;
    using base::shift;
    using base::shift_short;
    using base::length_sq_okay;
  public:
    typedef typename base::real real;
    typedef typename base::bvec bvec;
    typedef typename base::xvec xvec;
    typedef positional_offset_bits offsetbits;
    ///
    /// default ctor: open boundary in all dimensions
    ///
    PeriodicBox()
#if __cplusplus >= 201103L
    = default;
#else
    : base() {}
#endif
    ///
    /// ctor from given half range
    /// \param[in]  hp   half the period in each dimension
    /// \note if hp[d]=0, no periodic boundary is imposed in dimension d
    ///
    explicit PeriodicBox(xvec const&hp) noexcept : base(hp) {}
    ///
    /// given a velocity, give the maximum 1/step so that any particle will
    /// not move by more than a full period in any dimension
    ///
    real max_inv_step(xvec const&v) const noexcept
    {
      real p=0;
      if(PiDim[0]) update_max(p,InvFP[0]*abs(v[0]));
      if(PiDim[1]) update_max(p,InvFP[1]*abs(v[1]));
      return p;
    }
    ///
    /// is a position inside the periodic box?
    ///
    bool is_inside(xvec const&x) const noexcept
    { 
      return  !PiAny || ( ( !PiDim[0] || inside(x[0],HalfP[0]) ) &&
			  ( !PiDim[1] || inside(x[1],HalfP[1]) ) );
    }
    ///
    /// shift a position into the periodic box
    ///
    /// \param[in,out] x  position to map back periodically into box
    ///
    void make_periodic(xvec&x) const noexcept
    {
      if(!PiAny) return;
      if(PiDim[0]) shift(x[0],HalfP[0],FullP[0]);
      if(PiDim[1]) shift(x[1],HalfP[1],FullP[1]);
    }
    ///
    /// shift a position by at most one full period into the periodic box
    ///
    /// \param[in,out] x  position to map back periodically into box
    ///
    /// \note This routine fails to bring a position back into the periodic box
    ///       if it was away more than one full period.
    ///
    void make_periodic_short(xvec&x) const noexcept
    {
      if(!PiAny) return;
      if(PiDim[0]) shift_short(x[0],HalfP[0],FullP[0]);
      if(PiDim[1]) shift_short(x[1],HalfP[1],FullP[1]);
    }
    ///
    /// shortest periodic distance between two positions within periodic box
    ///
    /// \param[in] x  position inside box (assumed and not asserted)
    /// \param[in] y  position inside box (assumed and not asserted)
    /// \return x-y shifted into box
    ///
    xvec distance(xvec const&x, xvec const&y) const noexcept
    {
      xvec d=x-y;
      make_periodic_short(d);
      return d;
    }
    ///
    /// squared distance between two positions (not necessarily within box)
    ///
    /// \param[in] x  position inside box (assumed and not asserted)
    /// \param[in] y  position inside box (assumed and not asserted)
    /// \return periodic(x-y)^2
    ///
    real dist_sq(xvec const&x, xvec const&y) const noexcept
    {
      if(!PiAny) return WDutils::dist_sq(x,y);
      real q=0,d;
      d=x[0]-y[0]; if(PiDim[0]) shift(d,HalfP[0],FullP[0]); q+=d*d;
      d=x[1]-y[1]; if(PiDim[1]) shift(d,HalfP[1],FullP[1]); q+=d*d;
      return q;
    }
    ///
    /// is an offset bitset possibly correct
    ///
    bool is_okay(offsetbits off) const noexcept
    {
      return 
	( PiDim[0]? !((off& 1) && (off& 2)) : !(off& 3) ) &&
	( PiDim[1]? !((off& 4) && (off& 8)) : !(off&12) );
    }
    ///
    /// given an offset bitset, return the positional offset
    ///
    xvec offset(offsetbits off) const noexcept
    {
      return xvec(off& 1? -FullP[0] : off& 2? FullP[0] : real(0),
		  off& 4? -FullP[1] : off& 8? FullP[1] : real(0));
    }
    ///
    /// given @c offsetbits, return those corresponding to the reversed
    ///
    static offsetbits reversed(offsetbits old) noexcept
    {
       offsetbits off=0;
       if(old& 1) off |= 2; else if(old& 2) off |= 1;
       if(old& 4) off |= 8; else if(old& 8) off |= 4;
       return off;
    }
    ///
    /// render offset bitset such that offset(original) = -offset(final)
    ///
    static void reverse(offsetbits&off) noexcept
    { off = reversed(off); }
    ///
    /// find periodic ghosts positions which neighbour the box
    ///
    /// \param[in]  x    position within the box
    /// \param[in]  q    radius-squared of search sphere
    /// \param[out] off  offset bitsets for ghost positions
    /// \return          number of neighbouring ghosts
    ///
    /// \note We only check for the 3/7 (2D/3D) nearest ghosts. Further ghosts
    ///       can only be neighbouring the box if @a q > (half-period)^2 in
    ///       any dimension; see also length_okay() and length_sq_okay().
    ///
    /// \note The offset bitset has bit (2d) set if the offset is -FullP[d]
    ///       and the bit (2d+1) set if the offset if +FullP[d] (otherwise the
    ///       offset in dimension d is zero).
    ///
    unsigned find_ghosts(xvec const&x, real q, offsetbits off[7]) const
    {
      WDutilsAssert(length_sq_okay(q));
      unsigned n=0;
      if(PiAny) {
	real  qi[2] = {real(0)};
	int   in[2] = {1,1};
	if(PiDim[0]) { qi[0] = square(HalfP[0]-std::abs(x[0])); in[0]=2; }
	if(PiDim[1]) { qi[1] = square(HalfP[1]-std::abs(x[1])); in[1]=2; }
	for(int i0=0; i0!=in[0]; ++i0) {
	  real q0 = i0? qi[0] : 0;
	  if(q0>q) continue;
	  unsigned o0 = i0? x[0]>0? 1:2:0;
	  for(int i1=0; i1!=in[1]; ++i1) {
	    real q1 = q0 + (i1? qi[1] : 0);
	    if(q1>q) continue;
	    unsigned o1 = o0 | (i1? x[1]>0? 4:8:0);
	    if(o1) off[n++] = offsetbits(o1);
	  }
	}
      }
      return n;
    }
    ///
    /// distance-squared of ghost position from periodic box
    ///
    /// \param[in] x    position within the box
    /// \param[in] off  offset indicator
    /// \return    q    distance-squared of position x+offset(off) to box
    ///
    /// \note if x+offset(off) is within the box zero is returned
    ///
    /// \note @a find_ghosts() finds all ghosts with distance-squared > 0.
    ///
    real ghost_dist_sq(xvec const&x, offsetbits off) const noexcept
    {
      real q=0;
      if(PiAny) {
	if(PiDim[0] && (off& 3)) q+= square(HalfP[0]-std::abs(x[0]));
	if(PiDim[1] && (off&12)) q+= square(HalfP[1]-std::abs(x[1]));
      }
      return q;
    }
  };// class WDutils::PeriodicBox<2,X>
  ///
  /// specialisation for 2D
  ///
  template<typename _tX> class PeriodicBox<3,_tX>
  : public PeriodicBoxBase<3,_tX> 
  {
    typedef PeriodicBoxBase<3,_tX> base;
    using base::HalfP;
    using base::FullP;
    using base::InvFP;
    using base::MinHP;
    using base::PiDim;
    using base::PiAny;
    using base::inside;
    using base::shift;
    using base::shift_short;
    using base::length_sq_okay;
  public:
    typedef typename base::real real;
    typedef typename base::bvec bvec;
    typedef typename base::xvec xvec;
    typedef typename base::offsetbits offsetbits;
    ///
    /// default ctor: open boundary in all dimensions
    ///
    PeriodicBox()
#if __cplusplus >= 201103L
    = default;
#else
    : base() {}
#endif
    ///
    /// ctor from given half range
    /// \param[in]  hp   half the period in each dimension
    /// \note if hp[d]=0, no periodic boundary is imposed in dimension d
    ///
    explicit PeriodicBox(xvec const&hp) noexcept : base(hp) {}
    ///
    /// given a velocity, give the maximum 1/step so that any particle will
    /// not move by more than a full period in any dimension
    ///
    real max_inv_step(xvec const&v) const noexcept
    {
      real p=0;
      if(PiDim[0]) update_max(p,InvFP[0]*abs(v[0]));
      if(PiDim[1]) update_max(p,InvFP[1]*abs(v[1]));
      if(PiDim[2]) update_max(p,InvFP[2]*abs(v[2]));
      return p;
    }
    ///
    /// is a position inside the periodic box?
    ///
    bool is_inside(xvec const&x) const noexcept
    { 
      return  !PiAny || ( ( !PiDim[0] || inside(x[0],HalfP[0]) ) &&
			  ( !PiDim[1] || inside(x[1],HalfP[1]) ) &&
			  ( !PiDim[2] || inside(x[2],HalfP[2]) ) );
    }
    ///
    /// shift a position into the periodic box
    ///
    /// \param[in,out] x  position to map back periodically into box
    ///
    void make_periodic(xvec&x) const noexcept
    {
      if(!PiAny) return;
      if(PiDim[0]) shift(x[0],HalfP[0],FullP[0]);
      if(PiDim[1]) shift(x[1],HalfP[1],FullP[1]);
      if(PiDim[2]) shift(x[2],HalfP[2],FullP[2]);
    }
    ///
    /// shift a position by at most one full period into the periodic box
    ///
    /// \param[in,out] x  position to map back periodically into box
    ///
    /// \note This routine fails to bring a position back into the periodic box
    ///       if it was away more than one full period.
    ///
    void make_periodic_short(xvec&x) const noexcept
    {
      if(!PiAny) return;
      if(PiDim[0]) shift_short(x[0],HalfP[0],FullP[0]);
      if(PiDim[1]) shift_short(x[1],HalfP[1],FullP[1]);
      if(PiDim[2]) shift_short(x[2],HalfP[2],FullP[2]);
    }
    ///
    /// shortest periodic distance between two positions within periodic box
    ///
    /// \param[in] x  position inside box (assumed and not asserted)
    /// \param[in] y  position inside box (assumed and not asserted)
    /// \return x-y shifted into box
    ///
    xvec distance(xvec const&x, xvec const&y) const noexcept
    {
      xvec d=x-y;
      make_periodic_short(d);
      return d;
    }
    ///
    /// squared distance between two positions (not necessarily within box)
    ///
    /// \param[in] x  position inside box (assumed and not asserted)
    /// \param[in] y  position inside box (assumed and not asserted)
    /// \return periodic(x-y)^2
    ///
    real dist_sq(xvec const&x, xvec const&y) const noexcept
    {
      if(!PiAny) return WDutils::dist_sq(x,y);
      real q=0,d;
      d=x[0]-y[0]; if(PiDim[0]) shift(d,HalfP[0],FullP[0]); q+=d*d;
      d=x[1]-y[1]; if(PiDim[1]) shift(d,HalfP[1],FullP[1]); q+=d*d;
      d=x[2]-y[2]; if(PiDim[2]) shift(d,HalfP[2],FullP[2]); q+=d*d;
      return q;
    }
    ///
    /// is an offset bitset possibly correct
    ///
    bool is_okay(offsetbits off) const noexcept
    {
      return 
	( PiDim[0]? !((off& 1) && (off& 2)) : !(off& 3) ) &&
	( PiDim[1]? !((off& 4) && (off& 8)) : !(off&12) ) &&
	( PiDim[2]? !((off&16) && (off&32)) : !(off&48) );
    }
    ///
    /// given an offset bitset, return the positional offset
    ///
    xvec offset(offsetbits off) const noexcept
    {
      return xvec(off& 1? -FullP[0] : off& 2? FullP[0] : real(0),
		  off& 4? -FullP[1] : off& 8? FullP[1] : real(0),
		  off&16? -FullP[2] : off&32? FullP[2] : real(0));
    }
    ///
    /// given @c offsetbits, return those corresponding to the reversed
    ///
    static offsetbits reversed(offsetbits old) noexcept
    {
      offsetbits off=0;
      if(old& 1) off |= 2; else if(old& 2) off |= 1;
      if(old& 4) off |= 8; else if(old& 8) off |= 4;
      if(old&16) off |=32; else if(old&32) off |=16;
      return off;
    }
    ///
    /// render offset bitset such that offset(original) = -offset(final)
    ///
    static void reverse(offsetbits&off) noexcept
    { off = reversed(off); }
    ///
    /// find periodic ghosts positions which neighbour the box
    ///
    /// \param[in]  x    position within the box
    /// \param[in]  q    radius-squared of search sphere
    /// \param[out] off  offset bitsets for ghost positions
    /// \return          number of neighbouring ghosts
    ///
    /// \note We only check for the 7 nearest ghosts. Further ghosts can only be
    ///       neighbouring the box if @a q > (half-period)^2 in any dimension;
    ///       see also length_okay() and length_sq_okay().
    ///
    /// \note The offset bitset has bit (2d) set if the offset is -FullP[d]
    ///       and the bit (2d+1) set if the offset if +FullP[d] (otherwise the
    ///       offset in dimension d is zero).
    ///
    unsigned find_ghosts(xvec const&x, real q, offsetbits off[7]) const
    {
      WDutilsAssert(length_sq_okay(q));
      unsigned n=0;
      if(PiAny) {
	real  qi[3] = {real(0)};
	int   in[3] = {1,1,1};
	if(PiDim[0]) { qi[0] = square(HalfP[0]-std::abs(x[0])); in[0]=2; }
	if(PiDim[1]) { qi[1] = square(HalfP[1]-std::abs(x[1])); in[1]=2; }
	if(PiDim[2]) { qi[2] = square(HalfP[2]-std::abs(x[2])); in[2]=2; }
	for(int i0=0; i0!=in[0]; ++i0) {
	  real q0 = i0? qi[0] : 0;
	  if(q0>q) continue;
	  unsigned o0 = i0? x[0]>0? 1:2:0;
	  for(int i1=0; i1!=in[1]; ++i1) {
	    real q1 = q0 + (i1? qi[1] : 0);
	    if(q1>q) continue;
	    unsigned o1 = o0 | (i1? x[1]>0? 4:8:0);
	    for(int i2=0; i2!=in[2]; ++i2) {
	      real q2 = q1 + (i2? qi[2] : 0);
	      if(q2>q) continue;
	      unsigned o2 = o1 | (i2? x[2]>0? 16:32:0);
	      if(o2) off[n++] = offsetbits(o2);
	    }
	  }
	}
      }
      return n;
    }
    ///
    /// distance-squared of ghost position from periodic box
    ///
    /// \param[in] x    position within the box
    /// \param[in] off  offset indicator
    /// \return    q    distance-squared of position x+offset(off) to box
    ///
    /// \note if x+offset(off) is within the box zero is returned
    ///
    /// \note @a find_ghosts() finds all ghosts with distance-squared > 0.
    ///
    real ghost_dist_sq(xvec const&x, offsetbits off) const noexcept
    {
      real q=0;
      if(PiAny) {
	if(PiDim[0] && (off& 3)) q+= square(HalfP[0]-std::abs(x[0]));
	if(PiDim[1] && (off&12)) q+= square(HalfP[1]-std::abs(x[1]));
	if(PiDim[2] && (off&48)) q+= square(HalfP[2]-std::abs(x[2]));
      }
      return q;
    }
  };// class PeriodicBox<3,X>
#if __cplusplus < 201103L
#  undef noexcept
#endif
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_periodic_h
