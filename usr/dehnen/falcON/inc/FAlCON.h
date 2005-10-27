// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// FAlCON.h                                                                    |
//                                                                             |
// Copyright (C) 1999-2005 Walter Dehnen                                       |
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
//=============================================================================+
//                                                                             |
// FAlCON = Force ALgorithm with Complexity O(N)                               |
//                                                                             |
//=============================================================================+
//                                                                             |
// CONTENTS                                                                    |
// ========                                                                    |
//                                                                             |
// 1  Initialisation and related                                               |
//    1.1  Construction (initialisation)                                       |
//    1.2  Destruction (clearing up)                                           |
//    1.3  Meaning of the body flags                                           |
// 2  Generating and maintaining an oct-tree                                   |
// 3  Gravity approximation                                                    |
//    3.1  Approximating the accelerations and potentials                      |
//    3.2  Individual softening                                                |
//    3.3  A rude estimation of the mass- and number-density                   |
// 4  Full SPH support                                                         |
//    4.1  Adjusting the smoothing sizes                                       |
//    4.2  Sweep One                                                           |
//    4.3  Sweep Two                                                           |
// 5  Search for and counting of neighbours and collision partners             |
//    5.1  Collision partner search (for sticky particles) and counting        |
//    5.2  Neighbour or interaction partner search (SPH support) and counting  |
// 6  Other features                                                           |
// 7  Known bugs and problems                                                  |
//    7.1  Test-bodies are not possible                                        |
//    7.2  Bodies at identical positions                                       |
// References                                                                  |
// A  Non-public code and data                                                 |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// See the detailled comments below for how to use this code.                  |
// Some abbreviations:                                                         |
//                                                                             |
// R:       return value                                                       |
// I:       function argument is input                                         |
// O:       function argument is output                                        |
// I/O:     function argument is input & output                                |
// [I: ]    function argument is optional                                      |
//                                                                             |
// Note that in case of optional function arguments, one may only omit the     |
// last one (and the second but last if the last is already omitted, etc).     |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This code uses types defined in file types.h (via typedefs), which are:     |
//                                                                             |
// real     a single or double precision floating point number, depending on   |
//          the macro PRECISION                                                |
// areal    a single or double precision floating point number, depending on   |
//          the macro PRECISION                                                |
// vect     a tripel (for 3D) or pair (for 2D) of reals                        |
// indx     an unsigned integer type of size 16 bytes                          |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_FAlCON_h            // ensure these definitions are    |
#define falcON_included_FAlCON_h 1          // seen once only by the compiler  |
#ifndef falcON_included_types_h             //                                 |
#  include <public/types.h>                 // basic types etc.                |
#endif                                      //                                 |
#ifndef falcON_included_default_h           //                                 |
#  include <public/default.h>               // default parameters              |
#endif                                      //                                 |
#ifndef falcON_included_body_h              //                                 |
#  include <body.h>                         // body stuff                      |
#endif                                      //                                 |
//-----------------------------------------------------------------------------+
namespace falcON {                          // some forward declarations:      |
  class OctTree;                            //   oct tree structure            |
  class GravMAC;                            //   multipole acceptance criterion|
  class GravStats;                          //   statistics for gravity        |
  class GravEstimator;                      //   gravity approximation         |
  class PartnerEstimator;                   //   collision partner search      |
#ifdef falcON_SPH                           // walter's private addendum       |
  class SphEstimator;                       //   full-fledged SPH              |
  class EquationOfState;                    //   equation of state for SPH     |
#endif                                      //                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 1 INITIALISATION AND RELATED                                                |
// ============================                                                |
//                                                                             |
  class FAlCON {                            // Force ALgorithm with            |
  public:                                   // Complexity O(N)                 |
//                                                                             |
// 1.1 CONSTRUCTORS                                                            |
// ----------------                                                            |
// In order to use the code, one has first to construct an object of type      |
// 'FAlCON'. This is done in C++ by the statement                              |
//                                                                             |
//     falcON::FAlCON my_falcon( argument_list );                              |
//                                                                             |
// which invokes a constructor for class FAlCON. There are, in fact two        |
// constructors, one for the use with bodies as defined in file body.h:        |
//                                                                             |
    FAlCON (const bodies*,                     // I: bodies                    |
	    const real,                        // I: global/maximum eps        |
	    const real     =Default::theta,    //[I: tolerance parameter]      |
	    const kern_type=Default::kernel,   //[I: type of softening kernel] |
#ifdef falcON_INDI                             // non-public version only:     |
	    const bool     =false,             //[I: use individual eps?]      |
#endif                                         //                              |
	    const real     =one,               //[I: constant of gravity]      |
	    const MAC_type =theta_of_M,        //[I: type of MAC]              |
	    const int[4] =Default::direct      //[I: N_direct for gravity]     |
#ifdef falcON_SPH                              // SPH version only:            |
	   ,const int[3] =Default::SPHdirect   //[I: N_direct for SPH]         |
#endif
	    );
//                                                                             |
// Once a 'FAlCON' is constructed, one cannot change the bodies used,  but of  |
// course, one can change the data stored with these bodies. Using the flags,  |
// one can in this way obtain almost any desired behaviour, see below.         |
//                                                                             |
// The further arguments of the constructor are as follows.                    |
//                                                                             |
// EPS:   global  softening length, if SOFT==false                             |
//        ignored,                  if SOFT==true                              |
// THETA: theta or theta_min, depending on MAC.                                |
//        RECOMMENDED: THETA = 0.5 - 0.6                                       |
// KERN:  newton     : no softening                                            |
//        p0,p1,p2,p3: Plummer (P0) and related (Dehnen & Teuben, 2002)        |
// SOFT:  use individual softening length? see also approximate_gravity()      |
// MAC:   type of MAC (see above). Default is eq (13) of Dehnen (2002).        |
//                                                                             |
// All but the first two arguments may also be omitted, in which case the      |
// default values will be used.                                                |
//                                                                             |
// With the member functions                                                   |
//                                                                             |
    void reset_softening(const real,                     // I: fixed/max eps   |
			 const kern_type=Default::kernel)//[I: soft'ng kernel] |
                         const;                          //                    |
    void reset_opening(const real,                       // I: tolerance param |
		       const MAC_type=Default::mac)const;//[I: type of MAC]    |
    void reset_NewtonsG(const real) const;               // I: new G           |
//                                                                             |
// it is possible to change, after construction, the softening length and      |
// kernel, the opening criterion as well as Newton's constant of gravity.      |
//                                                                             |
//                                                                             |
// 1.2 Destructor                                                              |
// --------------                                                              |
// The class FAlCON has a destructor that will be called implicitly, whenever  |
// the scope of an object of type FAlCON terminates. The destructor can also   |
// be called explicitly, or via the delete command (if a tree was allocated    |
// by the new command).                                                        |
//                                                                             |
    ~FAlCON();                              //                                 |
//                                                                             |
//                                                                             |
// 1.3 Meaning of the body flags                                               |
// -----------------------------                                               |
//                                                                             |
// The flags are integers with the following meaning:                          |
//                                                                             |
// bit  value     meaning                                                      |
// ------------------------------------------------------------------------    |
//   1      1     this body is active, i.e. it wants update                    |
//   2      2     flagged for removal, see body.h                              |
//   3      4     this body is a SPH particle                                  |
//   4      8     this body is a sticky particle                               |
// i>4   2^(i-1)  not used                                                     |
//                                                                             |
// The default, flag=0, represents a plain body that is inactive, but still    |
// source of gravity.                                                          |
// The flag is obtained by setting the bits or, equivalently, adding the       |
// values. The flags are not used by the constructor and will unfold their     |
// effect only when the member functions, explained below, are used.           |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 2 GENERATING AND MAINTAINING AN OCT-TREE                                    |
// ========================================                                    |
//                                                                             |
// In order to establish a hierarchical tree structure from all bodies that    |
// are not flagged to be ignored, use                                          |
//                                                                             |
    void grow(int         const& =Default::Ncrit,  //[I: Ncrit]                |
	      const vect* const& =0);              //[I: pre-set root center   |
//                                                                             |
// which grows a new tree from scratch. Cells containing Ncrit or less bodies  |
// will NOT be splitted. Experiments showed that for a full force calculation  |
// Ncrit = 6-8 is about optimal and results in a few % decrease in CPU time    |
// consumption (compared to Ncrit=1) and about 20% reduction in memory.        |
//                                                                             |
// You may instead also re-use an old tree structure. In this case, nothing    |
// is done yet, but the routine approximate_gravity() below will re-compute    |
// the centers of mass and sizes of the cells (since the bodies have moved).   |
//                                                                             |
    void reuse();                             //                               |
//                                                                             |
// In this case, the logical linkage between cells and bodies is preserved     |
// and no new tree structure is established. However, the code accounts for    |
// the fact that the bodies might have moved. If the bodies have moved a lot,  |
// the sizes of the cells will be much larger than the physical size of the    |
// associated boxes, and the tree traversal will be very inefficient.          |
// However, if the bodies have moved only little, the force computation is     |
// hardly slowed down and re-using an old tree saves the costs for             |
// establishing a new one. Note that reuse() does not allow for any change in  |
// the flags indicating whether or not a body shall be ignored: changing this  |
// flag will have an effect only at the next call of grow().                   |
//                    +----------------------------+                           |
//                    | USE THIS OPTION CAREFULLY! |                           |
//                    +----------------------------+                           |
// You have been warned!                                                       |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 3 GRAVITY APPROXIMATION                                                     |
// =======================                                                     |
//                                                                             |
// 3.1 Approximating the Accelerations and Potentials                          |
// --------------------------------------------------                          |
// Once a tree structure is established, the following routine                 |
// - computes the cell's source properties (upward pass),                      |
// - performs all interactions (interaction phase), and                        |
// - evaluates the gravity for active bodies (evaluation phase).               |
// (Note that in earlier versions, the upward pass was already done at tree    |
// growth. The new version will therefore differ in the timing between tree    |
// growth and gravity approximation.)                                          |
//                                                                             |
    void approximate_gravity(bool     const& =true,   //[I: combine phases]    |
			     bool     const& =false   //[I: all or only active]|
#ifdef falcON_ADAP                                    // non-public  only:     |
			    ,real     const& =zero,   //[I: Nsoft adjust eps_i]|
			     unsigned const& =0u,     //[I: Nref  adjust eps_i]|
			     real     const& =zero,   //[I: eps_min]           |
			     real     const& =zero    //[I: max change of eps] |
#endif                                                //                       |
			     );                       //[I: N_direct]          |
//                                                                             |
// The first argument indicates whether the interaction and evaluation phases  |
// shall be interweaved, resulting in a reduced memory requirement for the     |
// Taylor-series coefficients.                                                 |
// If Newton's constant of gravity was set to 0, approximate_gravity() will    |
// issue a warning to stderr and set the acceleration and potential of active  |
// bodies (or of all, if the second argument is true) to zero.                 |
//                                                                             |
// IMPORTANT NOTE                                                              |
// Since Oct-2003, you MUST not change the bodies activity flag between tree   |
// growth (or re-growth, re-use) and a call to approximate_gravity. Whenever   |
// you change the flags, you MUST first (re-)grow the tree before you can      |
// approximate_gravity().                                                      |
//                                                                             |
//                                                                             |
// 3.2 Individual softening                                                    |
// ------------------------                                                    |
// For individual adaptive softening the routine can do more for you before    |
// the forces are actually computed:                                           |
// If Nsoft [2nd argument] is non-zero, it estimates for each active particle  |
// the local number density (using the number density of the smallest cell     |
// containing that particle and not less than Nref [3rd argument] bodies).     |
// If fac [4th arg] is zero, the bodies softening lengths are set such that,   |
// based on the estimated number density, their eps-spheres contain Nsoft      |
// bodies, but eps <= EPS [global parameter].                                  |
// If fac [4th arg] is non zero, Eps is computed in the same way, and the      |
// new softening is set to                                                     |
//        eps_new = Eps^2 / eps_old,                                           |
// with the restriction  eps_new in [eps_old/fac, eps_old*fac]                 |
// If eps_i is adjusted in this way, it will be copied back to the bodies.     |
//                                                                             |
// The last argument controls the usage of direct summation instead of the     |
// approximate Taylor-series based method. The elemnts of the array refer to   |
// the numbers N_cb^pre, N_cb^post, N_cc^post, and N_cs in Appendix B of       |
// Dehnen (2002).                                                              |
//                                                                             |
// NOTE. The default values for DIRECT and NCRIT (above) have been carefully   |
//       chosen for optimal efficiency. It is recommended not to use other     |
//       values, unless you perform some experiments on the influence on code  |
//       performance and accuracy.                                             |
//                                                                             |
// For test purposes, we also have a direct summation facility:                |
//                                                                             |
    void exact_gravity(bool const& =false          //[I: all of only active]   |
#ifdef falcON_ADAP                                 // non-public  only:        |
		      ,real     const& =zero,      //[I: Nsoft: adjust eps_i]  |
		       unsigned const& =0u,        //[I: Nref:  adjust eps_i]  |
		       real     const& =zero,      //[I: eps_min]              |
		       real     const& =zero       //[I: max change of eps]    |
#endif                                             //                          |
		       );                          //                          |
//                                                                             |
//                                                                             |
// 3.3 A crude Estimation of the Mass- and Number-Density                      |
// ------------------------------------------------------                      |
// There is also the possibility to obtain a rough estimate of the mass-,      |
// surface- or number density of bodies in the neighbourhood of every body     |
// flagged being active, via                                                   |
//                                                                             |
// estimate mass volume density                                                |
    void estimate_rho(unsigned const&,             // I: critical cell size    |
		      bool const& = 0);            //[I: all or active only?]  |
// estimate mass surface density                                               |
    void estimate_sd (unsigned const&,             // I: critical cell size    |
		      bool const& = 0);            //[I: all or active only?]  |
// estimate number volume density                                              |
    void estimate_n  (unsigned const&,             // I: critical cell size    |
		      bool const& = 0);            //[I: all or active only?]  |
//                                                                             |
// These estimates are simply the mean density within the smallest cells with  |
// more then NX (1st arg) bodies. Note, that when using test bodies (bodies    |
// with zero or tiny mass), this guess for the mass density can have terrible  |
// errors. Moreover, when the tree has not been grow()n but simply reuse()ed,  |
// these estimate will not change. Be careful using these functions.           |
// YOU HAVE BEEN WARNED.                                                       |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifdef falcON_SPH                           //                                 |
//                                                                             |
//                                                                             |
// 4 FULL SPH SUPPORT                                                          |
// ==================                                                          |
//                                                                             |
// Once a tree structure is established, you can use FAlCON to obtain          |
// estimates for SPH quantities. Currently however, this is restricted to      |
// Walter's private version of the code.                                       |
//                                                                             |
// This part of the code considers only bodies flagged as SPH particles.       |
//                                                                             |
//                                                                             |
// 4.1 Adjusting the smoothing sizes                                           |
// ---------------------------------                                           |
// The following routine iteratively adjusts the (smoothing) sizes such that   |
// the mu_i := h_i^3 * rho_i equals a given value Mu. The bodies size() fields |
// are updated.                                                                |
//                                                                             |
    void adjust_SPH_sizes(real     const&,         // I: desired Mu            |
			  real     const&,         // I: h_max for bodies      |
			  real     const&,         // I: rel error tolerated   |
			  bool     const& =0,      //[I: ALL (or only active)?]|
			  unsigned const& =10u);   //[I: max # iterations]     |
//                                                                             |
// The 1st arg gives the desired value for mu_i, the 2nd arg the maximum       |
// allowed size. The 3rd argument give the relative error in mu which can be   |
// tolerated for any body. The 4th arg indicates whether all SPH-bodies shall  |
// be considered active. The last arg gives the maximum number of iterations.  |
//                                                                             |
//                                                                             |
// 4.2 Sweep One                                                               |
// -------------                                                               |
//                                                                             |
    int SPH_sweep_one(real const&,
		      real const&,
		      real const&,
		      real const&,
		      bool const& =0);
//                                                                             |
// 4.3 Sweep Two                                                               |
// -------------                                                               |
//                                                                             |
    void SPH_sweep_two(const EquationOfState*const&,
		       real                  const&);
//                                                                             |
// 4.4 SPH Stuff                                                               |
// -------------                                                               |
//                                                                             |
    unsigned const&N_MuSmall() const;
    unsigned const&N_MuLarge() const;
    unsigned const&N_HatMax () const;
    unsigned       N_SPH_active(bool) const;
#endif                                      //                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 5 SEARCH FOR AND COUNTING OF COLLISION PARTNERS                             |
// ===============================================                             |
//                                                                             |
// After a tree has been grown, you can also use it to create interaction      |
// for SPH and sticky particles via the routine                                |
//                                                                             |
    typedef bodies::index indx_pair[2];      // type used below                |
//                                                                             |
    void make_iaction_list(indx_pair*,       // O: list of pairs of bodyindices|
			   unsigned  ,       // I: physical size of list       |
			   unsigned &,       // O: actual size of list         |
			   bool      ,       // I: use Max(h_i,h_j) OR h_i+h_j |
			   real      ,       // I: time step tau               |
			   bool      );      // I: count partners as well?     |
//                                                                             |
// In case of overflow, i.e. if the number of pairs found exceeds the size     |
// (2nd arg) of the list (1st arg), a warning is issued to stderr, but the     |
// search is not truncated, rather pairs are no longer copied into the         |
// interaction list. In this case, the value returned for the actual size      |
// (3rd arg) exceeds the maximum size (2nd arg), but reflects the actual       |
// number of interactions found.                                               |
//                                                                             |
//                                                                             |
// 5.1 Sticky-particle support                                                 |
// ---------------------------                                                 |
// In order to make a list of all pairs {i,j} of indices i,j < Nsph for which  |
//                                                                             |
//      (1) both flags indicate sticky particles,                              |
// and  (2) at least one is flagged being active,                              |
// and  (3) | (x_i+t*v_i)-(x_j+t*v_j) | < size_i+size_j  with t in [0,tau],    |
//                                                                             |
// use make_iaction_list() with the time step (5th arg) >= 0.                  |
//                                                                             |
// If the last arg is true, the number of interaction partners for each        |
// active sticky body is counted (simultaneously with making the list).        |
//                                                                             |
// If the first argument is zero (NULL pointer), but the lat one not, then we  |
// only count interaction partners, but will not compute an interactiom list.  |
//                                                                             |
//                                                                             |
// 5.2 SPH support: neighbour or interaction partner search and/or counting    |
// ------------------------------------------------------------------------    |
// In order to make a list of all pairs {i,j} of indices i,j < Nsph for which  |
//                                                                             |
//      (1) both flags indicate SPH particles,                                 |
// and  (2) at least one flagged being active,                                 |
// and  (3)     | x_i - x_j | < max(size_i,size_j)     IF Max==true            |
//          OR  | x_i - x_j | < size_i + size_j        IF Max==false           |
//                                                                             |
// use make_iaction_list() with negative time step (5th arg) and use the 4th   |
// argument to indicate which neighbourhood criterion you whish.               |
//                                                                             |
// If the last arg is true, the number of interaction partners for each active |
// SPH body is counted (simultaneously with making the list).                  |
//                                                                             |
// One may also just count interaction partners without creating an list of    |
// interactions. This is supported by the routine                              |
//                                                                             |
    void count_sph_partners(bool);
//                                                                             |
// with the argument being equivalent to the 4th arg of make_iaction_list().   |
//                                                                             |
//                                                                             |
// See file src/C/TestPair.cc for an example application for both supports.    |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 6 OTHER FEATURES                                                            |
// ================                                                            |
//                                                                             |
// There are many further routines by which the user can obtain information    |
// about the behaviour of the code.                                            |
// For the number of bodies used by FAlCON, use                                |
//                                                                             |
    unsigned No_bodies_used() const;        //                                 |
//                                                                             |
// For the number of cells in the tree, use                                    |
//                                                                             |
    unsigned No_cells_used() const;         //                                 |
//                                                                             |
// For the number of sets of Taylor series coefficients, use                   |
//                                                                             |
    unsigned const&No_coeffs_used() const;  //                                 |
//                                                                             |
// In order to dump (almost) all info of the bodies and cells to files, use    |
//                                                                             |
    void dump_nodes(const char* = 0,        //[I: dump cells to this file?]    |
		    const char* = 0) const; //[I: dump bodies to this file?]   |
//                                                                             |
// To display some statistics, use                                             |
//                                                                             |
    void stats     (std::ostream&) const;   // I: output stream                |
//                                                                             |
// To obtain the currently used MAC, use                                       |
//                                                                             |
    const MAC_type& MAC() const;            //                                 |
//                                                                             |
// To obtain description of the currently used MAC, use                        |
//                                                                             |
    const char* describe_MAC() const;       //                                 |
//                                                                             |
// To know whether or not individual softening lengths eps_i are used, use     |
//                                                                             |
#ifdef falcON_INDI                          //                                 |
    const bool& use_individual_eps() const; //                                 |
#endif                                      //                                 |
//                                                                             |
// To obtain the currently used kernel, use                                    |
//                                                                             |
    const kern_type& kernel() const;        //                                 |
//                                                                             |
// To obtain description of the currently used kernel, use                     |
//                                                                             |
    const char* describe_kernel() const;    //                                 |
//                                                                             |
// To obtain the currently used (global) softening length, use                 |
//                                                                             |
    const real& softening_length() const;   //                                 |
    const real& eps             () const;   //                                 |
//                                                                             |
// To obtain the current value of Newton's G, use                              |
//                                                                             |
    const real& NewtonsG() const;           //                                 |
//                                                                             |
// To obtain some interaction statistics, use                                  |
//                                                                             |
    unsigned BB_interactions() const;       // R: # of body-body IAs           |
    unsigned MB_interactions() const;       // R: # of many-body IAs           |
    unsigned CB_interactions() const;       // R: # of cell-body IAs           |
    unsigned CC_interactions() const;       // R: # of cell-cell IAs           |
    unsigned total_interactions() const;    // R: # of all IAs                 |
//                                                                             |
// To obtain some global tree properties, use                                  |
//                                                                             |
    const vect    &root_center() const;     // R: center of root               |
    const real    &root_radius() const;     // R: radius of root               |
    const int     &root_number() const;     // R: # bodies in root             |
    const real    &root_mass  () const;     // R: mass  of root                |
    const unsigned&root_depth () const;     // R: depth  of root               |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 7 KNOWN BUGS AND PROBLEMS                                                   |
// =========================                                                   |
//                                                                             |
// 7.1 Test-bodies are not possible                                            |
// --------------------------------                                            |
// A body that is loaded into the tree but has zero mass, will not acquire     |
// any acceleration. This is because the code computes first the force = mass  |
// times acceleration (it is symmetric and hence better suited for the         |
// computation of mutual interactions) and then divides by the mass to obtain  |
// the acceleration. The only possible work-around this problem is to set      |
// the mass of potential test bodies much smaller than the masses of source    |
// bodies. However, this must be done such that the gravity of test bodies     |
// is everywhere neglible compared to that of the source bodies.               |
// Note, however, that this work-around is wasteful: it computes the forces    |
// generated by the test bodies. (You have been warned!)                       |
//                                                                             |
//                                                                             |
// 7.2 Bodies at identical positions                                           |
// ---------------------------------                                           |
// The code cannot cope with more than Ncrit bodies at the same position       |
// (within the floating point accuracy). This situation will lead to an        |
// infinitely deep tree, i.e. the maximum allowed tree depth will be exceeded  |
// and the code aborts with an error message.                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// REFERENCES                                                                  |
// ==========                                                                  |
//                                                                             |
// Dehnen W., 2000, ApJ, 536, L39                                              |
// Dehnen W., 2001, MNRAS, 324, 273                                            |
// Dehnen W., 2002, J. Comp. Phys., 179, 27                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// A  NON-PUBLIC CODE AND CODE USED BY THE C AND FORTRAN SUPPORT               |
// =============================================================               |
//                                                                             |
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                                             |
// const protected methods                                                     |
//                                                                             |
  protected:                                      //                           |
    unsigned No_bodies        () const;           // R: # bodies               |
    OctTree* const& my_tree   () const            // R: my tree                |
    { return TREE; }                              //                           |
//-----------------------------------------------------------------------------+
//                                                                             |
// private methods                                                             |
//                                                                             |
  private:                                        //                           |
    FAlCON                    ();                 // not implemented           |
    FAlCON                    (const FAlCON&);    // not implemented           |
    FAlCON& operator=         (const FAlCON&);    // not implemented           |
//-----------------------------------------------------------------------------+
//                                                                             |
// data members                                                                |
//                                                                             |
    mutable GravStats     *STATS;                 // statistics                |
    const   bodies        *BODIES;                // sbodies to be used        |
    int                    Ncrit;                 // Ncrit                     |
    OctTree               *TREE;                  // tree to be used           |
    GravMAC               *GMAC;                  // theta(M)                  |
    GravEstimator         *GRAV;                  // gravity estimator         |
    PartnerEstimator      *PAES;                  // Partner finding & counting|
#ifdef falcON_SPH                                 //                           |
    SphEstimator          *SPHT;                  // SPH estimator             |
#endif                                            //                           |
  };                                              // END: class FAlCON         |
}                                                 // END: namespace falcON     |
//-----------------------------------------------------------------------------+
#include <public/FAlCON.cc>                       // include inline definitions|
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif                                            // falcON_included_FAlCON_h  |
