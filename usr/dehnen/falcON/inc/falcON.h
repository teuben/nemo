// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// falcON.h                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1999-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//=============================================================================+
//                                                                             |
// falcON = Force ALgorithm with Complexity O(N)                               |
//                                                                             |
//=============================================================================+
//                                                                             |
// header file for C++ users                                                   |
// (C and FORTRAN users, see files falcON_C.h and falcON.f, respectively)      |
//                                                                             |
// C++ source code implementing the code described by Dehnen (2000,2002).      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// See file readme.falcON.html for some guidelines on how to install, #include |
// and link this file and the code declared in it. You have the choice between |
// various options determining some properties of the code.                    |
//                                                                             |
//-----------------------------------------------------------------------------+
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
//    4.1  Partner Counting                                                    |
//    4.2  Adjusting the smoothing sizes                                       |
//    4.3  Sweep One                                                           |
//    4.4  Sweep Two                                                           |
// 5  Search for and counting of interaction partners                          |
//    5.1  SPH partner search                                                  |
//    5.2  sticky-particle partner search                                      |
//    5.3  SPH partner counting                                                |
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
// This code uses types defined in file auxx.h (via typedefs), which are:      |
//                                                                             |
// real     a single or double precision floating point number, depending on   |
//          the macro PRECISION                                                |
// areal    a single or double precision floating point number, depending on   |
//          the macro PRECISION                                                |
// vect     a tripel (for 3D) or pair (for 2D) of reals                        |
// indx     an unsigned short (for holding small unsigned integers)            |
// uint     an unsigned integer                                                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_falcON_h            // ensure these definitions are    |
#define falcON_included_falcON_h 1          // seen once only by the compiler  |
#ifndef falcON_included_auxx_h              //                                 |
#  include <public/auxx.h>                  // basic types etc.                |
#endif                                      //                                 |
#ifndef falcON_included_deft_h              //                                 |
#  include <public/deft.h>                  // default parameters              |
#endif                                      //                                 |
#ifndef falcON_included_nbio_h              //                                 |
#  include <public/nbio.h>                  // class nbdy::io                  |
#endif                                      //                                 |
//-----------------------------------------------------------------------------+
namespace nbdy {                            // some forward declarations:      |
  class bodies;                             //   bodies for serial code        |
  class ebodies;                            //   bodies for C or FORTRAN link  |
  class oct_tree;                           //   oct tree structure            |
  class grav_mac;                           //   multipole acceptance criterion|
  class grav_stat;                          //   statistics for gravity        |
  class grav_estimator;                     //   gravity approximation         |
  class stsp_estimator;                     //   collision partner search      |
#ifdef falcON_SPH                           // walter's private addendum       |
  class sph_estimator;                      //   full-fledged SPH              |
  class eq_of_state;                        //   equation of state for SPH     |
#endif                                      //                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 1 INITIALISATION AND RELATED                                                |
// ============================                                                |
//                                                                             |
  class falcON {                            // Force ALgorithm with            |
  public:                                   // Complexity O(N)                 |
//                                                                             |
// 1.1 CONSTRUCTORS                                                            |
// ----------------                                                            |
// In order to use the code, one has first to construct an object of type      |
// 'falcON'. This is done in C++ by the statement                              |
//                                                                             |
//     nbdy::falcON my_falcon( argument_list );                                |
//                                                                             |
// which invokes a constructor for class falcON. There are, in fact two        |
// constructors, one for the use with bodies as defined in file body.h:        |
//                                                                             |
    falcON (const bodies*,                     // I: bodies                    |
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
// and one for the use with plain arrays:                                      |
//                                                                             |
    falcON (const ebodies*,                    // I: bodies                    |
	    const real,                        // I: global/maximum eps        |
	    const real      =Default::theta,   //[I: tolerance parameter]      |
	    const kern_type =Default::kernel,  //[I: type of softening kernel] |
#ifdef falcON_INDI                             // non-public version only:     |
	    const bool      =false,            //[I: use individual eps?]      |
#endif                                         //                              |
	    const real      =one,              //[I: constant of gravity]      |
	    const MAC_type  =theta_of_M,       //[I: type of MAC]              |
	    const int[4] =Default::direct      //[I: N_direct for gravity]     |
#ifdef falcON_SPH                              // SPH version only:            |
	   ,const int[3] =Default::SPHdirect   //[I: N_direct for SPH]         |
#endif
	    );
//                                                                             |
// This second version is designed for the use with C and FORTRAN, see  files  |
// falcON_C.h and falcON.f. Once a 'falcON' is constructed, one can no longer  |
// switch between using  bodies or arrays.  Similarly, one cannot  change the  |
// arrays  used,  but,  of course,  one can change the  data stored in  these  |
// arrays.  Using the  flags,  one can in this way  obtain almost any desired  |
// behaviour, see below.                                                       |
//                                                                             |
// The  first  argument(s)  specify  the sink and  source  properties of  the  |
// bodies.  Each  body  has  the  sink  properties:  position (x,y,z),  mass,  |
// softening length (for the case of individual softening lengths), and flag.  |
// The source  properties are:  acceleration (ax,ay,az), potential, and mass-  |
// density.  If global softening is used, a  NULL pointer may be used instead  |
// of an array  holding eps_i. The same applies to the arrays for density and  |
// potential: if a NULL pointer is given, they will never be used.             |
//                                                                             |
// All but the first argument of both constructors are identical:              |
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
// The class falcON has a destructor that will be called implicitly, whenever  |
// the scope of an object of type falcON terminates. The destructor can also   |
// be called explicitly, or via the delete command (if a tree was allocated    |
// by the new command).                                                        |
//                                                                             |
    ~falcON();                              //                                 |
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
//   2      2     don't load this body into the tree, ignore it                |
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
// With                                                                        |
//                                                                             |
#if(0) // code currently broken
    void re_grow(int const&,                       // I: N_cut                 |
		 int const& =Default::Ncrit,       //[I: Ncrit]                |
		 const vect* const& =0);           //[I: pre-set root center   |
#endif
//                                                                             |
// a new tree is grown aided by the old tree, if existent. N_cut is the size   |
// of an old-tree cell that will be re-grown from scratch.                     |
// Ideally: N_crit << N_cut << N.                                              |
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
// However, if the bodies have moved only little, the force compuation is      |
// hardly slowed down and re-using an old tree saves the costs for             |
// establishing a new one. Note that reuse() does not allow for any change in  |
// the flags indicating whether or not a body shall be ignored: changing this  |
// flag will have an effect only at the next call of grow().                   |
//                    +---------------------------+                            |
//                    | USE THIS OPTION AREFULLY! |                            |
//                    +---------------------------+                            |
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
// - performs all interaction (interaction phase), and                         |
// - evaluates the gravity for active bodies (evaluation phase).               |
// (Note that in earlier versions, the upward pass was already done at tree    |
// growth. The new version will therefore differ in the timing between tree    |
// growth and gravity approximation.)                                          |
//                                                                             |
    void approximate_gravity(bool const& =true,       //[I: combine phases]    |
			     bool const& =false       //[I: all or only active]|
#ifdef falcON_ADAP                                    // non-public  only:     |
			    ,real const& =zero,       //[I: Nsoft adjust eps_i]|
			     uint const& =0u,         //[I: Nref  adjust eps_i]|
			     real const& =zero,       //[I: eps_min]           |
			     real const& =zero        //[I: max change of eps] |
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
// containing that particle with not less than Nref [3rd argument] bodies).    |
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
// approximate Taylor-series based method. The elements of the array refer to  |
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
    void exact_gravity(bool const& =false   //[I: all of only active]          |
#ifdef falcON_ADAP                          // non-public  only:               |
		      ,real const& =zero,   //[I: Nsoft: adjust eps_i]         |
		       uint const& =0u,     //[I: Nref:  adjust eps_i]         |
		       real const& =zero,   //[I: eps_min]                     |
		       real const& =zero    //[I: max change of eps]           |
#endif                                      //                                 |
		       );                   //                                 |
//                                                                             |
//                                                                             |
// 3.3 A crude Estimation of the Mass- and Number-Density                      |
// ------------------------------------------------------                      |
// There is also the possibility to obtain a rough estimate of the mass-,      |
// surface- or number density of bodies in the neighbourhood of every body     |
// flagged being active, via                                                   |
//                                                                             |
    void estimate_rho(                      // estimate mass volume density    |
		      uint const&,          // I: critical cell size           |
		      bool const& = 0);     //[I: all or active only?]         |
    void estimate_sd (                      // estimate mass surface density   |
		      uint const&,          // I: critical cell size           |
		      bool const& = 0);     //[I: all or active only?]         |
    void estimate_n  (                      // estimate number volume density  |
		      uint const&,          // I: critical cell size           |
		      bool const& = 0);     //[I: all or active only?]         |
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
// Once a tree structure is established, you can use falcON to obtain          |
// estimates for SPH quantities. Currently however, this is restricted to      |
// Walter's private version of the code.                                       |
//                                                                             |
// This part of the code considers only bodies flagged as SPH particles.       |
//                                                                             |
//                                                                             |
// 4.1 Partner Counting                                                        |
// --------------------                                                        |
// The following routine counts the number of interaction partners.            |
//                                                                             |
    void count_SPH_partners(bool const& = 0,       //[I: all (or only active)?]|
			    bool const& = 1);      //[I: update bodies?]       |
//                                                                             |
// The 1st arg indicates whether all SPH-bodies shall be considered active.    |
// The 2nd arg indicates whether or not the snum() field of the bodies shall   |
// be updated.                                                                 |
//                                                                             |
//                                                                             |
// 4.2 Adjusting the smoothing sizes                                           |
// ---------------------------------                                           |
// The following routine iteratively adjusts the (smoothing) sizes such that   |
// the number of interaction partners averages to a given value. On return,    |
// the size() and snum() fields of the bodies are updated.                     |
//                                                                             |
    void adjust_SPH_sizes(uint const&,            // I: desired # partners     |
			  real const&,            // I: h_max for bodies       |
			  bool const& =0,         //[I: ALL (or only active)?] |
			  uint        =0u,        //[I: # partners: min allwd] |
			  uint        =0u,        //[I: # partners: max allwd] |
			  uint const& =10u);      //[I: max # iterations]      |
//                                                                             |
// The 1st arg gives the desired average number K of partners, the 2nd arg     |
// the maximum allowed size. The 3rd arg indicates whether all SPH-bodies      |
// shall be considered active. The 4th and 5th arg give the max and min        |
// numbers of partners tolerated; if omitted 2/3K and 3/2K are used. The last  |
// arg gives the maximum number of iterations.                                 |
//                                                                             |
//                                                                             |
// 4.3 Sweep One                                                               |
// -------------                                                               |
// The following routine computes the number of interaction partners and the   |
// gas-density for each active SPH body. Usually, the body data are not yet    |
// updated, but to be updated after sweep two (see below).                     |
//                                                                             |
    void SPH_sweep_one(const eq_of_state*const&,  // I: equation of state      |
		       uint              const&,  // I: <# partners/body>      |
		       bool           const& =0,  //[I: ALL (or only active)?] |
		       bool           const& =0); //[I: update bodies' ?]      |
//                                                                             |
// The 1st arg gives the equation of state (see inc/sph/eofs.h), the 2nd arg   |
// indicates the mean number of interaction partners (this is used to estimate |
// the memory requirements). The 3rd arg indicates whether all SPH bodies      |
// shall beconsidered active. The last arg indicates whether or not the data   |
// fields of the bodies shall be updated.                                      |
// The routines uses as input the bodies' mass(), vprd(), size() fields and,   |
// if the equation of state requires it, also the uprd() field. For every      |
// active SPH body it will estimate the snum() and srho() data, but update     |
// them only if explicitly desired by the last arg. For inactive bodies, the   |
// srho() data are read into the tree, to be used in sweep two.                |
// Simultaneously with the partner counting and density estimation, this       |
// routine also produces an interaction list, which will be used in sweep two. |
//                                                                             |
//                                                                             |
// 4.4 Sweep Two                                                               |
// -------------                                                               |
// The following routine, which can only be used AFTER SPH_sweep_one(),        |
// estimates drho/dt (or grad v) and dv/dt and du/dt due to pressure gradient  |
// and artificial viscosity. du/dt is computed only if required by the eq of   |
// state (not for isothermal gas). Upon return, the snum(), srho(), drho(),    |
// and, possibly, udin() data fields of active SPH-bodies are updated and      |
// dv/dt added to the acc() data field.                                        |
//                                                                             |
    void SPH_sweep_two();
//                                                                             |
//  Currently, the code only implements the traditional artificial viscosity   |
//  with eps=0.01.                                                             |
//                                                                             |
#endif                                      //                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
// 5 SEARCH FOR AND COUNTING OF INTERACTION PARTNERS                           |
// =================================================                           |
//                                                                             |
// Once a tree structure is established, you can also use it to create         |
// interaction lists for neighbouring or colliding particles (useful for,      |
// e.g., SPH and sticky-particle simulations). This is provided by             |
//                                                                             |
    typedef uint elem_pair[2];              // element of an interaction list  |
    void                                    //                                 |
    make_iaction_list(elem_pair *,          // I: interaction list             |
		      uint const&,          // I: allocated size of list       |
		      uint      &,          // O: actual size of list          |
		      bool const& = 1,      //[I: 1/0: Max/Sum (see below)]    |
		      real const& =-one);   //[I: time step tau]               |
//                                                                             |
// In case of overflow, i.e. if the number of pairs found exceeds the size     |
// (2nd arg) of the list (1st arg), the routine issues a warning to stderr,    |
// but does not truncate the search, rather it ceases copying pairs into the   |
// interaction list. In this case, the value returned for the actual size      |
// (3rd arg) exceeds the maximum size (2nd arg), but reflects the actual       |
// number of interactions found.                                               |
//                                                                             |
//                                                                             |
// 5.1 SPH partner search                                                      |
// ----------------------                                                      |
// In order to make a list of all pairs {i,j} of bodies/indices for which      |
//                                                                             |
//      (1) both are flagged to be SPH particles,                              |
// and  (2) at least one is flagged being active,                              |
// and  (3)    | x_i - x_j | < max(size_i,size_j)   IF Max==true               |
//          OR | x_i - x_j | < size_i + size_j      IF Max==false              |
//                                                                             |
// use falcON::make_iaction_list() whereby omitting the time step (last arg).  |
// In the case of array comnunication, provide an array with sizes of the      |
// bodies. It is the responsibility of the user to ensure that this array is   |
// properly allocated (it is sufficient to have entries for all bodies         |
// flagged as SPH particles).                                                  |
//                                                                             |
//                                                                             |
// 5.2 Sticky-particle partner search                                          |
// ---------------------------------------------------                         |
// In order to make a list of all pairs {i,j} of bodies/indices for which      |
//                                                                             |
//      (1) both are flagged to be sticky particles,                           |
// and  (2) at least one is flagged being active,                              |
// and  (3) | (x_i+t*v_i)-(x_j+t*v_j) | < size_i+size_j  with t in [0,tau],    |
//                                                                             |
// use falcON::make_iaction_list() with the time step (last arg) >= 0. In the  |
// case of array communication, provide arrays with sizes and velocity         |
// components of the bodies. It is the responsibility of the user to ensure    |
// that these arrays are properly allocated (it is sufficient to have entries  |
// for all bodies flagged as sticky).                                          |
//                                                                             |
//                                                                             |
// 5.3 SPH partner counting                                                    |
// ------------------------                                                    |
// The following method counts for each active body flagged as sph the number  |
// of interaction partners, i.e. bodies that are also flagged as sph and       |
// for which | x_i - x_j | < size_i + size_j                                   |
//                                                                             |
    void count_sph_partners();              //                                 |
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
// For the number of bodies used by falcON, use                                |
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
    uint BB_interactions() const;           // R: # of body-body IAs           |
    uint MB_interactions() const;           // R: # of many-body IAs           |
    uint CB_interactions() const;           // R: # of cell-body IAs           |
    uint CC_interactions() const;           // R: # of cell-cell IAs           |
    uint total_interactions() const;        // R: # of all IAs                 |
//                                                                             |
// To obtain some global tree properties, use                                  |
//                                                                             |
    const vect& root_center() const;        // R: center of root               |
    const real& root_radius() const;        // R: radius of root               |
    const int & root_number() const;        // R: # bodies in root             |
    const real& root_mass  () const;        // R: mass  of root                |
    const uint& root_depth () const;        // R: depth  of root               |
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
// Dehnen W. & Teuben P.J., 2002, in preparation                               |
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
    const uint& No_bodies     () const;           // R: # bodies               |
    oct_tree* const& my_tree  () const            // R: my tree                |
    { return TREE; }                              //                           |
//-----------------------------------------------------------------------------+
//                                                                             |
// private methods                                                             |
//                                                                             |
  private:                                        //                           |
    falcON                    ();                 // not implemented           |
    falcON                    (const falcON&);    // not implemented           |
    falcON& operator=         (const falcON&);    // not implemented           |
//-----------------------------------------------------------------------------+
//                                                                             |
// data members                                                                |
//                                                                             |
    mutable grav_stat     *STATS;                 // statistics                |
    const   bodies        *BODIES;                // sbodies to be used        |
    const   ebodies       *ARRAYS;                // abodies to be used        |
    int                    Ncrit;                 // Ncrit                     |
    oct_tree              *TREE;                  // tree to be used           |
    grav_mac              *GMAC;                  // theta(M)                  |
    grav_estimator        *GRAV;                  // gravity estimator         |
    stsp_estimator        *STSP;                  // partner search            |
#ifdef falcON_SPH                                 //                           |
    sph_estimator         *SPHT;                  // SPH estimator             |
#endif                                            //                           |
  };                                              // END: class falcON         |
}                                                 // END: namespace nbdy       |
//-----------------------------------------------------------------------------+
#include <public/falcON.cc>                       // include inline definitions|
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif                                            // falcON_included_falcON_h  |
