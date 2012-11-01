// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/octtree.h
///
/// \brief  methods for building and walking an octtree in 2D or 3D and some
///         support for tree walking and tree using algorithms.
///
/// \author Walter Dehnen
///
/// \date   2009-2012
///
/// \version May-2009 WD  first tested version
/// \version Oct-2009 WD  new design using indices for leaves and cells
/// \version Nov-2009 WD  removed redundant template parameter
/// \version Jan-2010 WD  renamed methods in TreeAccess; added leaf's parent
/// \version Jan-2010 WD  neighbour search methods
/// \version Feb-2010 WD  new initialisation: removed need for OctalTree::Dot
/// \version Mar-2010 WD  class FastNeighbourFinder
/// \version Apr-2010 WD  class TreeWalkAlgorithm, tree pruning
/// \version Apr-2010 WD  faster and memory-leaner tree-building algorithm
/// \version Jun-2010 WD  16-byte alignement, using geometry.h
/// \version Jul-2012 WD  allowed for periodic boundaries
/// \version Aug-2012 WD  completely new, OMP parallel in non-public part
///                       requires C++11
/// \version Sep-2012 WD  alternative memory layout using std::vector
/// \version Oct-2012 WD  one memory block for OctalTree and InteractionTree
/// \version Oct-2012 WD  define only one TreeWalker instantination
/// \version Oct-2012 WD  prepare for AVX (block sizes and data alignment)
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009-2012 Walter Dehnen
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
#ifndef WDutils_included_octtree_h
#define WDutils_included_octtree_h

#if __cplusplus < 201103L
# error need C++11
#endif

#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_sstream
#  include <sstream>
#  define WDutils_included_sstream
#endif
#ifndef WDutils_included_iomanip
#  include <iomanip>
#  define WDutils_included_iomanip
#endif
#ifndef WDutils_included_periodic_h
#  include <periodic.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_io_h
#  include <io.h>
#endif
#ifndef WDutils_included_h
#  include <meta.h>
#endif

#if defined(_OPENMP) && defined(WDutilsDevel)
#  define  OCTALTREE_USE_OPENMP
#endif

#ifdef OCTALTREE_USE_OPENMP
# ifndef WDutils_included_parallel_h
#  include <parallel.h>
# endif

// shall we store the domain id for every cell? (it is hardly ever needed and
// can be obtained in O(1) time anyway)
# define OCTALTREE_HAVE_D0
# undef  OCTALTREE_HAVE_D0
#endif

// shall the tree data be allocated in one big block of memory or in individual
// chunks (using std::vector) ?
#define OCTALTREE_DATA_IN_ONE_BLOCK
//#undef OCTALTREE_DATA_IN_ONE_BLOCK

// shall the InteractionTree hold its own data block or share one with its base,
// the OctalTree? (takes only effect #ifdef OCTALTREE_DATA_IN_ONE_BLOCK)
#define INTERACTIONTREE_HOLDS_OWN_BLOCK
#undef  INTERATIONTREE_HOLDS_OWN_BLOCK

#undef INTERACTIONTREE_USES_BASE_BLOCK
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
# ifndef INTERACTIONTREE_HOLDS_OWN_BLOCK
#  define INTERACTIONTREE_USES_BASE_BLOCK
# endif
#else
# undef INTERACTIONTREE_HOLDS_OWN_BLOCK
# ifndef WDutils_included_vector
#  define  WDutils_included_vector "<vector>"
#  include <vector>
# endif
#endif

// shall we support user data directly from the tree?
// (required in falcON.2)
#define OCTALTREE_HAVE_DATA
// #undef  OCTALTREE_HAVE_DATA

// shall we perform some extra careful asserts, in particular on access to
// data which may not be available?
#if !defined(NDEBUG) && !defined(OCTALTREE_CAREFUL)
# define OCTALTREE_CAREFUL 1
#elif !defined(OCTALTREE_CAREFUL)
# define OCTALTREE_CAREFUL 0
#endif

#if OCTALTREE_CAREFUL

# define WDutilsTreeAssert(Condition) WDutilsAssert(Condition)
# define WDutilsTreeAssertE(Condition) WDutilsAssertE(Condition);

#else

# define WDutilsTreeAssert(Condition)
# define WDutilsTreeAssertE(Condition)

#endif

#ifdef WDutils_included_partree_h
# error "utils/octtree.h and utils/partree.h are incompatible"
#endif

////////////////////////////////////////////////////////////////////////////////
namespace { template<typename> struct TreeImplementer; }

namespace WDutils {
  template<typename, typename Enable=void> class TreeWalker;
  ///
  /// contains oct-tree-related functionality
  ///
  namespace octtree {
    template< template<typename,typename> class temp_type >
    struct temp_type_void_conv
    {
      template<typename _type>
      temp_type_void_conv(const temp_type<_type,void>&);
    };
  } // namespace octtree
  /// is a type @c Walker an instance of @c TreeWalker<>
  ///
  /// @c is_TreeWalker<Walker>::value is true iff @c Walker is (derived from a)
  /// @c TreeWalker<>
  template <class Walker> struct is_TreeWalker
  {
    static const bool value = 
      std::is_convertible<Walker,
			  octtree::temp_type_void_conv<TreeWalker> > :: value;
  };
  //
  namespace octtree {
    /// \name common types and constants used in trees
    //@{
    typedef uint32_t particle_key;
    typedef uint32_t count_type;
    typedef uint8_t depth_type;
    typedef uint8_t octant_type;
    typedef uint8_t local_count;
    typedef positional_offset_bits offbits;
#ifdef OCTALTREE_USE_OPENMP
    typedef uint8_t domain_id;
#endif// OCTALTREE_USE_OPENMP
    //@}
    template<int> struct block_flag_t;
    template<> struct block_flag_t<2> { typedef uint16_t type; };
    template<> struct block_flag_t<4> { typedef uint32_t type; };
    template<> struct block_flag_t<8> { typedef uint64_t type; };
    ///
    /// wrapper around count_type: either a cell or leaf in the tree
    ///
    template<bool _EXT, bool _CELL,
	     class = typename std::enable_if<!_EXT || !_CELL>::type>
    struct TreeNode {
      static const bool is_external =  _EXT;
      static const bool is_internal = !_EXT;
      static const bool is_cell     =  _CELL;
      static const bool is_leaf     = !_CELL;
      //
      count_type I;               ///< index within domain (can be changed)
      /// default
      constexpr TreeNode() = default;
      /// ctor
      constexpr explicit TreeNode(count_type i) noexcept : I(i) {}
      /// copy ctor
      constexpr TreeNode(TreeNode const&) = default;
      /// copy operator
      TreeNode&operator=(TreeNode const&) = default;
      /// prefix increment
      TreeNode&operator++() noexcept { ++I; return*this; }
      /// prefix decrement
      TreeNode&operator--() noexcept { --I; return*this; }
      /// postfix increment
      TreeNode operator++(int) noexcept { return TreeNode(I++); }
      /// postfix decrement
      TreeNode operator--(int) noexcept { return TreeNode(I++); }
      /// add: large increment
      TreeNode&operator+=(count_type i) noexcept { I+=i; return*this; }
      /// add: large decrement
      TreeNode&operator-=(count_type i) noexcept { I-=i; return*this; }
      /// obtain incremented
      TreeNode operator+ (count_type i) const noexcept { return TreeNode(I+i); }
      /// obtain decremented
      TreeNode operator- (count_type i) const noexcept { return TreeNode(I-i); }
      /// next node
      constexpr TreeNode next() const noexcept { return TreeNode(I+1); }
      /// is leaf K-algined?
      constexpr bool is_aligned(int K) const noexcept { return (I%K)==0; }
      /// alignment: last leaf that is K-aligned
      template<int K>
      constexpr TreeNode aligned() const noexcept
      {
	static_assert(K>0 && (K&(K-1))==0,"K not power of 2");
	return TreeNode(I&(~(K-1)));
      }
      /// \name comparisons
      //@{
      bool operator==(const TreeNode n) const noexcept { return I==n.I; }
      bool operator!=(const TreeNode n) const noexcept { return I!=n.I; }
      bool operator<=(const TreeNode n) const noexcept { return I<=n.I; }
      bool operator>=(const TreeNode n) const noexcept { return I>=n.I; }
      bool operator< (const TreeNode n) const noexcept { return I< n.I; }
      bool operator> (const TreeNode n) const noexcept { return I> n.I; }
      //@}
      /// an invalid node
      static TreeNode invalid() noexcept
      { return TreeNode(~count_type(0)); }
    };
    //
    typedef TreeNode<0,0> Leaf;              ///< type: internal leaf
    typedef TreeNode<1,0> ExtLeaf;           ///< type: external leaf
    typedef TreeNode<0,2> Cell;              ///< type: cell
    /// an invalid internal leaf
    inline Leaf invalid_leaf() noexcept
    { return Leaf::invalid(); }
    /// an invalid cell
    inline Cell invalid_cell() noexcept
    { return Cell::invalid(); }
    /// an invalid external leaf
    inline ExtLeaf invalid_extleaf() noexcept
    { return ExtLeaf::invalid(); }
    ///
    /// a range of 2^8 leaves
    ///
    struct LeafRange
    {
      enum {
	Shft = 8,                      ///< leaf -> range shift
	Size = 1<<Shft,                ///< # leaves / range
	Mask = Size-1                  ///< mask
      };
      //
      count_type I;                    ///< range index
      /// ctor
      explicit LeafRange(count_type i) : I(i) {}
      /// default ctor
      LeafRange(LeafRange const&) = default;
      /// default assignment
      LeafRange&operator=(LeafRange const&) = default;
      /// increment: become the next range, may be in another domain
      LeafRange&operator++() { ++I; return*this; }
      /// \name comparisons
      //@{
      bool operator==(const LeafRange r) const noexcept { return I==r.I; }
      bool operator!=(const LeafRange r) const noexcept { return I!=r.I; }
      bool operator<=(const LeafRange r) const noexcept { return I<=r.I; }
      bool operator>=(const LeafRange r) const noexcept { return I>=r.I; }
      bool operator< (const LeafRange r) const noexcept { return I< r.I; }
      bool operator> (const LeafRange r) const noexcept { return I> r.I; }
      //@}
      /// an invalid range
      static LeafRange invalid() noexcept
      { return LeafRange(~count_type(0)); }
      /// # ranges given the # leaves
      static count_type n_range(count_type nl) noexcept
      { return (nl+Mask) >> Shft; }
    private:
      LeafRange() = delete;
    };// struct LeafRange
#ifdef OCTALTREE_USE_OPENMP
    template<OMP::Schedule> struct _default_chunk {
      static unsigned size(unsigned n)
      { return std::max(1u,n/(32*OMP::TeamSize())); }
    };
    template<> struct _default_chunk<OMP::Static> {
      // note: 0 translates to non-overlapping static schedule in OMP::for_each
      static unsigned size(unsigned) { return 0; }
    };
    /// default chunk size: 0 for static n/(32*#threads) otherwise
    template<OMP::Schedule schedule>
    inline unsigned default_chunk(unsigned n)
    { return _default_chunk<schedule>::size(n); }
#endif
    /// collection of common types and constants used with trees
    struct static_tree_properties {
      using Leaf = octtree::Leaf;
      using Cell = octtree::Cell;
      template<bool EXT>
      using TreeLeaf = TreeNode<EXT,0>;
      using LeafRange = octtree::LeafRange;
      using ExtLeaf = octtree::ExtLeaf;
      using count_type = octtree::count_type;
      using particle_key = octtree::particle_key;
      using offbits = octtree::offbits;
      using depth_type = octtree::depth_type;
      using octant_type = octtree::octant_type;
      using local_count = octtree::local_count;
#ifdef OCTALTREE_USE_OPENMP
      using domain_id = octtree::domain_id;
#endif
    };
    /// block size for leaf data
#ifdef __AVX__
    constexpr count_type LeafBlockSize = 8;
#else
    constexpr count_type LeafBlockSize = 4;
#endif
    /// data alignement (bytes) used for all tree data in OctalTree<>
    constexpr count_type DataAlignment= 4*LeafBlockSize;
    /// the next multiple of @a DataAlignment to n*sizeof(T)
    template<typename T>
    inline constexpr size_t NextAligned(size_t n)
    { return WDutils::next_aligned<DataAlignment>(n*sizeof(T)); }
    /// useful as base struct for ensuring sensible block size
    template<int _BlockSize=4>
    struct static_assert_block_size
    {
      static const int BlockSize = _BlockSize;
      static_assert(0
#ifdef __SSE2__
		    || BlockSize==2
#endif
#ifdef __SSE__
		    || BlockSize==4
#endif
#ifdef __AVX__
		    || BlockSize==8
#endif
		    ,"Blocksize must be one of ["
#define _COLON
#ifdef __SSE2__
		    "2"
# undef  _COLON
# define _COLON ":"
#endif
#ifdef __SSE__
		    _COLON "4"
# undef  _COLON
# define _COLON ":"
#endif
#ifdef __AVX__
		    _COLON "8"
#endif
#undef  _COLON
		    "]");
    };
    /// std::vector with K-byte aligning allocator
    template<typename _Tp> using aligned_vector =
      std::vector<_Tp, AlignmentAllocator<_Tp, DataAlignment> >;
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    template<typename _Tp> using storage = _Tp*;
    template<typename _Tp>
    using aligned_storage = _Tp*;
    template<typename _Tp>
    inline bool have_data(const _Tp*X) { return X != 0; }
    template<typename _Tp>
    inline const _Tp*cp_data(const _Tp*X) { return X; }
#else // OCTALTREE_DATA_IN_ONE_BLOCK
    template<typename _Tp>
    using storage         = std::vector<_Tp>;
    template<typename _Tp>
    using aligned_storage = aligned_vector<_Tp>;
    template<typename _Tp>
    inline bool have_data(storage<_Tp> const&X) { return X.size() > 0; }
    template<typename _Tp> inline bool have_data(aligned_storage<_Tp> const&X)
    { return X.size() > 0; }
    template<typename _Tp>
    inline const _Tp*cp_data(aligned_storage<_Tp,_alignment> const&X)
    { return X.data(); }
#endif// OCTALTREE_DATA_IN_ONE_BLOCK
  } // namespace octtree
  //
  template<> struct traits< octtree::Leaf >
  { static const char*name() { return "octtree::Leaf"; } };
  template<> struct traits< octtree::ExtLeaf >
  { static const char*name() { return "octtree::ExtLeaf"; } };
  template<> struct traits< octtree::Cell >
  { static const char*name() { return "octtree::Cell"; } };
  template<> struct traits< octtree::LeafRange >
  { static const char*name() { return "octtree::LeafRange"; } };
  /// \name global functions on octtree::TreeNode<>
  /// \relates @c OctalTree
  //@{
  /// print name of external leaf
  /// \relates @c octtree::TreeNode<1,0>
  inline std::ostream& operator<<(std::ostream&s, const octtree::ExtLeaf l)
  {
    if(l == octtree::invalid_extleaf())
      return s << "E.invalid";
    return s << "E." 
	     << std::setw(7) << std::setfill('0') << l.I 
	     << std::setfill(' ');
  }
  /// print name of internal leaf
  /// \relates @c octtree::TreeNode<0,0>
  inline std::ostream& operator<<(std::ostream&s, const octtree::Leaf l)
  {
    if(l == octtree::invalid_leaf())
      return s << "L.invalid";
    return s << "L." 
	     << std::setw(7) << std::setfill('0') << l.I 
	     << std::setfill(' ');
  }
  /// print name of cell
  /// \relates @c octtree::TreeNode<0,1>
  inline std::ostream& operator<<(std::ostream&s, const octtree::Cell c)
  {
    if(c == octtree::invalid_cell())
      return s << "C.invalid";
    return s << "C." 
	     << std::setw(7) << std::setfill('0') << c.I 
	     << std::setfill(' ');
  }
  /// print name of leaf range
  /// \relates @c octtree::LeafRange
  inline std::ostream& operator<<(std::ostream&s, const octtree::LeafRange r)
  {
    if(r == octtree::LeafRange::invalid())
      return s << "R.invalid";
    return s << "R." 
	     << std::setw(7) << std::setfill('0') << r.I 
	     << std::setfill(' ');
  }
  /// name of TreeNode<> as C++ string
  /// \relates octtree::TreeNode
  template<bool E, bool C>
  inline std::string Name(const octtree::TreeNode<E,C> n)
  {
    std::ostringstream s;
    s << n;
    return s.str();
  }
  /// name of TreeNode<> as C string (for use in DebugInfo, Error, and Warning)
  template<bool E, bool C>
  inline const char*name(const octtree::TreeNode<E,C> n)
  { 
    return Name(n).c_str();
  }
  //@}
  /// useful as base class to ensure template parameters are suitable for @c
  /// OctalTree<>
  template<int _Dim, typename _PosType>
  class parameters_okay_for_OctalTree
  {
    static_assert(is_in_range<_Dim,2,3>::value,
		  "OctalTree<>: only 2D and 3D is allowed");
    static_assert(std::is_floating_point<_PosType>::value,
		  "OctalTree<>: position type must be floating point");
  };
  ///
  /// A spatial tree of square (2D) or cubic (3D) cells
  ///
  template<int _Dim, typename _PosType>
  class OctalTree :
    private parameters_okay_for_OctalTree<_Dim,_PosType>,
    public  octtree::static_tree_properties
  {
    template<typename, typename> friend class TreeWalker;
    friend struct ::TreeImplementer<OctalTree>;
    // disable default and copy ctor
    OctalTree           (const OctalTree&) = delete;
    OctalTree& operator=(const OctalTree&) = delete;
  public:
    /// \name public constants and types
    //@{
    static const int                 Dim = _Dim;    ///< # dimensions
    static const octtree::depth_type Nsub= 1<<Dim;  ///< # octants/cell
    static const int                 Noff= Nsub-1;  ///< Nsub-1
    /// maximum for @a NMAX
    static const octtree::depth_type MAXNMAX = 250;
    /// block size of leaf data
    static const count_type LeafBlockLim  = octtree::LeafBlockSize-1;
    static const count_type LeafBlockMask = ~LeafBlockLim;
    /// smallest multiple of LeafBlockSize larger or equal to n
    static count_type blocked(const count_type n) noexcept
    { return (n+LeafBlockLim)&LeafBlockMask; }
    //
    typedef _PosType pos_type;
    typedef Geometry::cube<Dim,pos_type> cube;
    typedef SSE::Extend16<cube> cube16;
    typedef typename cube::point point;
    typedef SSE::Extend16<point> point16;
    typedef PeriodicBox<Dim,pos_type> PerBoundary;
    template<typename _Tp> using storage        = octtree::storage<_Tp>;
    template<typename _Tp> using aligned_storage= octtree::aligned_storage<_Tp>;
    /// the maximum allowed tree depth
    /// \note Set to the number of binary digits in our floating point type.
    static const depth_type MaximumDepth=std::numeric_limits<pos_type>::digits;
    /// represents a particle during domain decomposition
    struct _Dot {
      point          X;                             ///< position
      particle_key   I;                             ///< particle identifier
      mutable void  *Next;                          ///< next dot in a list
    };
    typedef SSE::Extend16<_Dot> Dot;               ///< 16-aligned dot
    //@}
    ///
    /// Interface for initialising particles at tree building
    ///
    /// The interface here is intended to allow (1) for a general particle
    /// data layout in the user application, (2) for OMP parallelism during
    /// initialisation of particle positions, (3) removal and creation of
    /// particles, and (4) efficient implementation. These requirements seem
    /// to make it necessary to trust the user, in particular, with the
    /// particle order when re-initialisating the particle data.
    struct Initialiser {
      /// dtor
      virtual ~Initialiser() {}
      /// provide periodic boundary, if any
      virtual const PerBoundary*Boundary() const = 0;
      /// invalid entry for particle_key: all bits set
      static const particle_key InvalidKey = ~0;
      /// # external particles
      virtual count_type Nexternal() const = 0;
      /// initialise a external particle to load
      /// \param[in]  i  running index of external particle
      /// \param[out] x  position
      /// \param[out] k  particle key
      /// \note called @a Nexternal() times at construction or re-build.
      virtual void InitExtern(count_type i, point&x, particle_key&k)
	const = 0;
      /// 
      /// initialises @a Dot::X and @a Dot::I for all internal particles
      /// 
      /// \return            dots initialised (and allocated), must be != 0
      /// \note will be de-allocated using WDutils_DEL16
      ///
      /// \param[out] Ndot   # dots initialised (and allocated), must be > 0
#ifdef OCTALTREE_USE_OPENMP
      /// \param[in]  iT     local thread index
      /// \param[in]  nT     # threads
#endif
      /// \param[in]  Keys   if!=0 : global table of particle keys from old tree
      /// \param[in]  Nkey   size of table @a Keys
      /// \note In tree construction Keys==0, while in rebuild() Keys!=0
#ifdef OCTALTREE_USE_OPENMP
      ///       both times from within a parallel region
#endif
      /// \note In case of @a Keys!=0, a possible implementation if no particles
      ///       have been removed or created since the last tree build is as
      ///       follows (note this keeps the particle order, resulting in
      ///       fastest tree re-building)
      /// \code
      ///       count_type k;
      ///       OMP::PartitionSub(Nkey,iT,nT,k,Nloc);
      ///       Dot*const D0 = WDutils_NEW16(Dot,Nloc);
      ///       for(Dot*Di=D0,*DN=D0+Nloc; Di!=DN; ++k,++Di) {
      ///         Di->I = Keys[k];
      ///         Di->X = my_particle_position(Keys[k]);
      ///       }
      ///       return D0;
      /// \endcode
      ///       But if particles have been removed or created, the user has to
      ///       find a way to ensure that all particles are initialised.
      /// 
      /// 
      /// \note the implementation can assume that each Dot is 16-byte
      ///       aligned.
      virtual Dot*InitIntern(count_type&Nloc, 
#ifdef OCTALTREE_USE_OPENMP
			     depth_type iT, depth_type nT,
#endif
			     const particle_key*Keys, count_type Nkey)
	const=0;
      /// shall we load rungs and flags?
      virtual bool LoadRungFlag() const
      { return false; }
      /// load rungs and flags for all particles in a domain
      /// \param[in]  key    table: particle identifiers
      /// \param[out] rung   table: rung for particles associated with keys
      /// \param[out] flag   table: flag for particles associated with keys
      /// \param[in]  n      size of tables.
      /// \return            all particles flagged as active?
      virtual bool InitRungFlag(const particle_key*,
				float*rung, uint8_t*flag, count_type n) const
      {
	for(count_type i=0; i!=n; ++i) {
	  rung[i] = 0;
	  flag[i] = 0xff;
	}
	return true;
      }
    };// struct OctalTree<>::Initialiser 
#ifdef OCTALTREE_USE_OPENMP
    ///
    /// holding domain related data
    ///
    struct DomainData {
      std::vector<count_type> BRANCH;     ///< branch cells
      count_type              C0,CN,NC;   ///< cells: begin, end, number
      count_type              L0,LN,NL;   ///< leaves: begin, end, number
      depth_type              DEPTH;      ///< maximum depth of any branch
      /// depth of domain
      depth_type const&depth() const noexcept 
      { return DEPTH; }
      /// # branches in domain
      count_type n_branch() const noexcept 
      { return BRANCH.size(); }
      /// branch cell of given index within domain
      /// \note branch cells are in the top-tree
      Cell branch_no(count_type b) const
      { WDutilsAssert(b<n_branch()); return Cell(BRANCH[b]); }
      /// # leaves in domain
      count_type const&n_leaf() const noexcept
      { return NL; }
      /// begin of domain leaves
      Leaf begin_leaf() const noexcept
      { return Leaf(L0); }
      /// end of domain
      Leaf end_leaf() const noexcept
      { return Leaf(LN); }
      /// # non-branch cells in domain
      count_type const&n_nonbranch_cell() const noexcept
      { return NC; }
      /// begin of non-branch cells in domain
      Cell begin_nonbranch_cell() const noexcept
      { return Cell(C0); }
      /// end of non-branch cells in domain
      Cell end_nonbranch_cell() const noexcept
      { return Cell(CN); }
      /// last non-branch cell in domain
      Cell rbegin_nonbranch_cell() const noexcept
      { return -- end_nonbranch_cell(); }
      /// before first non-branch cell in domain
      Cell rend_nonbranch_cell() const noexcept
      { return -- begin_nonbranch_cell(); }
      /// loop domain leaves
      /// \param[in] f    function called for each leaf
      template<typename FuncOfLeaf>
      void loop_leaves(FuncOfLeaf f) const noexcept(noexcept(f))
      { for(auto l=begin_leaf(); l!=end_leaf(); ++l) f(l); }
      /// loop branch cells
      /// \param[in] f    function called for each branch cell in domain
      template<typename FuncOfCell>
      void loop_branches(FuncOfCell f) const noexcept(noexcept(f))
      { for(count_type b=0; b!=n_branch(); ++b) f(Cell(BRANCH[b])); }
      /// loop non-branch cells down
      /// \param[in] f    function called for each non-branch cell in domain
      template<typename FuncOfCell>
      void loop_nonbranch_cells_down(FuncOfCell f) const noexcept(noexcept(f))
      {
	for(Cell c=begin_nonbranch_cell(); c!=end_nonbranch_cell(); ++c)
	  f(c);
      }
      /// loop non-branch domain cells up
      /// \param[in] f    function called for each non-branch cell in domain
      template<typename FuncOfCell>
      void loop_nonbranch_cells_up(FuncOfCell f) const noexcept(noexcept(f))
      {
	for(Cell c=rbegin_nonbranch_cell(); c!=rend_nonbranch_cell(); --c)
	  f(c);
      }
      /// loop domain cells down (branch cells first), including branch cells
      /// \param[in] f    function called for each cell in domain
      /// \note useful for domain-parallel downward pass
      template<typename FuncOfCell>
      void loop_cells_down(FuncOfCell f) const noexcept(noexcept(f))
      {
	loop_branches(f);
	loop_nonbranch_cells_down(f);
      }
      /// loop domain cells up (branch cells last), including branch cells
      /// \param[in] f    function called for each cell in domain
      /// \note useful for domain-parallel upward pass
      template<typename FuncOfCell>
      void loop_cells_up(FuncOfCell f) const noexcept(noexcept(f))
      {
	loop_nonbranch_cells_up(f);
	loop_branches(f);
      }
    };// struct OctalTree<>::DomainData
    typedef const DomainData*cp_domain;
#endif
   private:
    /// \name data
    //@{
    const Initialiser*const            INIT;    ///< initialiser
    const PerBoundary*const            PERB;    ///< periodic boundary, if any
#ifdef OCTALTREE_USE_OPENMP
    const domain_id                    NDOM;    ///< # domains
    const DomainData *const            DOM;     ///< domain data
    count_type                         NTOPC;   ///< # top cells
#endif// OCTALTREE_USE_OPENMP
    count_type                         NLEAFEXT;///< # external leaves
    count_type                         NLEAF;   ///< # internal leaves
    count_type                         NCELL;   ///< # cells
    count_type                         NBUILD;  ///< # builds
    depth_type                         NMAX;    ///< N_max
    depth_type                         NMIN;    ///< N_min
#ifdef OCTALTREE_USE_OPENMP
    count_type                         TOL;     ///< tolerance at last build
#endif// OCTALTREE_USE_OPENMP
    const bool                         ASCC;    ///< avoid single-child cells?
    // tree data
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    raw_memory<octtree::DataAlignment> _BUF;    ///< memory holding tree data
#endif// OCTALTREE_DATA_IN_ONE_BLOCK
  private:
    // data for leaves
    aligned_storage<point16>           _XL[2];  ///< any leaf: position
    storage<particle_key>              _PL[2];  ///< any leaf: particle
    aligned_storage<float>             _RG[2];  ///< any leaf: rung
    aligned_storage<uint8_t>           _FL[2];  ///< any leaf: flag
    storage<count_type>                _PC;     ///< int leaf: parent cell
    // data for cells
    storage<depth_type>                _LE;     ///< cell: tree level
    storage<octant_type>               _OC;     ///< cell: octant in cell
    aligned_storage<cube16>            _XC;     ///< cell: cubic box
    storage<count_type>                _L0;     ///< cell: first leaf
    storage<local_count>               _NL;     ///< cell: # leaf kids
    storage<count_type>                _NM;     ///< cell: # leaves in total
    storage<count_type>                _C0;     ///< cell: index: 1st daughter
    storage<octant_type>               _NC;     ///< cell: # daughter cells
    storage<count_type>                _PA;     ///< cell: parent cell index
    storage<count_type>                _NA;     ///< cell: # active leaves
    storage<depth_type>                _DP;     ///< cell: tree depth
#ifdef OCTALTREE_USE_OPENMP
    storage<domain_id>                 _D0;     ///< cell: first domain
#endif// OCTALTREE_USE_OPENMP
    //@}
#ifdef OCTALTREE_HAVE_DATA
    static const int     NDAT=16;       ///< # pters to user data
    void                *DATA[2][NDAT]; ///< pters to user data
  public:
    static const int     FREE=0;        ///< first free slot in DATA[]
    /// set OctalTree::DATA[i] to given pointer
    /// \param[in] p  pointer to set data to
    template<bool EXT, int I>
    void SetData(void*p) noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::SetData<EXT,I>: I out of range");
      const_cast<void*&>(DATA[EXT][I]) = p;
    }
    /// return pointer set by SetData<EXT,I>
    template<bool EXT, int I>
    void*GetData() const noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::GetData<EXT,I>: I out of range");
      return DATA[EXT][I];	    
    }
    ///
# if OCTALTREE_CAREFUL
#  define ASSERT_DAT_OKAY(FUNC)						\
    if(DATA[EXT][I] == 0)						\
      WDutils_ErrorN("invalid " FUNC "<EXT=%d,I=%d,T=%s>(n=%d)\n",	\
		     EXT,I,nameof(T),n);
#  define ASSERT_VEC_OKAY(FUNC)						\
    if(DATA[EXT][I] == 0)						\
      WDutils_ErrorN("invalid " FUNC "<EXT=%d,I=%d,J=%d,T=%s>(n=%d)\n",	\
		     EXT,I,J,nameof(T),n);
# else
#  define ASSERT_DAT_OKAY(FUNC)
#  define ASSERT_VEC_OKAY(FUNC)
# endif// OCTALTREE_CAREFUL
    /// const access to DATA[I][i]
    template<bool EXT, int I, typename T>
    T const&dat(count_type n) const noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::dat<I>: I out of range");
      ASSERT_DAT_OKAY("dat");
      return static_cast<const T*>(DATA[EXT][I])[n];
    }
    /// non-const access to DATA[I][i]
    template<bool EXT, int I, typename T>
    T&dat_lv(count_type n) const noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::dat_lv<I>: I out of range");
      ASSERT_DAT_OKAY("dat)lv");
      return static_cast<T*>(DATA[EXT][I])[n];
    }
    /// const access to user vector datum
    template<bool EXT, int I, int J, typename T> 
    T const&vec(count_type n) const noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::vec<I,J>: I out of range");
      static_assert(0<=J && J<Dim, "OctalTree::vec<I,J>: J out of range");
      ASSERT_VEC_OKAY("vec");
      return static_cast<const T*const*>(DATA[EXT][I])[J][n];
    }
    /// non-const access to user vector datum
    template<bool EXT, int I, int J, typename T> 
    T&vec_lv(count_type n) const noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::vec_lv<I,J>: I out of range");
      static_assert(0<=J && J<Dim, "OctalTree::vec_lv<I,J>: J out of range");
      ASSERT_VEC_OKAY("vec_lv");
      return static_cast<T*const*>(DATA[EXT][I])[J][n];
    }
    /// pointer to user datum
    template<bool EXT, int I, typename T>
    T*p_dat(count_type n) const noexcept
    {
      static_assert(0<=I && I<NDAT,"OctalTree::p_dat<I>: I out of range");
      ASSERT_DAT_OKAY("p_dat");
      return static_cast<T*>(DATA[EXT][I])+n;
    }
    /// constant pointer to user datum
    template<bool EXT, int I, typename T>
    const T*cp_dat(count_type n) const
    {
      static_assert(0<=I && I<NDAT,"OctalTree::cp_dat<I>: I out of range");
      ASSERT_DAT_OKAY("cp_dat");
      return static_cast<const T*>(DATA[EXT][I])+n;
    }
    /// non-const pointer to vector component
    template<bool EXT, int I, int J, typename T>
    T*p_vec(count_type n) const
    {
      static_assert(0<=I && I<NDAT,"OctalTree::p_vec<I,J>: I out of range");
      static_assert(0<=J && J<Dim, "OctalTree::p_vec<I,J>: J out of range");
      ASSERT_VEC_OKAY("p_vec");
      return static_cast<T*const*>(DATA[EXT][I])[J]+n;
    }
    /// const pointer to vector component from free methods
    template<bool EXT, int I, int J, typename T>
    const T*cp_vec(count_type n) const
    {
      static_assert(0<=I && I<NDAT,"OctalTree::cp_vec<I,J>: I out of range");
      static_assert(0<=J && J<Dim, "OctalTree::cp_vec<I,J>: J out of range");
      ASSERT_VEC_OKAY("cp_vec");
      return static_cast<const T* const*>(DATA[EXT][I])[J]+n;
    }
# undef ASSERT_DAT_OKAY
# undef ASSERT_VEC_OKAY
#endif// OCTALTREE_HAVE_DATA
    /// \name tree building and related
    //@{
    /// data required in ctor.
    struct Data {
      depth_type NMAX;       ///< maximum # particles in unsplit octants
      depth_type NMIN;       ///< minimum # particles in cell
      bool       ASCC;       ///< avoid single-child cells?
      count_type TOL;        ///< tolerance for # leaves in domains
      Data()
	: NMAX( Nsub )
	, NMIN( 0 )
	, ASCC( true )
	, TOL ( 0 ) {}
      Data(Data &&) = default;
      Data(Data const&) = default;
      Data&operator=(Data &&) = default;
      Data&operator=(Data const&) = default;
    };
    /// \brief
    /// ctor: build tree from scratch
    ///
    /// \details
    /// The tree is build in several stages as follows.\n
    ///
    /// - The particle positions and keys for all tree leaves are initialised
#ifdef OCTALTREE_USE_OPENMP
    ///   from inside an OMP parallel region.
#endif
    ///
    /// - The root box is established and a 'box-dot' is built
#ifdef OCTALTREE_USE_OPENMP
    ///   by each thread from its share of all particles.
    ///
    /// - The global tree is decomposed into domains with (roughly) equal
    ///   number of particles along a space-filling curve. \n
    /// 
    /// - Each thread establishes its domain part of the global 'box-dot' tree
    ///   by merging the contributions from other threads. \n
#endif
    /// - The cell-leaf tree is established by mapping boxes and octants with at
    ///   least @a nmin particles to cells and dots to leaves such that any
    ///   cell's leaf and cell descendants are contiguous in memory.
    ///
#ifdef OCTALTREE_USE_OPENMP
    /// OMP parallelism is used in all stages, in fact all stages are within a
    /// single OMP parallel region. The number of domains equals the number of
    /// parallel threads and is given by @a RunInfo::omp_threads().
    ///
    /// \note Must not be called from within a parallel region.
#endif
    ///
    /// \param[in] _init   initialiser used to set particle data
    ///
    /// \note If periodic boundaries are used, as specified by
    ///       _init->Boundary(), and _init->CheckBoundary() returns true, then
    ///       we check each particle position for conformity to the boundary,
    ///       issue a warning if we find non-conforming positions, and use the
    ///       conformed position instead.
    ///
    /// \param[in] _data   data specifying tree
    ///
    /// \note We require @a _data.NMIN <= @a _data.NMAX <= 250 but map @a
    ///       _data.NMIN=0 to @a _data.NMIN=min(2,_data.NMAX). For @a
    ///       _data.NMIN=2 (and @a _data.NMAX>1) the cell-leaf tree reflects the
    ///       depth of the original box-dot tree. However, for @a _data.NMIN >
    ///       2, while the cell-leaf tree is not as deep, the tree order of the
    ///       deeper box-dot tree is preserved in the tree order of the
    ///       leaves. For @a _data.NMAX=_data.NMIN=1 the tree is build to
    ///       maximum depth when each leaf has its own final cell (and tree
    ///       building is most CPU time consuming), which is unlikely to be
    ///       required by most applications.
    ///
#ifdef OCTALTREE_USE_OPENMP
    /// \note _data.NMAX must not exceed the number of particles per domain.
    ///       We try to adjust _data.TOL (see below) to help meeting this
    ///       criterion, but we may not succeed.
    ///
    /// \note If @a _data.TOL==0 (default), we choose a reasonable value. Use @a
    ///       _data.TOL=1 if you want no tolerance, i.e. an completely even
    ///       distribution of leaves accross domains. Note that @a _data.TOL=1
    ///       enforces deep domain splits and hence domains with many small
    ///       branches at their edges.
#endif
    ///
    /// \note Single-child cells occur if all their particles are in just one
    ///       octant. If @a _data.ASCC=true (the default) such cells are
    ///       eliminated in favour of their only child cell. In this case,
    ///       parent and child cells may be more than one tree level apart and
    ///       the tree depth less than the maximum tree level of any cell.
    ///
    /// \note It is strongly recommended to keep @a _data.NMIN, and @a _data.TOL
    ///       at their default setting unless you have very specific reasons.
    OctalTree(const Initialiser*_init, Data const&_data) throw(exception);
    ///
    /// rebuild tree (after particles have changed position and/or number)
    ///
    /// The tree is build exactly in the same way as with the constructor,
    /// except for the initial @a Dot order, which is taken to be that of the
    /// old tree.
#ifdef OCTALTREE_USE_OPENMP
    ///
    /// \note Must not be called from within a parallel region.
#endif
    ///
    /// \param[in] _data  data specifying tree
    ///
    /// \note If @a _data.NMAX=0 we take the old values for both @a _data.NMIN
    ///       and @a _data.NMAX. Otherwise, we require @a _data.NMIN, @a
    ///       _data.NMAX <=250 and map @a _data.NMIN=0 to @a
    ///       _data.NMIN=min(2,_data.NMAX).
    ///
#ifdef OCTALTREE_USE_OPENMP
    ///
    /// \note If @a _data.TOL==0 (default), we choose a reasonable value. Use @a
    ///       _data.TOL=1 if you want no tolerance, i.e. an completely even
    ///       distribution of leaves accross domains. Note that @a _data.TOL=1
    ///       enforces deep domain splits and hence domains with many small
    ///       branches at their edges.
#endif
    void rebuild(Data const&_data) throw(exception);
    /// dtor
#ifdef OCTALTREE_USE_OPENMP
    /// \note Must not be called from within a parallel region.
#endif// OCTALTREE_USE_OPENMP
    virtual ~OctalTree();
    //@}
    ///
    /// interface supporting building of derived and dependent classes
    ///
    struct BuilderInterface {
      /// dtor
      virtual ~BuilderInterface() {}
      /// typically allocates global data
      /// \note called from one thread
      /// \param[in] fresh    is this a construction (or a re-build)?
      /// \param[in] ncell    # cells
      /// \param[in] nleaf    # leaves
      /// \param[in] nextleaf # external leaves
      virtual void Before(bool fresh, count_type ncell,
			  count_type nleaf, count_type nextleaf) const = 0;
      /// typically set derived data for sub-domain (including upward pass)
      /// \note called for each sub-domain after its built
      /// \param[in] tree  tree to be built
      /// \param[in] it    thread index (if using OMP parallelism)
      /// \note @a SetSub() is called when the domain in question has been
      ///       linked, but not necessarily any other domain.
      virtual void SetSub(const OctalTree*tree, count_type it) const = 0;
      /// typically set derived data for top-tree (including upward pass)
      /// \param[in] tree  tree being build or updated
      /// \note called once after all sub-domains have been set
      virtual void SetTop(const OctalTree*tree) const = 0;
    };
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    ///
    /// interface supporting building of derived types
    ///
    struct TreeBuilderInterface 
      : public BuilderInterface
    {
      virtual ~TreeBuilderInterface() {}
      /// \param[in] nl  # internal leaves
      /// \param[in] ne  # external leaves
      /// \return        # bytes required in leaf data
      virtual size_t bytes_for_leaves(count_type, count_type)
	const { return 0; }
      /// use memory for leaves
      /// \param[in] mem  memory holding @a bytes_for_leaves() bytes
      /// \param[in] nl   # internal leaves
      /// \param[in] ne   # external leaves
      virtual void set_leaf_memory(char*, count_type, count_type)
	const {}
      /// \param[in] nc  # cells
      /// \return        # bytes required in leaf data
      virtual size_t bytes_for_cells(count_type) const { return 0; }
      /// use memory for leaves
      /// \param[in] mem  memory holding @a bytes_for_cells() bytes
      /// \param[in] nc   # cells
      virtual void set_cell_memory(char*, count_type) const {}
    };
#else
    typedef BuilderInterface TreeBuilderInterface;
#endif
  protected:
    /// protected ctor: allows parallel construction of derived
    OctalTree(const Initialiser*init, Data const&data,
	      TreeBuilderInterface const&builder) throw(exception);
    /// protected re-build: allows parallel re-build of derived
    void rebuild(Data const&data, TreeBuilderInterface const&builder)
      throw(exception);
    /// protected up-date of derived only
    /// \param[in] fresh   argument for BuilderInterface::Before()
    /// \param[in] builder BuilderInterface for building
    /// \note may be called from within parallel region
    void update_application_data(bool fresh, BuilderInterface const&builder)
      const
#ifdef OCTALTREE_USE_OPENMP
      ;
#else
    { builder.Serial(fresh,this); }
#endif
#ifdef OCTALTREE_USE_OPENMP
    /// \name miscellaneous
    //@{
    /// try to adjust # threads to # domains
    /// \param[in] func  name of calling function
    /// \note must not be called from within parallel region
    void AdjustNoThreads(const char*func=0) const
    {
      if(OMP::MaxNumThreads() != int(NDOM)) {
	if(OMP::IsParallel())
	  WDutils_ErrorN("%s(): mismatch (%d thread, %d domains) "
			 "inside parallel region\n",
			 func? func : "OctalTree::AdjustNoThreads",
			 OMP::MaxNumThreads(),int(NDOM));
	WDutils_WarningN("%s(): %d threads but %d tree domains: "
			 "trying to adjust # threads\n",
			 func? func : "OctalTree::AdjustNoThreads",
			 OMP::MaxNumThreads(),int(NDOM));
	OMP::SetMaxNumThreads(int(NDOM));
      }
    }
    /// asserts that # threads == # domains
    /// \param[in] func  name of calling function
    /// \note must be called from within parallel region
    void AssertNoThreads(const char*func=0) const
    {
      if(OMP::TeamSize() != int(NDOM))
	WDutils_ErrorN("%s(): %d threads but %d domains\n",
		       func? func : "OctalTree::AdjustNoThreads",
		       OMP::TeamSize(),int(NDOM));
    }
    //@}
#endif// OCTALTREE_USE_OPENMP
  public:
    /// \name general data of tree
    //@{
    /// # builds
    count_type const&n_build() const noexcept
    { return NBUILD; }
    /// was @a ascc set during tree building?
    bool const&avoided_single_child_cells() const noexcept
    { return ASCC; }
    /// tree depth
    depth_type const&depth() const noexcept
    { return _DP[0]; }
    /// root cell
    static Cell root() noexcept
    { return Cell(0); }
    /// root box
    cube const&root_box() const noexcept
    { return _XC[0]; }
    /// root centre
    point const&root_centre() const noexcept
    { return _XC[0].X; }
    /// root radius
    pos_type const&root_radius() const noexcept
    { return _XC[0].H; }
#ifdef OCTALTREE_USE_OPENMP
    /// tolerance in domain splitting
    count_type const&tol() const noexcept
    { return TOL; }
#endif// OCTALTREE_USE_OPENMP
    /// N_max
    depth_type const&n_max() const noexcept
    { return NMAX; }
    /// N_min
    depth_type const&n_min() const noexcept
    { return NMIN; }
    /// periodic boundary conditions
    /// \note returns a null pointer if we have no periodic boundary
    const PerBoundary*const&periodic() const noexcept
    { return PERB; }
    //@}
    /// \name functionality related to leaves
    //@{
    /// # internal leaves
    count_type n_leaf() const noexcept
    { return NLEAF; }
    /// # leaves actually allocated, multiple of LeafBlockSize, >= @a n_leaf()
    count_type leaf_capacity() const noexcept
    { return blocked(n_leaf()); }
    /// # 4-blocks of internal leaves
    count_type n_4blck() const noexcept
    { return (NLEAF+3)>>2; }
    /// begin of all internal leaves
    static Leaf begin_leaf() noexcept
    { return Leaf(0); }
    /// end of all internal leaves
    Leaf end_leaf() const noexcept
    { return Leaf(n_leaf()); }
    /// end of complete leaf blocks of Blocksize
    template<int Blocksize>
    Leaf end_complete_block() const noexcept
    {
      static_assert(Blocksize==2 || Blocksize==4 || Blocksize==8,
		    "end_complete_block<>: Blocksize must be 2,4, or 8");
      return Leaf(n_leaf() & ~(Blocksize-1));
    }
    /// position of leaf, 16-byte aligned
    template<bool EXT>
    point const&position(const TreeLeaf<EXT> l) const noexcept
    { return _XL[EXT][l.I]; }
    /// position of leaf, 16-byte aligned
    template<bool EXT>
    point16 const&position16(const TreeLeaf<EXT> l) const noexcept
    { return _XL[EXT][l.I]; }
    /// index of particle associated with leaf
    template<bool EXT>
    particle_key const&particle(const TreeLeaf<EXT> l) const noexcept
    { return _PL[EXT][l.I]; }
    /// pointer to all particle keys
    template<bool EXT>
    const particle_key*particle_keys() const noexcept
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    { return _PL[EXT]; }
#else
    { return _PL[EXT].data(); }
#endif
    /// do we have rungs loaded
    bool have_rung() const noexcept
    { return octtree::have_data(_RG[0]); }
    /// const access to rung for leaf
    template<bool EXT>
    float const&rung(const TreeLeaf<EXT> l) const noexcept
    { WDutilsTreeAssertE(octtree::have_data(_RG[EXT])); return _RG[EXT][l.I]; }
    /// const to rungs for internal leaves
    const float*cp_rung(const Leaf l) const noexcept
    { WDutilsTreeAssertE(have_rung()); return octtree::cp_data(_RG[0])+l.I; }
    /// do we have flags loaded
    bool have_flag() const noexcept
    { return octtree::have_data(_FL[0]); }
    /// const access to flag for leaf
    template<bool EXT>
    uint8_t const&flag(const TreeLeaf<EXT> l) const noexcept
    { WDutilsTreeAssertE(octtree::have_data(_FL[EXT])); return _FL[EXT][l.I]; }
    /// pointer to flags for internal leaf 
    const uint8_t*cp_flag(const Leaf l) const
    { WDutilsTreeAssertE(have_flag()); return octtree::cp_data(_FL[0])+l.I; }
    //
  private:
    template<int _B> typename std::enable_if<_B==2,const uint16_t*>::type
    _cp_block_flag(const Leaf l) const noexcept
    {
      WDutilsTreeAssertE((l.I&1)==0);
      return reinterpret_cast<const uint16_t*>(cp_flag(l));
    }
    template<int _B> typename std::enable_if<_B==4,const uint32_t*>::type
    _cp_block_flag(const Leaf l) const noexcept
    { 
      WDutilsTreeAssertE((l.I&3)==0);
      return reinterpret_cast<const uint32_t*>(cp_flag(l));
    }
    template<int _B> typename std::enable_if<_B==8,const uint64_t*>::type
    _cp_block_flag(const Leaf l) const noexcept
    {
      WDutilsTreeAssertE((l.I&7)==0);
      return reinterpret_cast<const uint64_t*>(cp_flag(l));
    }
  public:
    /// pointer to combined flag of block of leaves
    /// \note BlockSize must be 2,4,8 and @a l be BlockSize aligned
    template<int BlockSize>
    const typename octtree::block_flag_t<BlockSize>::type*
    cp_block_flag(const Leaf l) const noexcept
    {
      static_assert(BlockSize==2 || BlockSize==4 || BlockSize==8,
		    "OctTree::cp_block_flag(): block size must be 2, 4, or 8");
      return _cp_block_flag<BlockSize>(l);
    }
    /// combined flag of block of leaves
    /// \note BlockSize must be 2,4,8 and @a l be BlockSize aligned
    template<int BlockSize>
    typename octtree::block_flag_t<BlockSize>::type
    block_flag(const Leaf l) const noexcept
    { 
      static_assert(BlockSize==2 || BlockSize==4 || BlockSize==8,
		    "OctTree::block_flag(): block size must be 2, 4, or 8");
      return *cp_block_flag<BlockSize>(l);
    }
  private:
    //
    template<bool _A> typename std::enable_if< _A,bool>::type
    _is_active(const Leaf ) const noexcept { return 1; }
    template<bool _A> typename std::enable_if<!_A,bool>::type
    _is_active(const Leaf l) const noexcept { return flag(l)!=0; }
    template<bool _A, int BlockSize> typename std::enable_if< _A,bool>::type
    _block_has_active(const Leaf ) const noexcept { return 1; }
    template<bool _A, int BlockSize> typename std::enable_if<!_A,bool>::type
    _block_has_active(const Leaf l) const noexcept
    { return block_flag<BlockSize>(l)!=0; }
  public:
    /// is internal leaf active?
    template<bool All>
    bool is_active(const Leaf l) const noexcept
    { return _is_active<All>(l); }
    /// is any of a block of @a BlockSize internal leaves active?
    template<bool All, int BlockSize>
    bool block_has_active(const Leaf l) const noexcept
    { return _block_has_active<All, BlockSize>(l); }
    /// parent cell of internal leaf
    Cell parent(const Leaf l) const noexcept
    { return Cell(_PC[l.I]); }
    /// is @a l a valid internal leaf?
    bool is_valid(const Leaf l) const noexcept
    { return l.I < n_leaf(); }
    /// # external leaves
    count_type const&n_extleaf() const noexcept
    { return NLEAFEXT; }
    /// # external leaves allocated, multiple of LeafBlockSize,
    /// >= @a n_extleaf_()
    count_type extleaf_capacity() const noexcept
    { return blocked(n_cell()); }
    /// begin of external leaves
    static ExtLeaf begin_extleaf() noexcept
    { return ExtLeaf(0); }
    /// end of external leaves
    ExtLeaf end_extleaf() const noexcept
    { return ExtLeaf(n_extleaf()); }
  private:
    //
    template<bool _A> typename std::enable_if< _A,bool>::type
    _is_active(const ExtLeaf ) const noexcept { return 1; }
    template<bool _A> typename std::enable_if<!_A,bool>::type
    _is_active(const ExtLeaf l) const noexcept { return flag(l)!=0; }
  public:
    /// is leaf active?
    template<bool All>
    bool is_active(const ExtLeaf l) const noexcept
    { return _is_active<All>(l); }
    /// is @a l a valid leaf?
    bool is_valid(const ExtLeaf l) const noexcept
    { return l.I < n_extleaf(); }
    /// total # leaf ranges
    count_type n_range() const noexcept
    { return LeafRange::n_range(n_leaf()); }
    /// range of given index
    static LeafRange range_No(count_type r) noexcept
    { return LeafRange(r); }
    /// begin of leaf ranges
    static LeafRange begin_range() noexcept
    { return LeafRange(0); }
    /// end of leaf ranges
    LeafRange end_range() const noexcept
    { return LeafRange(n_range()); }
    /// first leaf in range
    static Leaf begin_leaf(const LeafRange r) noexcept
    { return Leaf(r.I<<LeafRange::Shft); }
    /// end of leaves in range
    Leaf end_leaf(const LeafRange r) const noexcept
    { return Leaf(std::min(n_leaf(),(r.I<<LeafRange::Shft)+LeafRange::Size)); }
    /// range of given leaf
    static LeafRange range_of(const Leaf l) noexcept
    { return LeafRange(l.I>>LeafRange::Shft); }
    //@}

    /// \name cell related functionality
    //@{
    /// # cells
    count_type n_cell() const noexcept
    { return NCELL; }
    /// is a cell valid?
    bool is_valid(const Cell c) const noexcept
    { return c.I < n_cell(); }
    /// is a cell the root cell?
    static bool is_root(const Cell c) noexcept
    { return c.I==0; }
    /// begin of cells
    static Cell begin_cell() noexcept
    { return Cell(0); }
    /// end of cells
    Cell end_cell() const noexcept
    { return Cell(n_cell()); }
    /// first cell in reversed order: last cell
    Cell rbegin_cell() const noexcept
    { return Cell(n_cell()-1); }
    /// end cell in reversed order: invalid Cell
    static Cell rend_cell() noexcept
    { return --(Cell(0)); }
    /// tree level of cell
    depth_type const&level(const Cell c) const noexcept
    { return _LE[c.I]; }
    /// tree depth of cell
    depth_type const&depth(const Cell c) const noexcept
    { return _DP[c.I]; }
    /// octant of cell in parent
    octant_type const&octant(const Cell c) const noexcept
    { return _OC[c.I]; }
    /// cell's cubic box, 16-byte aligned
    cube const&box(const Cell c) const noexcept
    { return _XC[c.I]; }
    /// cell's geometric centre (of cubic box)
    point const&centre(const Cell c) const noexcept
    { return box(c).X; }
    /// radius (half-side-length of box) of cell
    pos_type const&radius(const Cell c) const noexcept
    { return box(c).H; }
    /// number of leaf kids
    local_count const&n_leaf_kids(const Cell c) const noexcept
    { return _NL[c.I]; }
    /// total number of leaves
    count_type const&number(const Cell c) const noexcept
    { return _NM[c.I]; }
    /// # active leaf in cell
    count_type const&n_active(const Cell c) const noexcept
    { WDutilsTreeAssertE(have_flag()); return _NA[c.I]; }
    /// # active internal leaves
    count_type const&n_active() const noexcept
    { WDutilsTreeAssertE(have_flag()); return _NA[0]; }
  private:
    //
    template<bool _A> typename std::enable_if< _A,bool>::type
    _has_active(const Cell) const noexcept
    { return 1; }
    //
    template<bool _A> typename std::enable_if<!_A,bool>::type
    _has_active(const Cell c) const noexcept
    { return n_active(c)>0; }
  public:
    /// does cell have active leaves?
    template<bool All>
    bool has_active(const Cell c) const noexcept
    { return _has_active<All>(c); }
    /// number of daughter cells
    octant_type const&n_cell(const Cell c) const noexcept
    { return _NC[c.I]; }
    /// number of daughter cells
    octant_type const&n_daughters(const Cell c) const noexcept
    { return _NC[c.I]; }
    /// has a cell no daughter cells?
    bool is_final(const Cell c) const noexcept
    { return n_cell(c)==0; }
    /// first leaf in cell
    Leaf begin_leaf(const Cell c) const noexcept
    { return Leaf(_L0[c.I]); }
    /// end of leaf children (loose leaves) in cell
    Leaf end_leaf_kids(const Cell c) const noexcept
    { return Leaf(_L0[c.I]+_NL[c.I]); }
    /// end of all leaves in cell
    Leaf end_leaf_desc(const Cell c) const noexcept
    { return Leaf(_L0[c.I]+_NM[c.I]); }
    /// last of all leaves in cell
    Leaf last_leaf_desc(const Cell c) const noexcept
    { return Leaf(_L0[c.I]+_NM[c.I]-1); }
    /// first daughter cell
    /// \note if @a c has no daughter cells, this returns @a c itself
    Cell begin_cell(const Cell c) const noexcept
    { return Cell(_C0[c.I]); }
    /// end of daughter cells
    /// \note if @a c has no daughter cells, this returns @a c itself
    Cell end_cell(const Cell c) const noexcept
    { return Cell(_C0[c.I]+_NC[c.I]); }
    /// before first daughter cell
    Cell rend_cell(const Cell c) const noexcept
    { return --begin_cell(c); }
    /// last daughter cell
    Cell last_cell(const Cell c) const noexcept
    { return --end_cell(c); }
    /// parent cell
    Cell parent(const Cell c) const noexcept
    { return Cell(_PA[c.I]); }
    /// does a cell contain a leaf ?
    bool contains(const Cell c, const Leaf l) const noexcept
    { return begin_leaf(c) <= l && l < end_leaf_desc(c); }
    /// is either cell ancestor of the other?
    bool is_ancestor(const Cell a, const Cell b) const noexcept
    { return maxnorm(centre(a)-centre(b)) < max(radius(a),radius(b)); }
    /// find smallest cell containing a given position
    /// \param[in]  x   given position
    /// \return         smallest cell containing @a x
    /// \note If @a x is not in the root box, the root box is returned.
    Cell smallest_cell_containing(point const&x) const noexcept
    { return Cell(scc(x)); }
    /// find internal leaf with given key and at given position
    /// \param[in] i   particle key to search for
    /// \param[in] x   position of particle with key @a i
    /// \return internal leaf with key @a i, if found, otherwise invalid leaf
    /// \note If the position @a x is not that of the leaf with key @a i, we
    ///       cannot find the leaf with this method.
    Leaf find_particle(particle_key i, point const&x) const noexcept
    {
      Cell c = smallest_cell_containing(x);
      for(auto l=begin_leaf(c); l!=end_leaf_kids(c); ++l)
	if(i==particle(l)) return l;
      return Leaf::invalid();
    }
    //@}
#ifdef OCTALTREE_USE_OPENMP
    /// \name domain and top-tree related methods 
    /// \note useful for parallel upward pass
    //@{
    /// # domains
    domain_id n_dom() const noexcept
    { return NDOM; }
    /// domain pter given domain id
    cp_domain domain(domain_id d) const
    { WDutilsAssert(d<n_dom()); return DOM+d; }
    /// # top-tree cells
    count_type n_topcell() const noexcept
    { return NTOPC; }
    /// begin of top-tree cells
    static Cell begin_topcell() noexcept
    { return Cell(0); }
    /// end of top-tree cell
    Cell end_topcell() const noexcept
    { return Cell(n_topcell()); }
    /// first top-tree cell in reversed order: last top-tree cell
    Cell rbegin_topcell() const noexcept
    { return Cell(n_topcell()-1); }
    /// end top-tree cell in reversed order
    static Cell rend_topcell() noexcept
    { return --(Cell(0)); }
#ifdef OCTALTREE_HAVE_D0
    /// first domain contributing to cell
    domain_id const&first_domain(const Cell c) const noexcept
    { return _D0[c.I]; }
    /// domain of a given leaf
    domain_id domain_of(const Leaf l) const noexcept
    { return first_domain(parent(l)); }
#else
    /// domain of a given leaf
    domain_id domain_of(const Leaf l) const noexcept
    {
      count_type dom = (l.I * NDOM) / NLEAF;
      while(dom        && l.I <  DOM[dom].L0) --dom;
      while(dom < NDOM && l.I >= DOM[dom].LN) ++dom;
      return dom;
    }
    /// first domain contributing to cell
    domain_id first_domain(const Cell c) const noexcept
    { return domain_of(begin_leaf(c)); }
#endif
    /// # domains contributing to cell (>1 only for non-branch top-tree cells)
    domain_id n_domain(const Cell c) const noexcept
    { 
      if(c.I >= NTOPC || _C0[c.I] >= NTOPC || _NC[c.I]==0) return 1;
      return 1 + domain_of(last_leaf_desc(c)) - first_domain(c);
    }
    /// is a cell a top-tree cell?
    bool is_top(const Cell c) const noexcept
    { return c.I < NTOPC; }
    /// is a cell a non-branch top-tree cell?
    bool is_nonbranch_top(const Cell c) const noexcept
    { return is_top(c) && n_cell(c) && is_top(begin_cell(c)); }
    /// is a cell a non-branch sub-domain cell?
    bool is_sub(const Cell c) const noexcept
    { return c.I >= NTOPC; }
    /// is a cell a branch cell?
    bool is_branch(const Cell c) const noexcept
    { return is_top(c) && (n_cell(c)==0 || !is_top(begin_cell(c))); }
    //@}
#endif// OCTALTREE_USE_OPENMP

    /// \name dump tree data
    //@{
  private:
    template<int _D> typename std::enable_if<_D==2,std::ostream&>::type
    _LeafHeader(std::ostream&) const;
    template<int _D> typename std::enable_if<_D==3,std::ostream&>::type
    _LeafHeader(std::ostream&) const;
    template<int _D> typename std::enable_if<_D==2,std::ostream&>::type
    _CellHeader(std::ostream&) const;
    template<int _D> typename std::enable_if<_D==3,std::ostream&>::type
    _CellHeader(std::ostream&) const;
  protected:
    /// header for leaf dump
    /// \note may be extended in derived class
    virtual std::ostream&LeafHeader(std::ostream&out) const
    { return _LeafHeader<Dim>(out); }
    /// dump leaf data
    /// \note may be extended in derived class
    virtual std::ostream&Dump(const Leaf, std::ostream&) const;
    /// dump external leaf data
    /// \note may be extended in derived class
    virtual std::ostream&Dump(const ExtLeaf, std::ostream&) const;
    /// header for cell dump
    /// \note may be extended in derived class
    virtual std::ostream&CellHeader(std::ostream&out) const
    { return _CellHeader<Dim>(out); }
    /// dump cell data
    virtual std::ostream&Dump(const Cell, std::ostream&) const;
  public:
    /// dump all leaf data
    void DumpLeaves(std::ostream&s) const
    {
      LeafHeader(s) << '\n';
      for(auto l=begin_extleaf(); l!=end_extleaf(); ++l) Dump(l,s) << '\n';
      for(auto l=begin_leaf   (); l!=end_leaf   (); ++l) Dump(l,s) << '\n';
      s.flush();
    }
    /// dump all cell data
    void DumpCells(std::ostream&s) const
    {
      CellHeader(s) << '\n';
      for(auto c=begin_cell(); c!=end_cell(); ++c) Dump(c,s) << '\n';
      s.flush();
    }
    //@}
  protected:
    /// initialiser
    const Initialiser*init() const noexcept
    { return INIT; }
    /// index of smallest cell containing a point
    count_type scc(point const&x) const noexcept;
  };// class OctalTree<>
  //
  template<> struct traits< OctalTree<2,float> >
  { static const char*name() { return "OctalTree<2,float>"; } };
  template<> struct traits< OctalTree<3,float> >
  { static const char*name() { return "OctalTree<3,float>"; } };
  template<> struct traits< OctalTree<2,double> >
  { static const char*name() { return "OctalTree<2,double>"; } };
  template<> struct traits< OctalTree<3,double> >
  { static const char*name() { return "OctalTree<3,double>"; } };
  //
  template<> struct traits< OctalTree<2,float>::Dot >
  { static const char*name() { return "OctalTree<2,float>::Dot"; } };
  template<> struct traits< OctalTree<3,float>::Dot >
  { static const char*name() { return "OctalTree<3,float>::Dot"; } };
  template<> struct traits< OctalTree<2,double>::Dot >
  { static const char*name() { return "OctalTree<2,double>::Dot"; } };
  template<> struct traits< OctalTree<3,double>::Dot >
  { static const char*name() { return "OctalTree<3,double>::Dot"; } };
  //
  namespace octtree {
    template <template<int, typename> class temp_int_type>
    struct temp_int_type_conv
    {
      template <int _int, typename _type>
      temp_int_type_conv(const temp_int_type<_int,_type>&);
    };
  } // namespace octtree
  /// is a type @c Tree an instance of @c OctalTree<>?
  ///
  /// @c is_OctalTree<Tree>::value is true iff @c Tree is (derived from an) @c
  /// OctalTree<>
  template <class Tree> struct is_OctalTree
  {
    static const bool value =
      std::is_convertible<Tree,
			  octtree::temp_int_type_conv<OctalTree> >::value;
  };
  //
  namespace octtree {
    ///
    /// tags for various additional particle properties
    ///
    enum PropTag {
      mass_tag = 1 << 0,     ///< tag: particle mass
      size_tag = 1 << 1,     ///< tag: particle size^2
      xvec_tag = 1 << 2      ///< tag: SSE position for particles
    };
    ///
    /// Interface for loading additional particle properties
    ///
    /// \note can be implemented (derived) together with
    ///       @c OctalTree<>::Initialiser
    template<typename _PropType>
    struct ParticlePropertyInitialiser
    {
      /// dtor
      virtual ~ParticlePropertyInitialiser() {}
      /// set size^2 for all particles in a domain
      /// \param[in]  key    table: particle identifiers
      /// \param[out] sq     table: size^2 for particles associated with keys
      /// \param[in]  n      size of tables.
      virtual void InitSizeQ(const particle_key*key,
			     _PropType*sq, count_type n) const = 0;
      /// set mass for all particles in a domain
      /// \param[in]  key    table: particle identifiers
      /// \param[out] mass   table: mass for particles associated with keys
      /// \param[in]  n      size of tables.
      virtual void InitMass(const particle_key*key,
			    _PropType*mass, count_type n) const = 0;
    };
    ////////////////////////////////////////////////////////////////////////////

    ///
    /// interaction tree data for non-SSE usage
    ///
    template<int _Dim, typename _PosType, typename _PropType>
#ifdef __SSE__
    class InteractionTreeDataBase
#else
    class InteractionTreeData
#endif
      : protected static_tree_properties 
    {
      /// ensure what is implemented in octtree.cc
      static_assert(is_in_range<_Dim,2,3>::value,
		    "InteractionTree<>: only 2D and 3D is allowed");
      static_assert(std::is_floating_point<_PosType>::value,
		    "InteractionTree<>: position type must be floating point");
      static_assert(std::is_floating_point<_PropType>::value,
		    "InteractionTree<>: property type must be floating point");
    protected:
      typedef OctalTree<_Dim,_PosType> OTree;
      static const int Dim = _Dim;
      typedef _PropType prop_type;
#ifndef __SSE__
      friend class TreeImplementer<InteractionTreeData>;
    private:
#endif
      template<typename _Tp>
      using aligned_storage = typename OTree::template aligned_storage<_Tp>;
#ifdef INTERACTIONTREE_HOLDS_OWN_BLOCK
      raw_memory<octtree::DataAlignment> _BUF;
#endif
      aligned_storage<prop_type>         _MASS[2];  ///< any leaf: mass
      aligned_storage<prop_type>         _HSQ [2];  ///< any leaf: size^2
    protected:
      /// ctor
#ifdef __SSE__
      InteractionTreeDataBase()
#else
	InteractionTreeData()
#endif
      { 
#ifdef OCTALTREE_USE_OPENMP
	WDutilsAssertE(!OMP::IsParallel());
#endif
      }
      /// dtor
#ifdef __SSE__
      ~InteractionTreeDataBase()
#else
      ~InteractionTreeData()
#endif
      {
#ifdef OCTALTREE_USE_OPENMP
	WDutilsAssertE(!OMP::IsParallel());
#endif
      }
#ifndef __SSE__
# ifdef OCTALTREE_DATA_IN_ONE_BLOCK
      /// # bytes needed for nl internal and ne external leaves
      size_t bytes_needed(int load, count_type nl, count_type ne) const;
      /// set the memory
      void set_memory(int load, char*mem, count_type nl, count_type ne);
# endif
# ifndef INTERACTIONTREE_USES_BASE_BLOCK
      /// (re-)allocate data
      void allocate(const int, const count_type, const count_type);
# endif
#endif
    public:
      /// do we have size^2 allocated and loaded?
      bool have_Hsq() const noexcept
      { return have_data(_HSQ[0]); }
      /// const access to size^2 for internal or external leaf
      template<bool EXT>
	prop_type const&Hsq(const TreeLeaf<EXT> l) const noexcept
      { WDutilsTreeAssertE(have_data(_HSQ[EXT])); return _HSQ[EXT][l.I]; }
      /// const pter to size^2 for internal leaves
      const prop_type*cp_Hsq(const Leaf l) const noexcept
      { WDutilsTreeAssertE(have_Hsq()); return cp_data(_HSQ[0])+l.I; }
      /// non-const access to size squared for leaf
      template<bool EXT>
      prop_type&Hsq_lv(const TreeLeaf<EXT> l) const noexcept
      {
	WDutilsTreeAssertE(have_data(_HSQ[EXT]));
	return const_cast<prop_type&>(_HSQ[EXT][l.I]);
      }
      /// do we have masses allocated and loaded?
      bool have_mass() const noexcept
      { return have_data(_MASS[0]); }
      /// const access to mass for internal or external leaf
      template<bool EXT>
      prop_type const&mass(const TreeLeaf<EXT> l) const noexcept
      { WDutilsTreeAssertE(have_data(_MASS[EXT])); return _MASS[EXT][l.I]; }
      /// const pter to masses for internal leaves
      const prop_type*cp_mass(const Leaf l) const noexcept
      { WDutilsTreeAssertE(have_mass()); return cp_data(_MASS[0])+l.I; }
    };// class octtree::InteractionTreeDataBase<>
#ifdef __SSE__
    ///
    /// interaction tree data supporting SSE w/o anchoring
    ///
    template<int _Dim, typename _PosType, typename _PropType,
	     bool _Anchoring>
    class InteractionTreeData
      : public InteractionTreeDataBase<_Dim,_PosType,_PropType>
    {
      friend class TreeImplementer<InteractionTreeData>;
      typedef InteractionTreeDataBase<_Dim,_PosType,_PropType> IBase;
      typedef typename IBase::OTree OTree;
      typedef typename IBase::prop_type prop_type;
      using IBase::Dim;
#ifdef INTERACTIONTREE_HOLDS_OWN_BLOCK
      using IBase::_BUF;
#endif
      using IBase::_MASS;
      using IBase::_HSQ;
    public:
      using IBase::have_Hsq;
      using IBase::Hsq;
      using IBase::Hsq_lv;
      using IBase::cp_Hsq;
      using IBase::have_mass;
      using IBase::mass;
      using IBase::cp_mass;
      typedef _PosType sse_pos_type;
      typedef void Anchor;
    protected:
      template<typename _Tp>
      using aligned_storage = typename OTree::template aligned_storage<_Tp>;
# ifndef OCTALTREE_DATA_IN_ONE_BLOCK
      aligned_storage<sse_pos_type> _xvec;
# endif
      sse_pos_type                 *_XVEC[Dim];   ///< int. leaf: SSE position
    public:
# ifdef OCTALTREE_DATA_IN_ONE_BLOCK
      /// # bytes needed for nl internal and ne external leaves
      size_t bytes_needed(int load, count_type nl, count_type ne) const;
      /// set the memory
      void set_memory(int load, char*mem, count_type nl, count_type ne);
# endif
# ifndef INTERACTIONTREE_USES_BASE_BLOCK
      /// (re-)allocate data
      void allocate(const int, const count_type, const count_type);
# endif
      /// do we have SSE positions allocated and set?
      bool have_xvec() const noexcept
      { return _XVEC[0]!=0; }
      /// const pter to sse positions of leaves
      template<int I>
      const sse_pos_type*cp_xvec(const Leaf l) const noexcept
      { WDutilsTreeAssertE(have_xvec()); return _XVEC[I]+l.I; }
    };// class octtree::InteractionTreeData< Anchoring = false >
    ///
    /// interaction tree data supporting SSE with anchoring
    ///
    template<int _Dim>
    class InteractionTreeData<_Dim,double,float,true>
      : public InteractionTreeDataBase<_Dim,double,float>
    {
      typedef InteractionTreeDataBase<_Dim,double,float> IBase;
      typedef typename IBase::OTree OTree;
      typedef typename IBase::prop_type prop_type;
      typedef typename OTree::point point;
      template<typename _Tp>
      using storage         = typename OTree::template storage<_Tp>;
      template<typename _Tp>
      using aligned_storage = typename OTree::template aligned_storage<_Tp>;
      using IBase::Dim;
#ifdef INTERACTIONTREE_HOLDS_OWN_BLOCK
      using IBase::_BUF;
#endif
      using IBase::_MASS;
      using IBase::_HSQ;
    public:
      using IBase::have_Hsq;
      using IBase::Hsq;
      using IBase::Hsq_lv;
      using IBase::cp_Hsq;
      using IBase::have_mass;
      using IBase::mass;
      using IBase::cp_mass;
      typedef float sse_pos_type;
      struct Anchor {
	point Z;                                 ///< anchor position
	Leaf  B,E;                               ///< begin,end of leaves
      };
    protected:
# ifndef OCTALTREE_DATA_IN_ONE_BLOCK
      aligned_storage<sse_pos_type> _xvec;
# endif
      sse_pos_type                 *_XVEC[Dim];  ///< int. leaf: SSE position
      storage<Anchor*>              _ANC;        ///< int. leaf: pter to anchor
      std::vector<Anchor>           _ANCH;       ///< anchor data
    public:
# ifdef OCTALTREE_DATA_IN_ONE_BLOCK
      /// # bytes needed for nl internal and ne external leaves
      size_t bytes_needed(int load, count_type nl, count_type ne) const;
      /// set the memory
      void set_memory(int load, char*mem, count_type nl, count_type ne);
# endif
# ifndef INTERACTIONTREE_USES_BASE_BLOCK
      /// (re-)allocate data
      void allocate(const int, const count_type, const count_type);
# endif
      /// do we have SSE positions allocated and set?
      bool have_xvec() const noexcept
      { return _XVEC[0]!=0; }
      /// const pter to sse positions of leaves
      template<int I>
      const sse_pos_type*cp_xvec(const Leaf l) const noexcept
      { WDutilsTreeAssertE(have_xvec()); return _XVEC[I]+l.I; }
      /// do we have anchors allocated and set?
      bool have_anchor() const noexcept
      { return have_data(_ANC); }
      /// const pter to anchor for leaf
      const Anchor*anchor(const Leaf l) const noexcept
      { return _ANC[l.I]; }
      //
      friend class TreeImplementer<InteractionTreeData>;
    };// class octtree::InteractionTreeData< Anchoring = true >
#endif// __SSE__
  } // namespace octtree
  //
  template<> struct traits< octtree::InteractionTreeDataBase<2,float,float> >
  { static const char*name()
    { return "InteractionTreeDataBase<2,float,float>"; }
  };
  template<> struct traits< octtree::InteractionTreeDataBase<3,float,float> >
  { static const char*name()
    { return "InteractionTreeDataBase<3,float,float>"; }
  };
  template<> struct traits< octtree::InteractionTreeDataBase<2,double,double> >
  { static const char*name()
    { return "InteractionTreeDataBase<2,double,double>"; }
  };
  template<> struct traits< octtree::InteractionTreeDataBase<3,double,double> >
  { static const char*name()
    { return "InteractionTreeDataBase<3,double,double>"; }
  };
#ifdef __SSE__
  template<>
  struct traits<octtree::InteractionTreeData<2,double,float,true>::Anchor>
  { static const char*name() { return "Anchor"; } };
  template<>
  struct traits<octtree::InteractionTreeData<3,double,float,true>::Anchor>
  { static const char*name() { return "Anchor"; } };
  template<> struct traits< octtree::InteractionTreeData<2,float,float,0> >
  { static const char*name() { return "InteractionTreeData<2,float,float,0>"; }
  };
  template<> struct traits< octtree::InteractionTreeData<3,float,float,0> >
  { static const char*name() { return "InteractionTreeData<3,float,float,0>"; }
  };
  template<> struct traits< octtree::InteractionTreeData<2,double,float,0> >
  { static const char*name() { return "InteractionTreeData<2,double,float,0>"; }
  };
  template<> struct traits< octtree::InteractionTreeData<3,double,float,0> >
  { static const char*name() { return "InteractionTreeData<3,double,float,0>"; }
  };
  template<> struct traits< octtree::InteractionTreeData<2,double,float,1> >
  { static const char*name() { return "InteractionTreeData<2,double,float,1>"; }
  };
  template<> struct traits< octtree::InteractionTreeData<3,double,float,1> >
  { static const char*name() { return "InteractionTreeData<3,double,float,1>"; }
  };
#endif
  ///
  /// oct-tree with particle masses, sizes, and SSE positions.
  ///
  /// The floating-point type for the additional particle properties mass and
  /// size is a template parameter, which must not be of higher precision than
  /// that used for positions in base class @c OctalTree<>.
  ///
  /// SSE positions for the leaves differ from their ordinary positions in the
  /// follwoing ways: (1) they are organised in component arrays (one for X, one
  /// for Y, and one for Z); (2) they are single precision, unless
  /// _Anchoring==false and _PropType==double; and (3) if anchoring is enabled
  /// and _PosType==double, they are relative to an anchor. SSE positions are
  /// not provided for external leaves.
  ///
  /// Anchoring is enabled if _Anchoring==true, _PropType==float, and
  /// _PosType==double. If anchoring is enabled and particle sizes are loaded,
  /// anchors and sse positions will be set. This is done such that the
  /// numerical precision of each particle's size relative to its sse position
  /// is not below a tolerance. Leaves that cannot be anchored will have null
  /// pter for their anchor and sse positions relative to their parent cell.
  ///
  /// Neighbour search and neighbour processing using anchored and unanchored
  /// @c InteractionTree is supported by the methods of <neighbour.h>
  ///
  template<int _Dim, typename _PosType, typename _PropType = float
#ifdef __SSE__
	   , bool _Anchoring = 
	   std::is_same<_PosType,double>::value &&
	   std::is_same<_PropType,float>::value
#endif
	   >
  class InteractionTree : 
    public octtree::InteractionTreeData<_Dim,_PosType,_PropType
#ifdef __SSE__
					,_Anchoring
#endif
					>,
    public OctalTree<_Dim,_PosType>
  {
  public:
    typedef OctalTree<_Dim,_PosType> OTree;
    using OTree::Dim;
    using OTree::Nsub;
    using OTree::Noff;
    using prop_type = _PropType;
#ifdef __SSE__
    static const bool Anchoring = _Anchoring; ///< anchoring enabled?
#endif// __SSE__
    typedef typename OTree::pos_type pos_type;
    typedef typename OTree::point point;
    using DataBase = octtree::InteractionTreeData<Dim,pos_type,prop_type
#ifdef __SSE__
					       ,_Anchoring
#endif
					       >;
#ifdef __SSE__
    typedef typename DataBase::Anchor Anchor;
    typedef typename DataBase::sse_pos_type sse_pos_type;
#endif// __SSE__
    typedef octtree::Leaf Leaf;
    typedef octtree::Cell Cell;
    typedef octtree::ExtLeaf ExtLeaf;
    typedef octtree::count_type count_type;
  protected:
    /// \name data dumping
    //@{
    /// header for leaf dump
    virtual std::ostream&LeafHeader(std::ostream&) const;
    /// dump internal leaf
    virtual std::ostream&Dump(const Leaf, std::ostream&) const;
    /// dump external leaf
    virtual std::ostream&Dump(const ExtLeaf, std::ostream&) const;
    /// enable Dump(Cell)
    using OTree::Dump;
    //@}
  public:
    /// particle property initialiser
    using PropInitialiser = octtree::ParticlePropertyInitialiser<prop_type>;
    /// Interface for obtaining particle sizes
    struct Initialiser : OTree::Initialiser, PropInitialiser {};
    /// data needed in ctor
    struct Data : OTree::Data
    {
      int   LOAD;               ///< bits from octtree::PropTag: props to load
      float DEL;                ///< relative precision of anchored positions.
      Data() : LOAD(0), DEL(0.0001f)  {}
      Data(Data&&) = default;
      Data(Data const&) = default;
      Data&operator=(Data&&) = default;
      Data&operator=(Data const&) = default;
    };
  protected:
    friend struct Builder;
    /// builder for additional data
    struct Builder : OTree::TreeBuilderInterface
    {
      //
      Builder(DataBase*, const Initialiser*, Data const&);
#ifdef INTERACTIONTREE_USES_BASE_BLOCK
      size_t bytes_for_leaves(count_type nl, count_type ne) const
      { return DATA->bytes_needed(LOAD,nl,ne); }
      void set_leaf_memory(char*buf, count_type nl, count_type ne) const
      { DATA->set_memory(LOAD,buf,nl,ne); }
#endif
      virtual void Before(bool, count_type, count_type, count_type) const;
      virtual void SetSub(const OTree*, count_type) const WD_HOT;
      virtual void SetTop(const OTree*) const WD_HOT;
      ///
    private:
      const PropInitialiser*INIT;
      DataBase*DATA;
      int      LOAD;
      float    DEL;
#ifdef __SSE__
      float   *QI;
    public:
      ~Builder() { if(QI) WDutils_DEL_A(QI); QI=0; }
#endif
    };
  public:
    /// initialiser
    const Initialiser*init() const noexcept
    { return static_cast<const Initialiser*>(OTree::init()); }
    /// ctor
    /// \param[in] _init   initialiser: set size^2, mass, flag
    /// \param[in] _data   data specifying details.
    ///
    /// \note If @a _data.DEL>0, the SSE positions are set to be relative to an
    ///       anchor. The anchors are cell centres chosen such that the
    ///       numerical precision is at least @a _data.DEL times the particle
    ///       size H. Thus @a _data.DEL must be (much) greater than 1.2e-7
    ///       (single precision) but (much) less then 1. The default is 1.e-4.
    ///
    InteractionTree(const Initialiser*_init, Data const&_data) throw(exception)
      : OTree(_init,_data, Builder(this,_init,_data)) {}
    /// re-build tree
    void rebuild(Data const&_data) throw(exception)
    { OTree::rebuild(_data, Builder(this,init(),_data)); }
    /// dtor
    virtual~InteractionTree() {}
  protected:
    /// protected ctor: allows parallel construction of derived
    InteractionTree(const Initialiser*ini, Data const&d,
		    Builder const&builder) throw(exception)
      : OTree(ini,d,builder) {}
    /// protected re-build: allows parallel re-build of derived
    void rebuild(Data const&d, Builder const&builder) throw(exception)
    { OTree::rebuild(d,builder); }
    ///
    template<typename, typename> friend class TreeWalker;
    friend struct TreeImplementer<DataBase>;
  };// class InteractionTree<>
  //
  template<> struct traits< InteractionTree<2,float,float> >
  { static const char*name() { return "InteractionTree<2,float,float>"; } };
  template<> struct traits< InteractionTree<3,float,float> >
  { static const char*name() { return "InteractionTree<3,float,float>"; } };
  template<> struct traits< InteractionTree<2,double,double> >
  { static const char*name() { return "InteractionTree<2,double,double>"; } };
  template<> struct traits< InteractionTree<3,double,double> >
  { static const char*name() { return "InteractionTree<3,double,double>"; } };
  template<> struct traits< InteractionTree<2,double,float,0> >
  { static const char*name() { return "InteractionTree<2,double,float,0>"; } };
  template<> struct traits< InteractionTree<3,double,float,0> >
  { static const char*name() { return "InteractionTree<3,double,float,0>"; } };
  template<> struct traits< InteractionTree<2,double,float,1> >
  { static const char*name() { return "InteractionTree<2,double,float,1>"; } };
  template<> struct traits< InteractionTree<3,double,float,1> >
  { static const char*name() { return "InteractionTree<3,double,float,1>"; } };
  //
  namespace octtree {
    template <template<int, typename, typename, bool> 
	      class temp_int_type_type_bool >
    struct temp_int_type_type_bool_conv
    {
      template <int _int, typename _type0, typename _type1, bool _bool >
      temp_int_type_type_bool_conv
      (const temp_int_type_type_bool<_int,_type0,_type1,_bool>&);
    };
    template <template<int, typename, typename> 
	      class temp_int_type_type >
    struct temp_int_type_type_conv
    {
      template <int _int, typename _type0, typename _type1 >
      temp_int_type_type_conv
      (const temp_int_type_type<_int,_type0,_type1>&);
    };
  } // namespace octtree
  /// is a type @c Tree an instance of @c InteractionTree<>
  ///
  /// @c is_InteractionTree<Tree>::value is true iff @c Tree is (derived from
  /// an) @c InteractionTree<>
  template <class Tree> struct is_InteractionTree
  {
    static const bool value =
#ifdef __SSE__
      std::is_convertible<
        Tree,
        octtree::temp_int_type_type_bool_conv<InteractionTree>
      > ::value;
#else
      std::is_convertible<
	Tree,
	octtree::temp_int_type_type_conv<InteractionTree>
      > ::value;
#endif
  };
  //
  namespace octtree {
    template<typename _Tree, bool = is_InteractionTree<_Tree>::value>
    struct itree_types;
    template<typename _Tree> struct itree_types<_Tree,false> {
      using prop_type = void;
#ifdef __SSE__
      using Anchor = void;
      using sse_pos_type = void;
#endif
    };
    template<typename _Tree> struct itree_types<_Tree,true> {
      using prop_type = typename _Tree::prop_type;
#ifdef __SSE__
      using Anchor = typename _Tree::Anchor;
      using sse_pos_type = typename _Tree::sse_pos_type;
#endif
    };
  } // namespace octtree
  ///
  /// access to tree
  ///
  /// rather than deriving from a Tree, we take a const Tree* member, but allow
  /// derived classes to construct and re-build synchronously with the tree.
  ///
  template<typename _Tree>
  class TreeWalker<_Tree,
		   typename enable_if<is_OctalTree<_Tree>::value>::type> : 
    public octtree::static_tree_properties,
    public octtree::itree_types<_Tree>
  {
  public:
    typedef _Tree Tree;
    static const bool IsITree = is_InteractionTree<Tree>::value;
    typedef octtree::itree_types<Tree> i_types;
    static const int Dim  = Tree::Dim;
    static const int Nsub = Tree::Nsub;
    static const int Noff = Tree::Noff;
    using typename i_types::prop_type;
#ifdef __SSE__
    using typename i_types::Anchor;
    using typename i_types::sse_pos_type;
#endif
    //
    template<bool _C, typename _T=void>
    using enable_if_type = typename std::enable_if<_C,_T>::type;
    //
    typedef typename Tree::BuilderInterface BuilderInterface;
    typedef typename Tree::pos_type pos_type;
    typedef typename Tree::cube cube;
    typedef typename Tree::cube16 cube16;
    typedef typename Tree::point point;
    typedef typename Tree::point16 point16;
    typedef typename Tree::PerBoundary PerBoundary;
    typedef typename Tree::cp_domain cp_domain;
    /// \name data
    //@{
  protected:
    const Tree*const TREE;               ///< 'our' tree
    //@}
  public:
    ///
    /// \name construction and re-building
    //@{
    /// ctor from a built tree: import the tree
    explicit TreeWalker(const Tree*t) noexcept
      : TREE(t) {}
    /// copy ctor
    TreeWalker(TreeWalker const&) = default;
    /// move ctor
    TreeWalker(TreeWalker &&) = default;
    /// copy operator
    TreeWalker&operator=(TreeWalker const&) = default;
    /// move operator
    TreeWalker&operator=(TreeWalker &&) = default;
    /// dtor
    virtual~TreeWalker() {}
  protected:
    /// protected update of application data
    /// \param[in] fresh    argument for builder::Before()
    /// \param[in] builder  Builder for building
    /// \note may be called from within parallel region
    void update(bool fresh, BuilderInterface const&builder)
    { TREE->update_application_data(fresh,builder); }
  private:
    /// no default ctor
    TreeWalker() = delete;
    //@}
#ifdef OCTALTREE_HAVE_DATA
  protected:
    /// \name access to additional tree data for derived classes
    //@{
    /// first available free datum, to be overridden in derived
    static const int FREE = Tree::FREE;
    /// end for free data: a derived class should check that FREE <= NDAT
    static const int NDAT = Tree::NDAT;
    /// set any_domain::DATA[i] to given pointer
    /// \param[in] p  pointer to set data to
    template<bool EXT, int I> void SetData(void*p) noexcept
    { const_cast<Tree*>(TREE)->template SetData<EXT,I>(p); }
    /// return pointer set by SetData<I>
    template<bool EXT, int I> void*GetData() const noexcept
    { return TREE->template GetData<EXT,I>(); }
    /// const access to DATA[I][i]
    template<bool EXT, int I, typename T>
    T const&dat(count_type n) const noexcept
    { return TREE-> template dat<EXT,I,T>(n); }
    /// non-const access to DATA[I][i]
    template<bool EXT, int I, typename T>
    T&dat_lv(count_type n) const noexcept
    { return TREE-> template dat_lv<EXT,I,T>(n); }
    /// const access to user vector datum
    template<bool EXT, int I, int J, typename T> 
    T const&vec(count_type n) const noexcept
    { return TREE-> template vec<EXT,I,J,T>(n); }
    /// non-const access to user vector datum
    template<bool EXT, int I, int J, typename T> 
    T&vec_lv(count_type n) const noexcept
    { return TREE-> template vec_lv<EXT,I,J,T>(n); }
    /// pointer to user datum
    template<bool EXT, int I, typename T>
    T*p_dat(count_type n) const noexcept
    { return TREE-> template p_dat<EXT,I,T>(n); }
    /// constant pointer to user datum
    template<bool EXT, int I, typename T>
    const T*cp_dat(count_type n) const noexcept
    { return TREE-> template cp_dat<EXT,I,T>(n); }
    /// non-const pointer to vector component
    template<bool EXT, int I, int J, typename T>
    T*p_vec(count_type n) const noexcept
    { return TREE -> template p_vec<EXT,I,J,T>(n); }
    /// const pointer to vector component from free methods
    template<bool EXT, int I, int J, typename T>
    T*cp_vec(count_type n) const noexcept
    { return TREE -> template cp_vec<EXT,I,J,T>(n); }
    //@}
#endif // OCTALTREE_HAVE_DATA
  public:
    /// \name general tree related data
    //@{
    /// tree
    const Tree*tree() const noexcept
    { return TREE; }
    /// was @a ascc set during tree building?
    bool avoided_single_child_cells() const noexcept
    { return TREE->avoided_single_child_cells(); }
#ifdef OCTALTREE_USE_OPENMP
    /// tolerance in domain splitting
    count_type const&tol() const noexcept
    { return TREE->tol(); }
#endif// OCTALTREE_USE_OPENMP
    /// N_max
    depth_type n_max() const noexcept
    { return TREE->n_max(); }
    /// N_min
    depth_type n_min() const noexcept
    { return TREE->n_min(); }
    /// # builds
    count_type n_build() const noexcept
    { return TREE->n_build(); }
    /// tree depth
    depth_type depth() const noexcept
    { return TREE->depth(); }
    /// root cell
    static Cell root() noexcept
    { return Tree::root(); }
    /// root box
    cube root_box() const noexcept
    { return TREE->root_box(); }
    /// root centre
    point root_centre() const noexcept
    { return TREE->root_centre(); }
    /// root radius
    pos_type root_radius() const noexcept
    { return TREE->root_radius(); }
    /// periodic boundary conditions
    /// \note returns a null pointer if we have no periodic boundary
    const PerBoundary*const&periodic() const noexcept
    { return TREE->periodic(); }
    //@}

    /// \name functionality related to leaves
    //@{
    /// # internal leaves
    count_type n_leaf() const noexcept
    { return TREE->n_leaf(); }
    /// # leaves actually allocated, multiple of LeafBlockSize, >= @a n_leaf()
    count_type leaf_capacity() const noexcept
    { return TREE->leaf_capacity(); }
    /// # 4-blocks of internal leaves
    count_type n_4blck() const noexcept
    { return TREE->n_4blck(); }
    /// leaf of given global index @a n in @a [0,n_leaf()[
    static Leaf leaf_no(count_type n) noexcept
    { return Leaf(n); }
    /// begin of all internal leaves
    static Leaf begin_leaf() noexcept
    { return Tree::begin_leaf(); }
    /// end of all internal leaves
    Leaf end_leaf() const noexcept
    { return TREE->end_leaf(); }
    /// end of complete leaf blocks of N
    template<int Blocksize>
    Leaf end_complete_block() const noexcept
    { return TREE->template end_complete_block<Blocksize>(); }
    /// position of internal leaf, 16-byte aligned
    template<bool EXT>
    point const&position(const TreeLeaf<EXT> l) const noexcept
    { return TREE->position(l); }
    /// position of internal leaf, 16-byte aligned
    template<bool EXT>
    point16 const&position16(const TreeLeaf<EXT> l) const noexcept
    { return TREE->position16(l); }
    /// index of particle associated with internal leaf
    template<bool EXT>
    particle_key const&particle(const TreeLeaf<EXT> l) const noexcept
    { return TREE->particle(l); }
    /// do we have rungs loaded
    bool have_rung() const noexcept
    { return TREE->have_rung(); }
    /// const access to rung for leaf
    template<bool EXT>
    float const&rung(const TreeLeaf<EXT> l) const noexcept
    { return TREE->rung(l); }
    /// const to rungs for internal leaves
    const float*cp_rung(const Leaf l) const noexcept
    { return TREE->cp_rung(l); }
    /// do we have flags loaded
    bool have_flag() const noexcept
    { return TREE->have_flag(); }
    /// const access to flag for leaf
    template<bool EXT>
    uint8_t const&flag(const TreeLeaf<EXT> l) const noexcept
    { return TREE->flag(l); }
    /// pointer to flags for internal leaf 
    const uint8_t*cp_flag(const Leaf l) const
    { return TREE->cp_flag(l); }
    /// find internal leaf with given key and at given position
    /// \param[in] i   particle key to search for
    /// \param[in] x   position of particle with key @a i
    /// \return internal leaf with key @a i, if found, otherwise invalid leaf
    /// \note If the position @a x is not that of the leaf with key @a i, we
    ///       cannot find the leaf with this method.
    Leaf find_particle(particle_key i, point const&x) const noexcept
    { return TREE->find_particle(i,x); }
    /// pointer to combined flag of block of leaves
    /// \note BlockSize must be 2,4,8 and @a l be BlockSize aligned
    template<int BlockSize>
    const typename octtree::block_flag_t<BlockSize>::type*
    cp_block_flag(const Leaf l) const noexcept
    { return TREE->template cp_block_flag<BlockSize>(l); }
    /// combined flag of block of leaves
    /// \note BlockSize must be 2,4,8 and @a l be BlockSize aligned
    template<int BlockSize>
    typename octtree::block_flag_t<BlockSize>::type
    block_flag(const Leaf l) const noexcept
    { return TREE->template block_flag<BlockSize>(l); }
    /// is leaf active?
    template<bool All, bool EXT>
    bool is_active(const TreeLeaf<EXT> l) const noexcept
    { return TREE->template is_active<All>(l); }
    /// is any of a block of 4 internal leaves active?
    template<bool All, int BlockSize>
    bool block_has_active(const Leaf l) const noexcept
    { return TREE->template block_has_active<All,BlockSize>(l); }
    /// parent cell of internal leaf
    Cell parent(const Leaf l) const noexcept
    { return TREE->parent(l); }
    /// is @a l a valid leaf?
    template<bool EXT>
    bool is_valid(const TreeLeaf<EXT> l) const noexcept
    { return TREE->is_valid(l); }
    /// # external leaves
    count_type const&n_extleaf() const noexcept
    { return TREE->n_extleaf(); }
    /// # external leaves allocated, multiple of LeafBlockSize,
    /// >= @a n_extleaves()
    count_type extleaf_capacity() const noexcept
    { return TREE->extleaf_capacity(); }
    /// begin of external leaves
    static ExtLeaf begin_extleaf() noexcept
    { return Tree::begin_extleaf(); }
    /// end of external leaves
    ExtLeaf end_extleaf() const noexcept
    { return TREE->end_extleaf(); }
    /// total # leaf ranges
    count_type n_range() const noexcept
    { return TREE->n_range(); }
    /// range of given index
    static LeafRange range_no(count_type r) noexcept
    { return LeafRange(r); }
    /// begin of leaf ranges
    static LeafRange begin_range() noexcept
    { return Tree::begin_range(); }
    /// end of leaf ranges
    LeafRange end_range() const noexcept
    { return TREE->end_range(); }
    /// first leaf in range
    static Leaf begin_leaf(const LeafRange r) noexcept
    { return Tree::begin_leaf(r); }
    /// end of leaves in range
    Leaf end_leaf(const LeafRange r) const noexcept
    { return TREE->end_leaf(r); }
    /// range of given leaf
    static LeafRange range_of(const Leaf l) noexcept
    { return Tree::range_of(l); }
    //
  private:
    template<bool i> enable_if_type< i,bool>
    _have_mass() const noexcept { return TREE->have_mass(); }
    template<bool i> enable_if_type<!i,bool>
    _have_mass() const noexcept { return 0; }
    //
    template<bool i, bool EXT> enable_if_type< i,prop_type> const&
    _mass(const TreeLeaf<EXT> l) const noexcept { return TREE->mass(l); }
    template<bool i, bool EXT> enable_if_type<!i> 
    _mass(TreeLeaf<EXT>) const noexcept {}
    //
    template<bool i, bool EXT> enable_if_type< i,prop_type> const*
    _cp_mass(const TreeLeaf<EXT> l) const noexcept { return TREE->cp_mass(l); }
    template<bool i, bool EXT> enable_if_type<!i> 
    _cp_mass(TreeLeaf<EXT>) const noexcept {}
    //
    template<bool i> enable_if_type< i,bool>
    _have_Hsq() const noexcept { return TREE->have_Hsq(); }
    template<bool i> enable_if_type<!i,bool>
    _have_Hsq() const noexcept { return 0; }
    //
    template<bool i, bool EXT> enable_if_type< i,prop_type> const&
    _Hsq(const TreeLeaf<EXT> l) const noexcept { return TREE->Hsq(l); }
    template<bool i, bool EXT> enable_if_type<!i> 
    _Hsq(TreeLeaf<EXT>) const noexcept {}
    //
    template<bool i, bool EXT> enable_if_type< i,prop_type> &
    _Hsq_lv(const TreeLeaf<EXT> l) const noexcept { return TREE->Hsq_lv(l); }
    template<bool i, bool EXT> enable_if_type<!i>
    _Hsq_lv(TreeLeaf<EXT>) const noexcept {}
    //
    template<bool i, bool EXT> enable_if_type< i,prop_type> const*
    _cp_Hsq(const TreeLeaf<EXT> l) const noexcept { return TREE->cp_Hsq(l); }
    template<bool i, bool EXT> enable_if_type<!i> 
    _cp_Hsq(TreeLeaf<EXT>) const noexcept {}
#ifdef __SSE__
    //
    template<bool i> enable_if_type< i,bool>
    _have_xvec() const noexcept { return TREE->have_xvec(); }
    template<bool i> enable_if_type<!i,bool>
    _have_xvec() const noexcept { return 0; }
    //
    template<bool i, int I> enable_if_type< i, sse_pos_type> const*
    _cp_xvec(const Leaf l) const noexcept
    { return TREE-> template cp_xvec<I>(l); }
    template<bool i, int I> enable_if_type<!i>
    _cp_xvec(Leaf) const noexcept {}
    //
    template<bool i> enable_if_type< i,bool>
    _have_anchor() const noexcept { return TREE->have_anchor(); }
    template<bool i> enable_if_type<!i,bool>
    _have_anchor() const noexcept { return 0; }
    //
    template<bool i> enable_if_type< i, Anchor> const*
    _anchor(const Leaf l) const noexcept { return TREE->anchor(l); }
    template<bool i> enable_if_type<!i>
    _anchor(Leaf) const noexcept {}
#endif
  public:
    // NOTE: these dummy methods are not implemented, but provide a public way
    //       to get the appropriate return type. They are necessary because a
    //       bug in GCC 4.7.0 prevents the use of private methods within the
    //       decltype specifier of the return type of public methods.
    template<bool i> static enable_if_type< i,prop_type> const&_prop();
    template<bool i> static enable_if_type<!i> _prop();
    //
    template<bool i> static enable_if_type< i,prop_type>&_prop_lv();
    template<bool i> static enable_if_type<!i> _prop_lv();
    //
    template<bool i> static enable_if_type< i,prop_type> const*_cp_prop();
    template<bool i> static enable_if_type<!i> _cp_prop();
    //
#ifdef __SSE__
    template<bool i> static enable_if_type< i,sse_pos_type> const*_cp_sse();
    template<bool i> static enable_if_type<!i> _cp_sse();
    //
    template<bool i> static enable_if_type< i,Anchor> const*_cp_anch();
    template<bool i> static enable_if_type<!i> _cp_anch();
#endif
    /// do we have masses?
    bool have_mass() const noexcept
    { return _have_mass<IsITree>(); }
    /// const access to leaf mass
    template<bool EXT>
    auto mass(const TreeLeaf<EXT> l) const noexcept
      -> decltype(_prop<IsITree>()) { return _mass<IsITree>(l); }
    /// const pointer to leaf mass
    template<bool EXT>
    auto cp_mass(const TreeLeaf<EXT> l) const noexcept
      -> decltype(_cp_prop<IsITree>()) { return Tree::_cp_mass<IsITree>(l); }
    /// do we have sizes?
    bool have_Hsq() const noexcept
    { return _have_Hsq<IsITree>(); }
    /// const access to leaf size^2
    template<bool EXT>
    auto Hsq(const TreeLeaf<EXT> l) const noexcept
      -> decltype(_prop<IsITree>()) { return _Hsq<IsITree>(l); }
    /// non-const access to leaf size^2
    template<bool EXT>
    auto Hsq_lv(TreeLeaf<EXT> const&l) const noexcept
      -> decltype(_prop_lv<IsITree>()) { return _Hsq_lv<IsITree>(l); }
    /// const pointer to leaf size^2
    template<bool EXT>
    auto cp_Hsq(const TreeLeaf<EXT> l) const noexcept
      -> decltype(_cp_prop<IsITree>()) { return _cp_Hsq<IsITree>(l); }
#ifdef __SSE__
    /// do we have sse positions?
    bool have_xvec() const noexcept
    { return _have_xvec<IsITree>(); }
    /// const pter to sse positions of leafs
    template<int I>
    auto cp_xvec(const Leaf l) const noexcept
      -> decltype(_cp_sse<IsITree>()) { return _cp_xvec<IsITree,I>(l); }
    /// do we have anchors
    bool have_anchor() const noexcept
    { return _have_anchor<IsITree>(); }
    /// const access to leaf anchor
    auto anchor(const Leaf l) const noexcept
      -> decltype(_cp_anch<IsITree>()) { return _anchor<IsITree>(l); }
#endif
    //@}

    /// \name functionality related to cells
    //@{
    /// # cells
    count_type n_cell() const noexcept
    { return TREE->n_cell(); }
    /// is a cell valid?
    bool is_valid(const Cell c) const noexcept
    { return TREE->is_valid(c); }
    /// is a cell the root cell?
    static bool is_root(const Cell c) noexcept
    { return Tree::is_root(c); }
    /// begin of cells
    static Cell begin_cell() noexcept
    { return Tree::begin_cell(); }
    /// end of cells
    Cell end_cell() const noexcept
    { return TREE->end_cell(); }
    /// first cell in reversed order: last cell
    Cell rbegin_cell() const noexcept
    { return TREE->rbegin_cell(); }
    /// end cell in reversed order: invalid Cell
    static Cell rend_cell() noexcept
    { return Tree::rend_cell(); }
    /// tree level of cell
    depth_type const&level(const Cell c) const noexcept
    { return TREE->level(c); }
    /// tree depth of cell
    depth_type const&depth(const Cell c) const noexcept
    { return TREE->depth(c); }
    /// octant of cell in parent
    octant_type const&octant(const Cell c) const noexcept
    { return TREE->octant(c); }
    /// cell's cubic box, 16-byte aligned
    cube const&box(const Cell c) const noexcept
    { return TREE->box(c); }
    /// cell's geometric centre (of cubic box)
    point const&centre(const Cell c) const noexcept
    { return TREE->centre(c); }
    /// radius (half-side-length of box) of cell
    pos_type const&radius(const Cell c) const noexcept
    { return TREE->radius(c); }
    /// number of leaf kids
    local_count const&n_leaf_kids(const Cell c) const noexcept
    { return TREE->n_leaf_kids(c); }
    /// total number of leaves
    count_type const&number(const Cell c) const noexcept
    { return TREE->number(c); }
    /// # active leaf in cell
    count_type const&n_active(const Cell c) const noexcept
    { return TREE->n_active(c); }
    /// total # active internal leaves
    count_type const&n_active() const noexcept
    { return TREE->n_active(); }
    /// does cell have active leaves?
    template<bool All>
    bool has_active(const Cell c) const noexcept
    { return TREE->template has_active<All>(c); }
    /// number of daughter cells
    octant_type const&n_cell(const Cell c) const noexcept
    { return TREE->n_cell(c); }
    /// number of daughter cells
    octant_type const&n_daughters(const Cell c) const noexcept
    { return TREE->n_daughters(c); }
    /// has a cell no daughter cells?
    bool is_final(const Cell c) const noexcept
    { return TREE->is_final(c); }
    /// first leaf in cell
    Leaf begin_leaf(const Cell c) const noexcept
    { return TREE->begin_leaf(c); }
    /// end of leaf children (loose leaves) in cell
    Leaf end_leaf_kids(const Cell c) const noexcept
    { return TREE->end_leaf_kids(c); }
    /// end of all leaves in cell
    Leaf end_leaf_desc(const Cell c) const noexcept
    { return TREE->end_leaf_desc(c); }
    /// last of all leaves in cell
    Leaf last_leaf_desc(const Cell c) const noexcept
    { return TREE->last_leaf_desc(c); }
    /// first daughter cell
    /// \note if @a c has no daughter cells, this returns @a c itself
    Cell begin_cell(const Cell c) const noexcept
    { return TREE->begin_cell(c); }
    /// end of daughter cells
    /// \note if @a c has no daughter cells, this returns @a c itself
    Cell end_cell(const Cell c) const noexcept
    { return TREE->end_cell(c); }
    /// before first daughter cell
    Cell rend_cell(const Cell c) const noexcept
    { return TREE->rend_cell(c); }
    /// last daughter cell
    Cell last_cell(const Cell c) const noexcept
    { return TREE->last_cell(c); }
    /// parent cell
    Cell parent(const Cell c) const noexcept
    { return TREE->parent(c); }
    /// does a cell contain a leaf ?
    bool contains(const Cell c, const Leaf l) const noexcept
    { return TREE->contains(c,l); }
    /// is either cell ancestor of the other?
    bool is_ancestor(const Cell a, const Cell b) const noexcept
    { return TREE->is_ancestor(a,b); }
    /// find smallest cell containing a given position
    /// \param[in]  x   given position
    /// \return         smallest cell containing @a x
    /// \note If @a x is not in the root box, the root box is returned.
    Cell smallest_cell_containing(point const&x) const noexcept
    { return TREE->smallest_cell_containing(x); }
    //@}

#ifdef OCTALTREE_USE_OPENMP
    /// \name domain and top-tree related methods 
    //@{
    /// # domains
    domain_id n_dom() const noexcept
    { return TREE->n_dom(); }
    /// domain pter given domain id
    cp_domain domain(domain_id d) const
    { return TREE->domain(d); }
    /// # top-tree cells
    count_type n_topcell() const noexcept
    { return TREE->n_topcell(); }
    /// begin of top-tree cells
    static Cell begin_topcell() noexcept
    { return Tree::begin_topcell(); }
    /// end of top-tree cell
    Cell end_topcell() const noexcept
    { return TREE->end_topcell(); }
    /// first top-tree cell in reversed order: last top-tree cell
    Cell rbegin_topcell() const noexcept
    { return TREE->rbegin_topcell(); }
    /// end top-tree cell in reversed order
    static Cell rend_topcell() noexcept
    { return Tree::rend_topcell(); }
    /// first domain contributing to cell
    auto first_domain(const Cell c) const noexcept
      -> decltype(this->tree()->first_domain(c))
    { return TREE->first_domain(c); }
    /// domain of a given leaf
    domain_id domain_of(const Leaf l) const noexcept
    { return TREE->domain_of(l); }
    /// # domains contributing to cell (>1 only for non-branch top-tree cells)
    domain_id n_domain(const Cell c) const noexcept
    { return TREE->n_domain(c); }
    /// is a cell a top-tree cell?
    bool is_top(const Cell c) const noexcept
    { return TREE->is_top(c); }
    /// is a cell a non-branch top-tree cell?
    bool is_nonbranch_top(const Cell c) const noexcept
    { return TREE->is_nonbranch_top(c); }
    /// is a cell a non-branch sub-domain cell?
    bool is_sub(const Cell c) const noexcept
    { return TREE->is_sub(c); }
    /// is a cell a branch cell?
    bool is_branch(const Cell c) const noexcept
    { return TREE->is_branch(c); }
    //@}
#endif // OCTALTREE_USE_OPENMP

    /// \name loops over cells, leaves, and leaf ranges
    //@{
    /// loop all internal leaves
    /// \param[in] f  function called for each leaf
    template<typename FuncOfLeaf>
    void loop_leaves(FuncOfLeaf f) const noexcept(noexcept(f))
    { for(auto l=begin_leaf(); l!=end_leaf(); ++l) f(l); }
#ifdef OCTALTREE_USE_OPENMP
    /// useful chunk size for parallelism on leaves
    template<OMP::Schedule schedule>
    unsigned leaf_chunk() const noexcept
    { return octtree::default_chunk<schedule>(n_leaf()); }
    /// loop all internal leaves in parallel
    /// \param[in] f  function called in parallel for each leaf
    /// \param[in] c  size of chunk for parallel loop
    /// \note to be called from within a OMP parallel region
    /// \note implies a barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_leaves_omp(FuncOfLeaf f, unsigned c)
      const noexcept(noexcept(f))
    {
      OMP::for_each<schedule>([f](unsigned l)
			      { f(reinterpret_cast<Leaf const&>(l)); },
			      0,n_leaf(),c);
    }
    /// loop all internal leaves in parallel without implicit barrier
    /// \param[in] f  function called in parallel for each leaf
    /// \param[in] c  size of chunk for parallel loop
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_leaves_omp_nowait(FuncOfLeaf f, unsigned c)
      const noexcept(noexcept(f))
    {
      OMP::for_each_nowait<schedule>([f](unsigned l)
				     { f(reinterpret_cast<Leaf const&>(l)); },
				     0,n_leaf(),c);
    }
    /// loop all internal leaves in parallel
    /// \param[in] f  function called in parallel for each leaf
    /// \note For static schedule we use non-overlapping chunks, while for
    ///       dynamic and guided schedule we use chunks of size @a
    ///       n_leaf()/(32*OMP::TeamSize().
    /// \note to be called from within a OMP parallel region
    /// \note implies a synchronisation barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_leaves_omp(FuncOfLeaf f)
      const noexcept(noexcept(f))
    {
      OMP::for_each<schedule>([f](unsigned l)
			      { f(reinterpret_cast<Leaf const&>(l)); },
			      0,n_leaf(),leaf_chunk<schedule>());
    }
    /// loop all internal leaves in parallel without implicit barrier
    /// \param[in] f  function called in parallel for each leaf
    /// \note For static schedule we use non-overlapping chunks, while for
    ///       dynamic and guided schedule we use chunks of size @a
    ///       n_leaf()/(32*OMP::TeamSize().
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_leaves_omp_nowait(FuncOfLeaf f)
      const noexcept(noexcept(f))
    {
      OMP::for_each_nowait<schedule>([f](unsigned l)
				     { f(reinterpret_cast<Leaf const&>(l)); },
				     0,n_leaf(),leaf_chunk<schedule>());
    }
#endif// OCTALTREE_USE_OPENMP
    /// loop all 4-blocks of internal leaves
    /// \param[in] f  function called for each 4-block of leaves
    /// \note the last 4-block may not be complete, though data are allocated
    template<typename FuncOfLeaf>
    void loop_blocks(FuncOfLeaf f) const noexcept(noexcept(f))
    { for(auto l=begin_leaf(); l<end_leaf(); l+=4) f(l); }
    /// loop all 4-blocks of internal leaves in range
    /// \param[in] f  function called for each 4-block of leaves
    /// \note the last 4-block may not be complete, though data are allocated
    template<typename FuncOfLeaf>
    void loop_blocks(const LeafRange r, FuncOfLeaf f)
      const noexcept(noexcept(f))
    { for(auto l=begin_leaf(r); l<end_leaf(r); l+=4) f(l); }
#ifdef OCTALTREE_USE_OPENMP
    /// useful chunk size for parallelism on 4-blocks of leaves
    template<OMP::Schedule schedule>
    unsigned block_chunk()  const noexcept
    { return octtree::default_chunk<schedule>(n_4blck()); }
    //    /// loop all 4-blocks of internal leaves in parallel
    /// \param[in] f  function called in parallel for each 4-block of leaves
    /// \param[in] c  size of chunk for parallel loop
    /// \note to be called from within a OMP parallel region
    /// \note implies a synchronisation barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_blocks_omp(FuncOfLeaf f, unsigned c)
      const noexcept(noexcept(f))
    {
      OMP::for_each<schedule>([f](unsigned l) { f(Leaf(l<<2)); },
			      0, n_4blck(), c);
    }
    /// loop all 4-blocks of internal leaves in parallel without implicit
    /// barrier
    /// \param[in] f  function called in parallel for each 4-block of leaves
    /// \param[in] c  size of chunk for parallel loop
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_blocks_omp_nowait(FuncOfLeaf f, unsigned c)
      const noexcept(noexcept(f))
    {
      OMP::for_each_nowait<schedule>([f](unsigned l) { f(Leaf(l<<2)); },
				     0, n_4blck(), c);
    }
    /// loop all 4-blocks of internal leaves in parallel
    /// \param[in] f  function called in parallel for each 4-block of leaves
    /// \note For static schedule we use non-overlapping chunks, while for
    ///       dynamic and guided schedule we use chunks of size @a
    ///       n_leaf()/(32*OMP::TeamSize().
    /// \note to be called from within a OMP parallel region
    /// \note implies a synchronisation barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_blocks_omp(FuncOfLeaf f)
      const noexcept(noexcept(f))
    {
      OMP::for_each<schedule>([f](unsigned l) { f(Leaf(l<<2)); },
			      0, n_4blck(), block_chunk<schedule>());
    }
    /// loop all 4-blocks of internal leaves in parallel without implicit
    /// barrier
    /// \param[in] f  function called in parallel for each 4-block of leaves
    /// \note For static schedule we use non-overlapping chunks, while for
    ///       dynamic and guided schedule we use chunks of size @a
    ///       n_leaf()/(32*OMP::TeamSize().
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_blocks_omp_nowait(FuncOfLeaf f)
      const noexcept(noexcept(f))
    {
      OMP::for_each_nowait<schedule>([f](unsigned l) { f(Leaf(l<<2)); },
				     0, n_4blck(), block_chunk<schedule>());
    }
#endif// OCTALTREE_USE_OPENMP
    /// loop all leaf ranges
    /// \param[in] f  function called for each range
    template<typename FuncOfLeafRange>
    void loop_ranges(FuncOfLeafRange f) const noexcept(noexcept(f))
    { for(auto r=begin_range(); r!=end_range(); ++r) f(r); }
    /// loop all internal leaves in range
    /// \param[in] f  function called for each cell
    template<typename FuncOfLeaf>
    void loop_leaves(const LeafRange r, FuncOfLeaf f)
      const noexcept(noexcept(f))
    { for(auto l=begin_leaf(r); l!=end_leaf(r); ++l) f(l); }
#ifdef OCTALTREE_USE_OPENMP
    /// loop all leaf ranges in parallel
    /// \param[in] f  function called in parallel for each leaf
    /// \param[in] c  size of chunk for parallel loop
    /// \note @a c=0 implies chunk=1 for dynamic and guided schedule and no
    ///       specified chunk size for static schedule.
    /// \note to be called from within a OMP parallel region
    /// \note implies a synchronisation barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeafRange>
    void loop_ranges_omp(FuncOfLeafRange f, unsigned c=0)
      const noexcept(noexcept(f))
    {
      OMP::for_each<schedule>([f](unsigned r) {
	  f(reinterpret_cast<LeafRange const&>(r)); },
	0,n_range(),c);
    }
    /// loop all leaf ranges in parallel without implicit barrier
    /// \param[in] f  function called in parallel for each leaf
    /// \param[in] c  size of chunk for parallel loop
    /// \note @a c=0 implies chunk=1 for dynamic and guided schedule and no
    ///       specified chunk size for static schedule.
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeafRange>
    void loop_ranges_omp_nowait(FuncOfLeafRange f, unsigned c=0)
      const noexcept(noexcept(f))
    {
      OMP::for_each_nowait<schedule>([f](unsigned r) {
	  f(reinterpret_cast<LeafRange const&>(r)); },
	0,n_range(),c);
    }
    /// loop leaves using parallelism on ranges
    /// \param[in] f  function called in parallel for each leaf
    /// \note to be called from within a OMP parallel region
    /// \note implies a synchronisation barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_leaves_ranged(FuncOfLeaf f, unsigned c=0)
      const noexcept(noexcept(f))
    {
      loop_ranges_omp<schedule>([this,f](const LeafRange r)
				{ this->loop_leaves(r,f); },
				c);
    }
    /// loop leaves using parallelism on ranges without implicit barrier
    /// \param[in] f  function called in parallel for each leaf
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_leaves_ranged_nowait(FuncOfLeaf f, unsigned c=0)
      const noexcept(noexcept(f))
    {
      loop_ranges_omp_nowait<schedule>([this,f](const LeafRange r)
				       { this->loop_leaves(r,f); },
				       c);
    }
    /// loop all 4-blocks of internal leaves using parallelism on ranges
    /// \param[in] f  function called in parallel for each 4-block of leaves
    /// \param[in] c  size of chunk for parallel loop
    /// \note to be called from within a OMP parallel region
    /// \note implies a synchronisation barrier before returning
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_blocks_ranged(FuncOfLeaf f, unsigned c=0)
      const noexcept(noexcept(f))
    {
      loop_ranges_omp<schedule>([this,f](const LeafRange r)
				{ this->loop_blocks(r,f); },
				c);
    }
    /// loop all 4-blocks of internal leaves using parallelism on ranges
    /// without implicit barrier
    /// \param[in] f  function called in parallel for each 4-block of leaves
    /// \param[in] c  size of chunk for parallel loop
    /// \note to be called from within a OMP parallel region
    template<OMP::Schedule schedule, typename FuncOfLeaf>
    void loop_blocks_ranged_nowait(FuncOfLeaf f, unsigned c=0)
      const noexcept(noexcept(f))
    {
      loop_ranges_omp_nowait<schedule>([this,f](const LeafRange r)
				       { this->loop_blocks(r,f); },
				       c);
    }
#endif// OCTALTREE_USE_OPENMP
    /// loop all external leaves
    /// \param[in] f  function called for each external leaf
    template<typename FuncOfExtLeaf>
    void loop_extleaves(FuncOfExtLeaf f) const noexcept(noexcept(f))
    { for(auto l=begin_extleaf(); l!=end_extleaf(); ++l) f(l); }
    /// loop all cells down (root first)
    /// \param[in] f  function called for each cell
    /// \note useful for downward pass: f(c) is called after f(parent(c))
    template<typename FuncOfCell>
    void loop_cells_down(FuncOfCell f) const noexcept(noexcept(f))
    { for(auto c=begin_cell(); c!=end_cell(); ++c) f(c); }
    /// loop all cells up (root last)
    /// \param[in] f  function called for each cell
    /// \note useful for upward pass: f(c) is called before f(parent(c))
    template<typename FuncOfCell>
    void loop_cells_up(FuncOfCell f) const noexcept(noexcept(f))
    { for(auto c=rbegin_cell(); c!=rend_cell(); --c) f(c); }
#ifdef OCTALTREE_USE_OPENMP
    /// loop non-branch top cells down (root first)
    /// \param[in] f  function called in parallel for each non-branch top cell
    /// \note useful for domain-parallel downward pass
    /// \note to be called from a single thread.
    template<typename FuncOfCell>
    void loop_topcells_down(FuncOfCell f) const noexcept(noexcept(f))
    {
      for(Cell c=begin_topcell(); c!=end_topcell(); ++c)
	if(is_nonbranch_top(c)) f(c);
    }
    /// loop non-branch top cells up (root last)
    /// \param[in] f  function called in parallel for each non-branch top cell
    /// \note useful for domain-parallel upward pass
    /// \note to be called from a single thread.
    template<typename FuncOfCell>
    void loop_topcells_up(FuncOfCell f) const noexcept(noexcept(f))
    {
      for(Cell c=rbegin_topcell(); c!=rend_topcell(); --c)
	if(is_nonbranch_top(c)) f(c);
    }
    /// loop all cells down (root first)
    /// \param[in] f  function called for each cell
    /// \note @a f(c) is called after @a f(parent(c)) (useful for downward pass)
    /// \note to be called from within OMP parallel region with as many threads
    ///       as there are domains
    template<typename FuncOfCell>
    void loop_cells_down_omp(FuncOfCell f) const noexcept(noexcept(f))
    {
      WDutilsAssertE(n_dom() == OMP::TeamSize());
#pragma omp single
      loop_topcells_down(f);
      domain(OMP::Rank())->loop_cells_down(f);
    }
    /// loop all cells up (root last)
    /// \param[in] f  function called for each cell
    /// \note f(c) is called before f(parent(c)) (useful for upward pass)
    /// \note to be called from within openMP parallel region with as many
    ///       threads as there are domains
    template<typename FuncOfCell>
    void loop_cells_up_omp(FuncOfCell f) const noexcept(noexcept(f))
    {
      WDutilsAssertE(n_dom() == OMP::TeamSize());
      domain(OMP::Rank())->loop_cells_up(f);
#pragma omp barrier
#pragma omp single
      loop_topcells_up(f);
    }
#endif
    /// loop daughter cells of a given cell
    /// \param[in] c  cell the daughters of which are looped over
    /// \param[in] f  function called for each cell
    template<typename FuncOfCell>
    void loop_daughters(const Cell c, FuncOfCell f) const noexcept(noexcept(f))
    { for(auto cc=begin_cell(c); cc!=end_cell(c); ++cc) f(cc); }
    /// find first daughter of a given cell matching a predicate
    /// \param[in] c  cell to find daughter for
    /// \param[in] p  predicate, must return bool
    /// \return first cell for which p(cell) returns true, of invalid cell
    template<typename Predicate>
    Cell find_first_daughter(const Cell c, Predicate p)
      const noexcept(noexcept(p))
    {
      static_assert(std::is_same<bool,
		    typename std::result_of<Predicate(Cell)>::type>::value,
		    "find_first_daughter() requires predicate returning bool");
      for(auto cc=begin_cell(c); cc!=end_cell(c); ++cc)
	if(p(cc)) return cc;
      return Cell::invalid();
    }
    /// loop daughter cells of a given cell in reverse order
    /// \param[in] c  cell the daughters of which are looped over
    /// \param[in] f  function called for each cell
    template<typename FuncOfCell>
    void loop_daughters_reverse(const Cell c, FuncOfCell f) const
      noexcept(noexcept(f))
    { for(auto cc=last_cell(c); cc!=rend_cell(c); ++cc) f(cc); }
    /// loop all pairs of daughter cells
    /// \param[in] c   cell the daughters of which are looped over
    /// \param[in] fc  function called for each daughter
    /// \param[in] fp  function called for each distinct pair of daughters
    template<typename FuncOfCell, typename FuncOfCellPair>
    void loop_daughter_pairs(const Cell c, FuncOfCell fc, FuncOfCellPair fp)
      const noexcept(noexcept(fc) && noexcept(fp))
    {
      for(auto ci=begin_cell(c); ci<end_cell(c); ++ci) {
	fc(ci);
	for(auto cj=ci.next(); cj<end_cell(c); ++cj) fp(ci,cj);
      }
    }
    /// loop leaf kids of a given cell
    /// \param[in] c  cell the leaf kids of which are looped over
    /// \param[in] f  function called for each leaf kid
    template<typename FuncOfLeaf>
    void loop_leaf_kids(const Cell c, FuncOfLeaf f) const noexcept(noexcept(f))
    { for(auto l=begin_leaf(c); l!=end_leaf_kids(c); ++l) f(l); }
    /// find first leaf kid, if any, satisfying predicate
    /// \param[in] c  cell the leaf kids of which are looped over
    /// \param[in] p  predicate: return first leaf kid for which f(l)==true
    template<typename Predicate>
    Leaf find_first_leaf_kid(const Cell c, Predicate p)
      const noexcept(noexcept(p))
    {
      for(auto l=begin_leaf(c); l!=end_leaf_kids(c); ++l)
	if(p(l)) return l;
      return Leaf::invalid();
    }
    /// loop leaf kids of a given cell using different function for first
    /// \param[in] c   cell the leaf kids of which are looped over
    /// \param[in] f0  function called for first leaf
    /// \param[in] f   function called for any further leaf child
    /// \note @a f0(begin_leaf(c)) is called regardless whether @c Cell @c has
    ///       leaf children or not
    /// \note useful for reducing over a cell
    template<typename FuncOfFirstLeaf, typename FuncOfLeaf>
    void loop_leaf_kids(const Cell c, FuncOfFirstLeaf f0, FuncOfLeaf f)
      const noexcept(noexcept(f0) && noexcept(f))
    { 
      auto l = begin_leaf(c);
      f0(l);
      for(++l; l<end_leaf_kids(c); ++l) f(l);
    }
    /// loop distinct leaf kids pairs
    /// \param[in] c  cell the leaf kids of which are looped over
    /// \param[in] f  function called for each distinct leaf kid pair
    template<typename FuncOfLeafPair>
    void loop_leaf_kid_pairs(const Cell c, FuncOfLeafPair f)
      const noexcept(noexcept(f))
    {
      for(auto li=begin_leaf(c); li<end_leaf_kids(c); ++li)
	for(auto lj=li.next(); lj<end_leaf_kids(c); ++lj)
	  f(li,lj);
    }
    /// loop all leaf descendants of a given cell
    /// \param[in] c  cell the leaf descendants of which are looped over
    /// \param[in] f  function called for each leaf
    template<typename FuncOfLeaf>
    void loop_all_leaves(const Cell c, FuncOfLeaf f) const noexcept(noexcept(f))
    { for(auto l=begin_leaf(c); l!=end_leaf_desc(c); ++l) f(l); }
    /// find first leaf descendent, if any, satisfying predicate
    /// \param[in] c  cell the leaf kids of which are looped over
    /// \param[in] p  predicate: return first leaf for which f(l)==true
    template<typename Predicate>
    Leaf find_first_leaf_desc(const Cell c, Predicate p)
      const noexcept(noexcept(p))
    {
      for(auto l=begin_leaf(c); l!=end_leaf_desc(c); ++l)
	if(p(l)) return l;
      return Leaf::invalid();
    }
    /// loop all distinct pairs of leaf descendants of given cell
    /// \param[in] c  cell the leaf descendants of which are looped over
    /// \param[in] f  function called for each distinct leaf pair
    template<typename FuncOfLeafPair>
    void loop_all_leaf_pairs(const Cell c, FuncOfLeafPair f)
      const noexcept(noexcept(f))
    {
      for(Leaf li=begin_leaf(c); li!=last_leaf_desc(c); ++li)
	for(Leaf lj=li.next(); lj!=end_leaf_desc(c); ++lj)
	  f(li,lj);
    }
    //@}

    ///
    /// \name support for upward passing
    //@{
  protected:
    /// passing data up for one tree cell
    ///
    /// \param[in] cell  cell to update
    /// \param[in] data  any data additional required or passed up
    ///
    /// \note see @a pass_up() for the requirements on @c Passer.
    template<typename Passer>
    void pass_up_cell(const Cell cell, typename Passer::Data*data=0) const
    {
      if(n_leaf_kids(cell)) {
	auto l=begin_leaf(cell);
	Passer passer(tree(),l,data);
	for(++l; l!=end_leaf_kids(cell); ++l)
	  passer.update(tree(),l,data);
	for(auto c=begin_cell(cell); c!=end_cell(cell); ++c)
	  passer.update(tree(),c,data);
	passer.set(tree(),cell,data);
      } else {
	auto c=begin_cell(cell);
	Passer passer(tree(),c,data);
	for(++c; c!=end_cell(cell); ++c)
	  passer.update(tree(),c,data);
	passer.set(tree(),cell,data);
      }
    }
  public:
#ifdef OCTALTREE_USE_OPENMP
    /// pass data up a domain
    /// \param[in] dom    id of domain to pass data up for
    /// \param[in] data   any data required for upward pass
    template<typename Passer>
    void pass_up_sub(domain_id dom, typename Passer::Data*data=0) const
    {
      domain(dom)->loop_cells_up([&](const Cell c)
				 { this->pass_up_cell<Passer>(c,data); }
				 );
    }
    /// pass data up the top-tree
    /// \note to be called *after* all domains have been passed
    /// \param[in] tree   tree to use (this is a static method)
    /// \param[in] data   any data required for upward pass
    template<typename Passer>
    void pass_up_top(typename Passer::Data*data=0) const
    {
      loop_topcells_up([&](const Cell c)
		       { this->pass_up_cell<Passer>(c,data); }
		       );
      Passer::final(tree(),data);
    }
#endif
    ///
    /// passing data up the tree
    ///
    /// \param[in] data  any data required by up-ward passing
    ///
    /// \note The type @c Passer must implement the following concept
    /// \code
    /// struct Passer {
    ///   typename Data;
    ///   Passer     (const Tree*, Leaf, Data*); // ctor from leaf kid
    ///   Passer     (const Tree*, Cell, Data*); // ctor from daughter cell
    ///   void update(const Tree*, Leaf, Data*); // update w.r.t. leaf kid
    ///   void update(const Tree*, Cell, Data*); // update w.r.t. daughter cell
    ///   void set   (const Tree*, Cell, Data*); // set cell
    ///   static void final(const Tree*, Data*); // called at the end
    /// }; \endcode
    ///
#ifdef OCTALTREE_USE_OPENMP
    template<typename Passer>
    void pass_up_omp(typename Passer::Data*data=0) const
    {
      if(OMP::IsParallel()) {
	// 1.1 already within a parallel team
	AssertNoThreads("pass_up_omp()");
	pass_up_sub<Passer>(OMP::Rank(),data);
# pragma omp barrier
# pragma omp single
	pass_up_top<Passer>(data);
      } else {
	// 1.2 envoke parallel team here
# pragma omp parallel
	{
	  AssertNoThreads("pass_up_omp()");
	  pass_up_sub<Passer>(OMP::Rank(),data);
	}
	pass_up_top<Passer>(data);
      }
    }
    /// serial version of @a pass_up_omp()
#endif// OCTALTREE_USE_OPENMP
    template<typename Passer>
    void pass_up_serial(typename Passer::Data*data=0) const
    {
      loop_cells_up([&](const Cell c)
		       { this->pass_up_cell<Passer>(c,data); }
		       );
      Passer::final(tree(),data);
    }
    //@}

    ///
    /// \name dump tree data
    //@{
  protected:
    /// header for leaf dump
    /// \note may be extended in derived class
    virtual std::ostream&LeafHeader(std::ostream&s) const
    { return TREE->LeafHeader(s); }
    /// dump leaf data
    /// \note may be extended in derived class
    virtual std::ostream&Dump(const Leaf l, std::ostream&s) const
    { return TREE->Dump(l,s); }
    /// dump external leaf data
    /// \note may be extended in derived class
    virtual std::ostream&Dump(const ExtLeaf l, std::ostream&s) const
    { return TREE->Dump(l,s); }
    /// header for cell dump
    /// \note may be extended in derived class
    virtual std::ostream&CellHeader(std::ostream&s) const
    { return TREE->CellHeader(s); }
    /// dump cell data
    virtual std::ostream&Dump(const Cell c, std::ostream&s) const
    { return TREE->Dump(c,s); }
  public:
    /// dump all leaf data
    //  NOTE: we cannot call the tree's method (need to catch virtual Dump())
    void DumpLeaves(std::ostream&s) const
    {
      LeafHeader(s) << '\n';
      loop_extleaves([&](const ExtLeaf l){ this->Dump(l,s) << '\n'; } );
      loop_leaves   ([&](const Leaf    l){ this->Dump(l,s) << '\n'; } );
      s.flush();
    }
    /// dump all cell data
    //  NOTE: we cannot call the tree's method (need to catch virtual Dump())
    void DumpCells(std::ostream&s) const
    {
      CellHeader(s) << '\n';
      loop_cells_down([&](const Cell c){ this->Dump(c,s) << '\n'; } );
      s.flush();
    }
    //@}

#ifdef OCTALTREE_USE_OPENMP
    /// \name miscellaneous
    //@{
    /// try to adjust # threads to # domains
    /// \param[in] func  name of calling function
    /// \note must not be called from within parallel region
    void AdjustNoThreads(const char*func=0) const
    { TREE->AdjustNoThreads(func); }
    /// asserts that # threads == # domains
    /// \param[in] func  name of calling function
    /// \note must be called from within parallel region
    void AssertNoThreads(const char*func=0) const
    { TREE->AssertNoThreads(func); }
    //@}
#endif
  };// TreeWalker<Tree derived from OctalTree>




















  //
  template<typename _Tree, class Enable=void> class StandardTreeWalk;
  ///
  /// The standard tree-walking algorithm.
  ///
  template<typename _Tree>
  class StandardTreeWalk<_Tree,
			 typename enable_if<is_OctalTree<_Tree>::value>::type>
    : public TreeWalker<_Tree>
  {
    typedef TreeWalker<_Tree>   Base;
    typedef typename Base::Tree Tree;
    typedef typename Base::Cell Cell;
    typedef typename Base::Leaf Leaf;
    typedef typename Base::ExtLeaf ExtLeaf;
    /// datum: a stack of cells
    mutable Stack<Cell> cell_stack;
    /// process many leaves using Proc::process(Leaf,Leaf)
    template<bool HasMany, typename Proc>
    typename std::enable_if< HasMany>::type
    process_many(const Proc*proc, Leaf l, const Leaf lN)
    { proc->process(l,lN); }
    /// process many leaves using Proc::process(Leaf)
    template<bool HasMany, typename Proc>
    typename std::enable_if<!HasMany>::type
    process_many(const Proc*proc, Leaf l, const Leaf lN)
    { for(; l<lN; ++l) proc->process(l); }
    //  template magic for determining whether
    //  \code void Proc::process(Leaf,Leaf) const \endcode
    //  exists
    template<typename T, T>
    struct sig_check : std::true_type {};
    //
    template<typename T, typename = std::true_type>
    struct proc_has_many : std::false_type {};
    //
    template<typename T> struct proc_has_many
    < T, std::integral_constant<bool, sig_check<void(T::*)(Leaf, Leaf) const,
						&T::process>::value>
    > : std::true_type {};
  public:
    /// ctor: initialise cell stack
    /// \param[in] tree  tree to walk
    explicit StandardTreeWalk(const Tree*tree)
      : Base(tree)
      , cell_stack(Tree::Nsub*this->depth()) {}
    ///
    /// tree walk
    ///
    /// \param[in] proc  (pter to) processor of tree nodes
    /// \note @c Proc must implement the following concept
    /// \code
    ///   // try to process cell, return true in case of success
    ///   bool Proc::process(Tree::Cell c) const;
    ///   // processs a single internal leaf
    ///   void Proc::process(Tree::Leaf l) const;
    ///   // optional: processs many internal leaves. if this function is not
    ///   // implemented, we process each leaf individually.
    ///   void Proc::process(Tree::Leaf l0, Tree::Leaf lN) const;
    ///   // should we also process external leaves (if any)?
    ///   bool Proc::do_external() const;
    ///   // processs a single external leaf
    ///   void Proc::process(Tree::ExtLeaf l) const;
    /// \endcode
    template<typename Proc>
    void operator()(const Proc*proc) const
    {
      // does @c Proc have member function @a void process(Leaf,Leaf) const ?
      static const bool ProcHasMany = proc_has_many<Proc>::value;
      // reset the cell stack
      cell_stack.reset();
      // ensure that the cell stack is big enough (tree depth may have changed)
      if(cell_stack.capacity()  < Tree::Nsub*this->depth())
	cell_stack.reset_capacity(Tree::Nsub*this->depth());
      // try processing the root cell, put it on stack if it needs opening
      if(!proc->process(this->root()))
	cell_stack.push(this->root());
      // while the stack is not empty: process cells
      while(!cell_stack.is_empty()) {
	const Cell c = cell_stack.pop();
	// process leaf kids of cell c
	if(this->n_leaf_kids(c))
	  process_many<ProcHasMany>(proc, this->begin_leaf(c),
					  this->end_leaf_kids(c));
	// process (or stack) daughters of cell c
	if(this->n_cell(c))
	  this->loop_daughters_reverse([&](const Cell cc){
	      if(!proc->process(cc))
		cell_stack.push(cc);
	    });
      }
      // process external leaves, if any
      if(proc->do_external() && this->n_extleaf())
	this->loop_extleaves([&](const ExtLeaf el){ proc->process(el); });
    }
  };// class StandardWalk<Tree>
  //
  template<typename _Tree >
  struct traits< StandardTreeWalk<_Tree> >
  { static const char*name()
    { return message("StandardTreeWalk<%s>",nameof(_Tree)); }
  };
} //  namespace WDutils
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_octtree_h
