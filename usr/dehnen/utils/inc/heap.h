// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/heap.h
///
/// \brief  provides support for min and max heap structures
///                                                                             
/// \author Walter Dehnen
///                                                                             
/// \date   2007,2009
///                                                                             
/// \version 14-06-2007 WD created from scratch 
/// \version 31-08-2007 WD transferred to utils
/// \version 08-12-2009 WD removed argument @a n from Walk::up()
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007,2009 Walter Dehnen
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
#ifndef WDutils_included_heap_h
#define WDutils_included_heap_h

#ifndef WDutils_included_algorithm
#  include <algorithm>
#  define WDutils_included_algorithm
#endif
#ifndef WDutils_included_functional
#  include <functional>
#  define WDutils_included_functional
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif

namespace WDutils {
  struct WalkBase {
    static unsigned parent(unsigned i) { return (i-1)>>1; }
    static void to_parent(unsigned&i) { --i; i>>=1; }
    static unsigned child (unsigned i) { return 1+(i<<1); }
    static void to_child(unsigned&i) { i<<=1; ++i; }
  };
  //
  //  struct WalkSwap
  //
  /// provides the elementary heap operations @a down() and @a up(), using @c
  /// std::swap(), as static member methods, to be used as template parameter
  /// for class @c HeapAlgorithms.
  ///
  /// \note The array holding the heap is addressed from 0 to n-1. This makes
  ///       the relations between parent and child indices somewhat arkward, but
  ///       fits much better with the C/C++ convention of array addressing.
  ///
  /// \note The template parameter @a Comparator for down() and up() specifies
  ///       the heap order. Use, for instance, @c std::less for a max heap and
  ///       @c std::greater for a min heap.
  struct WalkSwap : public WalkBase {
    /// downwards pass (heapify algorithm) using @c std::swap()
    /// \param[in,out] a        array[0..@a n -1] to be heapified
    /// \param[in]     p        index of element to be walked down the tree
    /// \param[in]     n        size of array
    /// \param[in]     compare  Comparator (see \<functional\>), e.g. std::less
    template<typename T, class Comparator>
    static void down(T*a, unsigned p, unsigned n, Comparator compare)
    {
      for(unsigned c=child(p); c<n; p=c,to_child(c)) {
	if((c+1)<n && compare(a[c],a[c+1])) ++c;
	if(compare(a[p],a[c])) std::swap(a[p],a[c]);
	else break;
      }
    } 
    /// upwards pass using @c std::swap()
    /// \param[in,out] a        array[0..@a n -1] to be heapified
    /// \param[in]     c        index of element to be walked up the tree
    /// \param[in]     compare  Comparator (see \<functional\>), e.g. std::less
    template<typename T, class Comparator>
    static void up(T*a, unsigned c, Comparator compare)
    {
      for(unsigned p=parent(c); c; c=p,to_parent(p)) {
	if(compare(a[p],a[c])) std::swap(a[p],a[c]);
	else break;
      }
    } 
  };
  //
  //  struct WalkNoSwap
  //
  /// provides the elementary heap operations @a down() and @a up() without @c
  /// std::swap() as static member methods, to be used as template parameter
  /// for class @c HeapAlgorithms.
  ///
  /// \note Instead of @c std::swap(), we use the copy constructor and assign
  ///       (=) operator of type @a T. For most types, i.e. when std::swap()
  ///       is not specialised, this algorithm is faster than that implemented
  ///       in @c WalkSwap. Thus, @c WalkSwap should be used only if the type
  ///       heaped has no copy constructor and assign operator or if @c
  ///       std::swap() for this type is more efficient (as for some container
  ///       classes).
  ///
  /// \note The array holding the heap is addressed from 0 to @a n -1. This
  ///       makes the relations between parent and child indices somewhat
  ///       arkward, but fits much better with the C/C++ convention of array
  ///       addressing.
  ///
  /// \note The template parameter @a Comparator for @a down() and @a up()
  ///       specifies the the heap order. Use, for instance, std::less for a
  ///       max heap and std::greater for a min heap.
  struct WalkNoSwap : public WalkBase {
    /// downwards pass (heapify algorithm) not using std::swap(), but copy
    /// \param[in,out] a        array[0..@a n -1] to be heapified
    /// \param[in]     p        index of element to be walked down the tree
    /// \param[in]     n        size of array
    /// \param[in]     compare  Comparator (see \<functional\>), e.g. std::less
    template<typename T, class Comparator>
    static void down(T*a, unsigned p, unsigned n, Comparator compare) {
      if(p>=(n>>1)) return;
      T tmp=a[p];
      unsigned i=p;
      for(unsigned c=child(p); c<n; p=c,to_child(c)) {
	if((c+1)<n && compare(a[c],a[c+1])) ++c;
	if(compare(tmp,a[c])) a[p]=a[i=c];
	else break;
      }
      a[i]=tmp;
    }
    /// upwards pass not using std::swap()
    /// \param[in,out] a array[0..@a n -1] to be heapified
    /// \param[in]     c        index of element to be walked up the tree
    /// \param[in]     compare  Comparator (see \<functional\>), e.g. std::less
    template<typename T, class Comparator>
    static void up(T*a, unsigned c, Comparator compare) {
      if(c<1) return;
      T tmp=a[c];
      unsigned i=c;
      for(unsigned p=parent(c); c; c=p,to_parent(p)) {
	if(compare(a[p],tmp)) a[c]=a[i=p];
	else break;
      }
      a[i]=tmp;
    } 
  };
  //
  // class HeapAlgorithms
  //
  /// provides (static member) methods for supporting heaps
  /// \note @a Walk to be either WalkSwap or WalkNoSwap
  /// \note @a Comparator to be usually std::less (for a max-heap) or
  ///       std::greater (for a min-heap)
  template<class Walk, template<typename T> class Comparator>
  class HeapAlgorithms {
  public:
    /// downwards pass
    /// \param[in,out] a        array[0..@a n -1] to be heapified
    /// \param[in]     p        index of element to be walked down the tree
    /// \param[in]     n        size of array
    template<typename T>
    static void down(T*a, unsigned p, unsigned n)
    { Walk::down(a,p,n,Comparator<T>()); }
    /// upwards pass
    /// \param[in,out] a        array[0..@a n -1] to be heapified
    /// \param[in]     c        index of element to be walked up the tree
    template<typename T>
    static void up(T*a, unsigned c)
    { Walk::up(a,c,Comparator<T>()); }
    /// put randomly ordered array into heap order
    /// \param[in,out] a  array[0..@a n -1]; input: random; output: heap order
    /// \param[in]     n  size of array
    template<typename T> 
    static void build(T*a, unsigned n)
    { for(unsigned i=n>>1; i; --i) down(a,i-1,n); }
    /// put array back into heap order after the top element has been replaced.
    /// This is equivalent to, but faster than, pop() followed by push().
    /// \param[in,out] a  array[0..@a n -1] in heap order, except for element 0 
    /// \param[in]     n  size of array
    template<typename T>
    static void after_top_replace(T*a, unsigned n)
    { down(a,0,n); }
    /// put array back into heap order after the last element has been replaced
    /// \note used in push()
    /// \param[in,out] a  array[0..@a n -1] in heap order, except for element @a
    ///                   n -1
    /// \param[in]     n  size of array
    template<typename T>
    static void after_last_replace(T*a, unsigned n)
    { if(n) up(a,n-1); }
    /// put array back into heap order after the @a i th element (only) has
    /// been replaced
    /// \param[in,out] a  array[0..@a n -1] in heap order, except for element @a
    ///                   i
    /// \param[in]     i  element that has been replaced, must be in [0,@a n -1]
    /// \param[in]     n  size of array
    template<typename T>
    static void after_replace(T*a, unsigned i, unsigned n)
    {
      if(i==0)
	return after_top_replace(a,n);
      unsigned p=Walk::parent(i);
      Comparator<T> compare;
      if(compare(a[p],a[i]))
	Walk::up(a,i,compare);
      else {
	unsigned c=Walk::child(i);
	if((c  <n && compare(a[i],a[c  ])) || 
	   (c+1<n && compare(a[i],a[c+1]))    )
	  Walk::down(a,i,n,compare);
      }
    }
    /// insert new element into heap
    /// \note array must allow for one additional element
    /// \param[in]     x  new element to be added
    /// \param[in,out] a  array[0..@a n -1] in heap order
    /// \param[in,out] n  size of array, increased by one on output
    template<typename T>
    static void push(T const&x, T*a, unsigned&n)
    {
      a[n++] = x;
      after_last_replace(a,n);
    }
    /// remove top element from heap, but preserve heap order
    /// \note number of elements is reduced by one
    /// \param[in]     x  top element of original heap
    /// \param[in,out] a  array[0..@a n -1] in heap order
    /// \param[in,out] n  size of array, reduced by one on output
    template<typename T>
    static void pop(T&x, T*a, unsigned&n)
    {
      if(n) {
	x = a[0];
	a[0] = a[n-1];
	after_top_replace(a,--n);
      }
    }
    /// transform heap structured array into ordered array.
    /// For a max-heap, i.e. if @a Comparator is std::less, the array will be
    /// ascendingly ordered, while a min-heap will be descendingly ordered
    /// \param[in,out] a  array[0..@a n -1]; input heap; output ordered
    /// \param[in]     n  size of array, assuming @a n > 0
    template<typename T>
    static void sort(T*a, unsigned n)
    {
      if(n)
	for(--n; n; --n) {
	  std::swap(a[n],a[0]);
	  down(a,0,n);
	}
    }
    /// find the index of the smallest element in a heap which is larger than x
    /// (for @a Comparator std::less).
    /// \param[in] x  key < top of heap
    /// \param[in] a  array[0..@a n -1] in heap order
    /// \param[in] n  size of array
    /// \return index of element closest to x, but upwards in Comparator order
    template<typename T>
    static unsigned next_up(T const&x, const T*a, unsigned n)
    {
      WDutilsAssert(compare(x,a[0]));     // ensure x not larger than top
      for(unsigned i,j,p=0;;) {
	i=j=Walk::child(p);
	if(i>=n) return p;                // no children -> return p
	if((i+1)<n) {                     // 2 children? i=larger, j=smaller
	  if(compare(a[i],a[i+1])) ++i; else ++j;
	}
	if(!compare(x,a[i])) return p;    // x >= larger -> return p
	p = compare(x,a[j])? j:i;         // x < smaller -> take p=smaller
      }                                   //   otherwise -> take p=larger
    }
    //--------------------------------------------------------------------------
    template<typename T> 
    static void build(WDutils::Array<T>&a)
    { build(a.array(),a.size()); }
    template<typename T> 
    static void after_top_replace(Array<T>&a)
    { after_top_replace(a.array(),a.size()); }
    template<typename T> 
    static void after_last_replace(Array<T>&a)
    { after_last_replace(a.array(),a.size()); }
    template<typename T> 
    static void after_replace(Array<T>&a, unsigned i)
    { after_replace(a.array(),i,a.size()); }
    template<typename T> 
    static void sort(Array<T>&a)
    { sort(a.array(),a.size()); }
  };// class HeapAlgorithms<Walk,Comparator>
  //
  /// support for max heaps, using std::swap()
  typedef HeapAlgorithms<WalkSwap,std::less> MaxHeapSwap;
  /// support for max heaps, avoiding std::swap()
  typedef HeapAlgorithms<WalkNoSwap,std::less> MaxHeap;
  /// support for min heaps, using std::swap()
  typedef HeapAlgorithms<WalkSwap,std::greater> MinHeapSwap;
  /// support for min heaps, avoiding std::swap()
  typedef HeapAlgorithms<WalkNoSwap,std::greater> MinHeap;
  //
  /// \name some simple sort algorithms put together from the above.
  //@{
  /// ascending heapsort in place, avoiding std::swap()
  template<typename T>
  void HeapSortAsc(T*a, unsigned n)
  {
    if(n>0) {
      MaxHeap::build(a,n);
      MaxHeap::sort(a,n);
    }
  }
  /// descending heapsort in place, avoiding std::swap()
  template<typename T>
  void HeapSortDesc(T*a, unsigned n)
  {
    if(n>0) {
      MinHeap::build(a,n);
      MinHeap::sort(a,n);
    }
  }
  /// ascending heapsort in place, using std::swap()
  template<typename T>
  void HeapSortAscSwap(T*a, unsigned n)
  {
    if(n>0) {
      MaxHeapSwap::build(a,n);
      MaxHeapSwap::sort(a,n);
    }
  }
  /// descending heapsort in place, using std::swap()
  template<typename T>
  void HeapSortDescSwap(T*a, unsigned n)
  {
    if(n>0) {
      MinHeapSwap::build(a,n);
      MinHeapSwap::sort(a,n);
    }
  }
  //@}
  // ///////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_heap_h
