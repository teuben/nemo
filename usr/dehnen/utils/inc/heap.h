// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   utils/inc/heap.h                                                    
///                                                                             
/// \brief  provides support for min and max heap structures                    
///                                                                             
/// \author Walter Dehnen                                                       
///                                                                             
/// \date   2007                                                                
///                                                                             
/// \version 14-06-2007 WD created from scratch                                 
/// \version 31-08-2007 WD transferred to utils                                 
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2007 Walter Dehnen                                             
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
    static int parent(int i) { return (i-1)>>1; }
    static void to_parent(int&i) { --i; i>>=1; }
    static int child (int i) { return 1+(i<<1); }
    static void to_child(int&i) { i<<=1; ++i; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  struct WalkSwap                                                           
  //                                                                            
  /// provides the elementary heap operations down() and up() using std::swap() 
  ///                                                                           
  /// To be used as template parameter for class HeapAlgorithms.\n              
  /// \note The array holding the heap is addressed from 0 to n-1. This makes   
  ///       the relations between parent and child indices somewhat arkward, but
  ///       fits much better with the C/C++ convention of array addressing.     
  /// \note The Comparator template parameter for down() and up() specifies the 
  ///       the heap order. Use, for instance, std::less for a max heap and     
  ///       std::greater for a min heap.                                        
  // ///////////////////////////////////////////////////////////////////////////
  struct WalkSwap : public WalkBase {
    /// downwards pass (heapify algorithm) using swap()
    /// \param a array[0..n-1] to be heapified
    /// \param p index of element to be walked down the tree
    /// \param n size of array
    /// \param compare Comparator (see <functional>), e.g. std::less
    template<typename T, class Comparator>
    static void down(T*a, int p, int n, Comparator compare) {
      for(int c=child(p); c<n; p=c,to_child(c)) {
	if((c+1)<n && compare(a[c],a[c+1])) ++c;
	if(compare(a[p],a[c])) std::swap(a[p],a[c]);
	else break;
      }
    } 
    /// upwards pass using swap()
    /// \param a array[0..n-1] to be heapified
    /// \param c index of element to be walked up the tree
    /// \param n size of array
    /// \param compare Comparator (see <functional>), e.g. std::less
    template<typename T, class Comparator>
    static void up(T*a, int c, int n, Comparator compare) {
      for(int p=parent(c); c; c=p,to_parent(p)) {
	if(compare(a[p],a[c])) std::swap(a[p],a[c]);
	else break;
      }
    } 
  };// struct WalkSwap
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  struct WalkNoSwap                                                         
  //                                                                            
  /// provides the elementary heap operations down() and up() without swap()    
  ///                                                                           
  /// To be used as template parameter for class HeapAlgorithms<>.\n            
  /// Instead of swap(), we use the copy constructor and assign (=) operator of 
  /// the value_type. If these don't exist (are private) or if swap() is more   
  /// efficient, use class WalkSwap.\n                                          
  /// \note The array holding the heap is addressed from 0 to n-1. This makes   
  ///       the relations between parent and child indices somewhat arkward, but
  ///       fits much better with the C/C++ convention of array addressing.     
  /// \note The Comparator template parameter for down() and up() specifies the 
  ///       the heap order. Use, for instance, std::less for a max heap and     
  ///       std::greater for a min heap.                                        
  // ///////////////////////////////////////////////////////////////////////////
  struct WalkNoSwap : public WalkBase {
    /// downwards pass (heapify algorithm) not using swap(), but copy
    /// \param a array[0..n-1] to be heapified
    /// \param p index of element to be walked down the tree
    /// \param n size of array
    /// \param compare Comparator (see <functional>), e.g. std::less
    template<typename T, class Comparator>
    static void down(T*a, int p, int n, Comparator compare) {
      if(p>=(n>>1)) return;
      T tmp=a[p];
      int i=p;
      for(int c=child(p); c<n; p=c,to_child(c)) {
	if((c+1)<n && compare(a[c],a[c+1])) ++c;
	if(compare(tmp,a[c])) a[p]=a[i=c];
	else break;
      }
      a[i]=tmp;
    }
    /// upwards pass not using swap()
    /// \param a array[0..n-1] to be heapified
    /// \param c index of element to be walked up the tree
    /// \param n size of array
    /// \param compare Comparator (see <functional>), e.g. std::less
    template<typename T, class Comparator>
    static void up(T*a, int c, int n, Comparator compare) {
      if(c<1) return;
      T tmp=a[c];
      int i=c;
      for(int p=parent(c); c; c=p,to_parent(p)) {
	if(compare(a[p],tmp)) a[c]=a[i=p];
	else break;
      }
      a[i]=tmp;
    } 
  };// struct WalkNoSwap
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class HeapAlgorithms                                                       
  //                                                                            
  /// provides (static member) methods for supporting heaps                     
  /// \param Walk either WalkSwap or WalkNoSwap                                 
  /// \param Comparator usually either std::less or std::greater                
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  template<class Walk, template<typename T> class Comparator>
  class HeapAlgorithms {
  public:
    /// put randomly ordered array into heap order
    /// \param a array[0..n-1]; input: random; output: heap order
    /// \param n size of array
    template<typename T> 
    static void build(T*a, int n) {
      for(int i=n>>1; i; --i)
	Walk::down(a,i-1,n,Comparator<T>());
    }
    /// put array back into heap order after the top (first) element has been
    /// replaced.
    /// This is equivalent to, but much faster than, pop() followed by push().
    /// \param a array[0..n-1] in heap order, except for element 0
    /// \param n size of array
    template<typename T>
    static void after_top_replace(T*a, int n) {
      Walk::down(a,0,n,Comparator<T>());
    }
    /// put array back into heap order after the last element has been replaced
    /// \note used in push()
    /// \param a array[0..n-1] in heap order, except for element n-1
    /// \param n size of array
    template<typename T>
    static void after_last_replace(T*a, int n) {
      Walk::up(a,n-1,n,Comparator<T>());
    }
    /// put array back into heap order after the ith element (only) has
    /// been replaced
    /// \param a array[0..n-1] in heap order, except for ith element
    /// \param i element that has been replaced, must be in [0,n-1]
    /// \param n size of array
    template<typename T>
    static void after_replace(T*a, int i, int n) {
      int p=Walk::parent(i);
      Comparator<T> compare;
      if(compare(a[p],a[i]))
	Walk::up(a,i,n,compare);
      else {
	int c=Walk::child(i);
	if(c  <n && compare(a[i],a[c]) || 
	   c+1<n && compare(a[i],a[c+1]))
	  Walk::down(a,i,n,compare);
      }
    }
    /// insert new element into heap
    /// \note array must allow for one additional element
    /// \param x new element to be added
    /// \param a array[0..n-1] in heap order
    /// \param n size of array, increased by one on output
    template<typename T>
    static void push(T const&x, T*a, int&n) {
      a[n++] = x;
      after_last_replace(a,n);
    }
    /// remove top element from heap, but preserve heap order
    /// \note number of elements is reduced by one
    /// \param x top element of original heap
    /// \param a array[0..n-1] in heap order
    /// \param n size of array, reduced by one on output
    template<typename T>
    static void pop(T&x, T*a, int&n) {
      x = a[0];
      a[0] = a[n-1];
      after_top_replace(a,--n);
    }
    /// transform heap structured array into ordered array.
    /// For a max-heap, ie. if Comparator is std::less, the array will be
    /// ascendingly ordered, while a min-heap will be descendingly ordered
    /// \param a array[0..n-1]; input heap; output ordered
    /// \param n size of array
    template<typename T>
    static void sort(T*a, int n) {
      for(int i=n-1; i; --i) {
	swap(a[i],a[0]);
	Walk::down(a,0,--n,Comparator<T>());
      }
    }
    /// find next one up
    /// \param x key < top of heap
    /// \param a array[0..n-1] in heap order
    /// \param n size of array
    /// \return index of element closest to x, but upwards in Comparator order
    template<typename T>
    static int next_up(T const&x, const T*a, int n) {
      if(!compare(x,a[0])) return -1;     // x larger than top: illegal
      for(int i,j,p=0;;) {
	i=j=Walk::child(p);
	if(i>=n) return p;                // no children -> return p
	if((i+1)<n)                       // 2 children? i=larger, j=smaller
	  if(compare(a[i],a[i+1])) ++i; else ++j;
	if(!compare(x,a[i])) return p;    // x >= larger -> return p
	p = compare(x,a[j])? j:i;         // x < smaller -> take p=smaller
      }                                   //   otherwise -> take p=larger
    }
    //--------------------------------------------------------------------------
    template<typename T> 
    static void build(Array<T>&a) {
      build(a.array(),a.size()); }
    template<typename T> 
    static void after_top_replace(Array<T>&a) {
      after_top_replace(a.array(),a.size()); }
    template<typename T> 
    static void after_last_replace(Array<T>&a) {
      after_last_replace(a.array(),a.size()); }
    template<typename T> 
    static void after_replace(Array<T>&a, int i) {
      after_replace(a.array(),i,a.size()); }
    template<typename T> 
    static void sort(Array<T>&a) {
      sort(a.array(),a.size()); }
  };// class HeapAlgorithms<Walk,Comparator>
  // ///////////////////////////////////////////////////////////////////////////
  /// support for max heaps, using swap()
  typedef HeapAlgorithms<WalkSwap,std::less> MaxHeapSwap;
  /// support for max heaps, avoiding swap()
  typedef HeapAlgorithms<WalkNoSwap,std::less> MaxHeap;
  /// support for min heaps, using swap()
  typedef HeapAlgorithms<WalkSwap,std::greater> MinHeapSwap;
  /// support for min heaps, avoiding swap()
  typedef HeapAlgorithms<WalkNoSwap,std::greater> MinHeap;
  // ///////////////////////////////////////////////////////////////////////////
  /// \name some simple sort algorithms put together from the above.
  //@{
  /// ascending heapsort in place, avoiding swap()
  template<typename T>
  void HeapSortAsc(T*a, int n) {
    MaxHeap::build(a,n);
    MaxHeap::sort(a,n);
  }
  /// descending heapsort in place, avoiding swap()
  template<typename T>
  void HeapSortDesc(T*a, int n) {
    MinHeap::build(a,n);
    MinHeap::sort(a,n);
  }
  /// ascending heapsort in place, using swap()
  template<typename T>
  void HeapSortAscSwap(T*a, int n) {
    MaxHeapSwap::build(a,n);
    MaxHeapSwap::sort(a,n);
  }
  /// descending heapsort in place, using swap()
  template<typename T>
  void HeapSortDescSwap(T*a, int n) {
    MinHeapSwap::build(a,n);
    MinHeapSwap::sort(a,n);
  }
  //@}
  // ///////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_heap_h
