
/* 

NAME
  merge sort

DESCRIPTION

  Divide up the array into halves, then quarters, then eighths and so
  on until it can be divided no more.  Perform a "tape merge" on
  adjacent sections, building up larger ordered subsections until
  you're done.

  This version has some nifty features:

  1. It avoids mempcy for arrays of 4 byte or 8 bytes, so sorting
     pointers, ints, longs, float, doubles, or appropriately size tiny
     structures works quickly.  Someday I will add merge2 and merge16
     and round out the collection, but beyond that the savings are
     questionable.

  2. After the merges of subsections become larger than a certain
     number of records, a quick check is made to see if the last item
     of the first subsection is less than the first item of the second
     subsection; in this case the array is already in order there so
     there's no point in doing a merge of the two subsections.  This
     feature can be turned off by setting _maybe_sorted to 0.

  3. Array indexes and pointer math have been simplified to a degree.
     A sophisticated compiler would have no trouble performing strength
     reduction, but this way some improvement is available even if you
     don't cc -O.

  4. If it can't malloc the space it needs, it gives up and calls qsort().

AUTHORSHIP

  Knuth mentions Von Neumann as the first person who implemented a
  decent merge sort.  Merge sort as an algorithm probably predates 
  electronic computers.

REFERENCES

  John Von Neumann, Collected Works 5, 1963
  Knuth, Art of Computer Programming Vol 3, pg 159-168

COMPLEXITY

  O(n log n)

  This code is recursive, although the algorithm can be (and often is)
  implemented without recursion.

  Unlike many sort algorithms, this takes 2n memory.

PORTABILITY PROBLEMS

  My choice of sizes for special sort functions and for the
  MAYBE_THRESHOLD aren't necessarily ideal for all computers
  for all time.  
  
  The alignment check is hopelessly RISC-centric, but easy to remove.

COPYRIGHT

  Copyright 1992 Michael Lee.

  (1) Permission to use, copy, distribute, and make changes is granted
  providing that (a) that you do so without charging anyone a fee and
  (b) this copyright notice must be included verbatim in all copies and 
  distributions.  

  (2) There is no warranty for this program, to the extent permitted by
  applicable law.  This program is provided "AS IS" without warranty of
  any kind, either expressed or implied, including, but not limited to,
  the implied warranties of merchantability and fitness for a particular 
  purpose.  The entire risk as to the quality and performance of this 
  program is with the user.  Should this program prove defective, the 
  user assumes the cost of all necessary servicing, repair or correction.

  (3) In no event unless required by applicable law will the copyright
  holder be liable to the user for damages, including any general,
  special, incidental or consequential damages arising out of the use
  or inability to use this program (including but not limited to loss of
  data or data being rendered inaccurate or losses sustained by the
  user or third parties or a failure of this program to operate with any
  other programs), even if such holder or third party has been advised
  of the possibility of such damages.

  (4) Object code produced by a compiler from this code may be 
  incorporated into commercial products without being subject to 
  restrictions (1)(a) or (1)(b).

*/

#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <string.h>

#include "sorting.h"

#ifdef PARANOID
#include <assert.h>
#endif

int _maybe_sorted = 1;

static int merge4();
static int merge8();
static int mergen();

int merge_sort(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int width;
  int (*cmpr_func)();
{
  if (width == sizeof(chunk4) && (int) base % sizeof(chunk4) == 0)
    return merge4((chunk4 *) base, nelem, cmpr_func);
  else if (width == sizeof(chunk8) && (int) base % sizeof(chunk8) == 0) 
    return merge8((chunk8 *) base, nelem, cmpr_func);
  else 
    return mergen(base, nelem, width, cmpr_func);
}


static int merge4(base, nelem, cmpr_func)
  chunk4 * base;
  int nelem;
  int (*cmpr_func)();
{
  int split;
  int a, b;
  chunk4 * c;
  int in_order;
  static chunk4 * out = NULL;
  int mine = 0;

  if (out == NULL)
  {
    out = (chunk4 *) malloc(nelem * sizeof(chunk4));
    if (out == NULL)
    {
      qsort((char *) base, nelem, sizeof(chunk4), cmpr_func);
      return 0;
    }
    mine = 1;
  }

  split = (nelem+1) / 2;

  if (split > 1) 
    (void) merge4(base, split, cmpr_func);
  if (nelem - split > 1) 
    (void) merge4(base+split, nelem-split, cmpr_func);

  if (_maybe_sorted && nelem > MAYBE_THRESHOLD &&
      (*cmpr_func)(base+split, base+split-1) >= 0) 
  {
    if (mine)
    {
      free((char *) out);
      out = NULL;
    }
    return 0;
  }

  a = 0; 
  b = split;

  c = out;

  while(a < split)
  {
    if (b >= nelem)
      in_order = 1;
    else 
      in_order = ((*cmpr_func)(base+a, base+b) <= 0);
    
    if (in_order)
    {
      *c = *(base + a);
      c ++;
      a ++;
    }
    else
    {
      *c = *(base + b);
      c ++;
      b ++;
    }

#ifdef PARANOID
    assert(c - out <= nelem);
#endif

  }

#ifdef PARANOID
  assert(c - out <= nelem);
#endif

  for (a = 0; a < c - out; a ++) base[a] = out[a];

  if (mine)
  {
    free((char *) out);
    out = NULL;
  }
    
  return 0;
}



static int merge8(base, nelem, cmpr_func)
  chunk8 * base;
  int nelem;
  int (*cmpr_func)();
{
  int split;
  int a, b;
  chunk8 * c;
  int in_order;
  static chunk8 * out = NULL;
  int mine = 0;

  if (out == NULL)
  {
    out = (chunk8 *) malloc(nelem * sizeof(chunk8));
    if (out == NULL)
    {
      qsort((char *) base, nelem, sizeof(chunk8), cmpr_func);
      return 0;
    }
    mine = 1;
  }

  split = (nelem+1) / 2;

  if (split > 1) 
    (void) merge8(base, split, cmpr_func);
  if (nelem - split > 1) 
    (void) merge8(base+split, nelem-split, cmpr_func);

  if (_maybe_sorted && nelem > MAYBE_THRESHOLD &&
      (*cmpr_func)(base+split, base+split-1) >= 0) 
  {
    if (mine) 
    {
      free((char *) out);
      out = NULL;
    }
    return 0;
  }

  a = 0; 
  b = split;

  c = out;

  while(a < split)
  {
    if (b >= nelem)
      in_order = 1;
    else 
      in_order = ((*cmpr_func)(base+a, base+b) <= 0);
    
    if (in_order)
    {
      *c = *(base + a);
      c ++;
      a ++;
    }
    else
    {
      *c = *(base + b);
      c ++;
      b ++;
    }

#ifdef PARANOID
    assert(c - out <= nelem);
#endif

  }

#ifdef PARANOID
  assert(c - out <= nelem);
#endif

  for (a = 0; a < c - out; a ++) base[a] = out[a];
    
  if (mine)
  {
    free((char *) out);
    out = NULL;
  }

  return 0;
}


static int mergen(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int (*cmpr_func)();
  int width;
{
  int split, nw, sw;
  int a, b;
  char * c;
  int in_order;
  static char * out = NULL;
  int mine = 0;

  nw = nelem*width;

  if (out == NULL)
  {
    out = (char *) malloc(nw);
    if (out == NULL)
    {
      qsort(base, nelem, width, cmpr_func);
      return 0;
    }
    mine = 1;
  }

  split = (nelem+1) / 2;
  sw = split*width;

  if (split > 1) 
    (void) mergen(base, split, width, cmpr_func);
  if (nelem - split > 1) 
    (void) mergen(base+sw, nelem-split, width, cmpr_func);

  if (_maybe_sorted && nelem > MAYBE_THRESHOLD &&
      (*cmpr_func)(base+sw, base+sw-width) >= 0) 
  {
    if (mine) 
    {
      free(out);
      out = NULL;
    }
    return 0;
  }

  a = 0; 
  b = sw;
  c = out;

  while(a < sw)
  {
    if (b >= nw)
      in_order = 1;
    else 
      in_order = ((*cmpr_func)(base+a, base+b) <= 0);
    
    /* todo: try coalescing adjacent iterations into the
       same call to memcpy when in_order is the same */

    if (in_order)
    {
      memcpy(c, base+a, width);
      c += width;
      a += width;
    }
    else
    {
      memcpy(c, base+b, width);
      c += width;
      b += width;
    }

#ifdef PARANOID
    assert(c - out <= nelem * width);
#endif

  }

#ifdef PARANOID
  assert(c - out <= nelem * width); 
#endif

  memcpy(base, out, c - out);

  if (mine)
  {
    free(out);
    out = NULL;
  }

  return 0;
}
