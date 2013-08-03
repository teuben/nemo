
/*

NAME
  quicksort, or partition exchange sort

DESCRIPTION

  Sort by partitioning, then sort the partitions, and so on.

AUTHORSHIP

  C. A. R. Hoare.
  This particular implementation is taken from K&R2.

REFERENCES

  Hoare, Comp. J. 5, 1962
  Knuth, Vol 3, page 114ff
  Kernighan & Ritchie The C Programming Language, Revised Edition, pg 87
  Standish, Data Structure Techniques pg 23-27
  Sedgewick, Implementing Quicksort Programs CACM 21:10 October 1978

COMPLEXITY

  Best case O(n log n)
  Worse case O(n^2)
  Recursive.

PORTABILITY PROBLEMS

  Same problems as with merge_sort.

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
#include <assert.h>

#include "sorting.h"

static char * copy_buffer = NULL;

static int quick_sortn();
static int quick_sort4();
static int quick_sort8();

int quick_sort(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int width;
  int (*cmpr_func)();
{
  copy_buffer = malloc(width); 
  assert(copy_buffer != NULL);

  if (width == sizeof(chunk4) && (long) base % sizeof(chunk4) == 0) 
    quick_sort4((chunk4 *) base, cmpr_func, 0, nelem-1);
  else if (width == sizeof(chunk8) && (long) base % sizeof(chunk8) == 0) 
    quick_sort8((chunk8 *) base, cmpr_func, 0, nelem-1);
  else 
    quick_sortn(base, width, cmpr_func, 0, nelem-1);

  free(copy_buffer);
  copy_buffer = NULL;
}


static quick_sortn(base, width, cmpr_func, left, right)
  char * base;
  int width;
  int (*cmpr_func)();
  int left, right;
{
  int i, last;
  int half = (left+right)/2;

  if (left >= right) return;

  memcpy(copy_buffer, base+width*left, width);
  memcpy(base+width*left, base+width*half, width);
  memcpy(base+width*half, copy_buffer, width);

  last = left;
  for (i = left + 1; i <= right; i++)
    if ((*cmpr_func)(base+width*i, base+width*left) < 0)
    {
      last ++;
      if (last == i) continue;
      memcpy(copy_buffer, base+width*last, width);
      memcpy(base+width*last, base+width*i, width);
      memcpy(base+width*i, copy_buffer, width);
    }

  memcpy(copy_buffer, base+width*left, width);
  memcpy(base+width*left, base+width*last, width);
  memcpy(base+width*last, copy_buffer, width);

  quick_sortn(base, width, cmpr_func, left, last-1);
  quick_sortn(base, width, cmpr_func, last+1, right);
}


static quick_sort4(base, cmpr_func, left, right)
  chunk4 * base;
  int (*cmpr_func)();
  int left, right;
{
  int i, last; 
  chunk4 temp;

  if (left >= right) return;

  temp = base[left];
  base[left] = base[(left+right)/2];
  base[(left+right)/2] = temp;

  last = left;
  for (i = left + 1; i <= right; i++)
    if ((*cmpr_func)(&base[i], &base[left]) < 0)
    {
      last ++;
      if (last == i) continue;
      temp = base[last];
      base[last] = base[i];
      base[i] = temp;
    }

  temp = base[left];
  base[left] = base[last];
  base[last] = temp;

  quick_sort4(base, cmpr_func, left, last-1);
  quick_sort4(base, cmpr_func, last+1, right);
}


static quick_sort8(base, cmpr_func, left, right)
  chunk8 * base;
  int (*cmpr_func)();
  int left, right;
{
  int i, last; 
  chunk8 temp;

  if (left >= right) return;

  temp = base[left];
  base[left] = base[(left+right)/2];
  base[(left+right)/2] = temp;

  last = left;
  for (i = left + 1; i <= right; i++)
    if ((*cmpr_func)(&base[i], &base[left]) < 0)
    {
      last ++;
      if (last == i) continue;
      temp = base[last];
      base[last] = base[i];
      base[i] = temp;
    }

  temp = base[left];
  base[left] = base[last];
  base[last] = temp;

  quick_sort8(base, cmpr_func, left, last-1);
  quick_sort8(base, cmpr_func, last+1, right);
}

