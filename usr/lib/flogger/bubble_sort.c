
/* 

NAME
  bubble sort

DESCRIPTION

  Traipse up and down the records until the entire array is in order
  exchanging records with their neighbors if the adjacent records are
  out of order with respect to each other.  A flag keeps track of
  whether exchanges were done; when an entire pass is made and no
  exchanges were performed, the sort is complete.  Also, after each
  pass, the element at the end of the pass is guaranteed to be in it's
  final place so the next pass excludes it.

  This algorithm reverses direction on each pass, so that a single item
  out of order won't force worst-case behavior.  Knuth refers to this
  as the "cocktail shaker sort."

AUTHORSHIP

  Unknown to me.

REFERENCES

  Knuth Vol 3, page 106-111

COMPLEXITY

  Best case O(n)
  Worst case O(n^2)
  Iterative.

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

int bubble_sort(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int width;
  int (*cmpr_func)();
{
  int top = nelem - 1;
  int bottom = 0;
  int i, did_swap = 1;
  char * temp;

  temp = malloc(width);
  assert(temp != NULL);

  while (top > bottom)
  {
    did_swap = 0;
    for (i = bottom; i < top; i++)
      if ((*cmpr_func)(base+i*width, base+(i+1)*width) > 0)
      {
        memcpy(temp, base+i*width, width);
        memcpy(base+i*width, base+(i+1)*width, width);
        memcpy(base+(i+1)*width, temp, width);
        did_swap = 1;
      }

    if (!did_swap) break;
    top--;

    did_swap = 0;
    for (i = top - 1; i >= bottom; i--)
      if ((*cmpr_func)(base+i*width, base+(i+1)*width) > 0)
      {
        memcpy(temp, base+i*width, width);
        memcpy(base+i*width, base+(i+1)*width, width);
        memcpy(base+(i+1)*width, temp, width);
        did_swap = 1;
      }

    if (!did_swap) break;
    bottom ++;
  }

  free(temp);
  return 0;

}
