
/* 

NAME
  insertion sort 

DESCRIPTION

  Sorts by inserting each successive record into its proper place in
  the preceeding, already sorted records.   Perform a binary search on
  the preceeding records to find where to insert the current record,
  then shift all the records over to make room in the right place.

AUTHORSHIP

  This is properly called the binary insertion sort and credit is
  given for first publication to John Mauchly, 1946.

REFERENCES

  Knuth, Art of Computer Programming Vol 3, page 83.

COMPLEXITY

  O(n log n) comparisons 
  O(0.25 * n^2) memory operations
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

int insertion_sort(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int width;
  int (*cmpr_func)();
{
  int i, top, middle, bottom;
  int found, test;
  char * temp;

  temp = malloc(width);
  assert(temp != NULL);

  for (i = 1; i < nelem; i++)
  {
    /* binary search looking for place between base[0] and base[i-1] where
       you can insert base[i] into in order. */

    bottom = 0;
    top = i-1;
    middle = 0;
    found = 0;

    while (top >= bottom && ! found)
    {
      middle = (top+bottom) / 2;
      test = (*cmpr_func)(base+middle*width, base+i*width);

      if (test > 0)
        top = middle - 1;
      else if (test < 0)
      {
        middle ++;
        bottom = middle;
      }
      else
      {
        middle ++;
        found = 1;
      }
    }

    /* make room at base[middle] for base[i] */

    if (i != middle)
    {
      memcpy(temp, base+i*width, width);
      memmove(base+middle*width+width, base+middle*width, (i-middle)*width);
      memcpy(base+middle*width, temp, width);
    }
  }

  if (temp != NULL) free(temp);

  return 0;

}
