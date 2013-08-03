
/*

NAME
  shell sort 

DESCRIPTION

  An exhange sort algorithm which exchanges over larger distances
  than bubble sort, thus reducing the number of exchanges needed.

AUTHORSHIP

  D. L. Shell

REFERENCES

  D. L. Shell, CACM 2, July 1959
  Kernighan & Ritchie, C Programming Language, Second Edition, pg 62
  Knuth, Art of Computer Programming Vol 3, pg 84-95

COMPLEXITY

  I'm not exactly sure, but I think it's O(N^1.5) 
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

int shell_sort(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int width;
  int (*cmpr_func)();
{
  int gap, i, j;
  char * temp;

  temp = malloc(width);
  assert(temp != NULL);

  for (gap = nelem / 2; gap > 0; gap /= 2)
  {
    for (i = gap; i < nelem; i++)
    {
      for (j = i-gap; 
           j >=0 && (*cmpr_func)(base+j*width, base+(j+gap)*width) > 0;
           j -= gap)
      {
        memcpy(temp, base+j*width, width);
        memcpy(base+j*width, base+(j+gap)*width, width);
        memcpy(base+(j+gap)*width, temp, width);
      }
    }
  }

  if (temp != NULL) free(temp);

  return 0;

}
