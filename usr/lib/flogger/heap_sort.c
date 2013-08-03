
/* 

NAME
  heap sort 

DESCRIPTION

  See some of the comments below for a partial description.

AUTHORSHIP

  This version was written by mouse@larry.mcrcim.mcgill.edu
  who has gracefully given me permission to include it in
  this package of sort algorithms. 

  The original heapsort is credited to Williams by Knuth and
  some improvements are credited to Floyd by Standish.

REFERENCES

  J. W. J. Williams, CACM 7 1964
  R. W. Floyd, Algorithm 245: Treesort 3, CACM 7:12 December, 1964
  Knuth, Art of Computer Programming Vol 3, page 145ff
  Thomas A. Standish, Data Structure Techniques, page 91-92

COMPLEXITY

  O(n log n)
  Iterative.

COPYRIGHT

  Questions about the copyright status of should be addressed
  to the author.

*/


#include <stdio.h>
#include <malloc.h>
#include <string.h>

/* this is an interface routine that takes a parameter
   list like qsort's and then calls the real heapsort. */
   
int heap_sort(base, nelem, width, cmpr_func)
  char * base;
  int nelem;
  int width;
  int (*cmpr_func)();
{
  char ** temp_in;
  char * temp_out;
  int i;

  temp_in = (char **) malloc(nelem * sizeof(char *));
  for (i = 0; i < nelem; i++)
    temp_in[i] = base+i*width;

  heapsort(temp_in, nelem, cmpr_func);

  temp_out = malloc(nelem * width);
  for (i = 0; i < nelem; i++)
    memcpy(temp_out+i*width, temp_in[i], width);

  memcpy(base, temp_out, nelem*width);

  free((char *) temp_in);
  free(temp_out);

  return 0;
}

/*
From uunet!Thunder.McRCIM.McGill.EDU!mouse Fri Nov 20 04:41:30 1992
Date: Fri, 20 Nov 92 05:13:19 -0500
From: der Mouse  <uunet!Thunder.McRCIM.McGill.EDU!mouse>
To: mikey@ontek.com
Subject: Re: qucik, merge, shell and bubble sort
Status: R

> I'll also be adding heapsort and shellsort algorithms; if anyone has
> a qsort-like drop-in replacement for either or both of those I'd
> appreciate a copy.

Note the interface difference (an extra level of pointers, which is a
bit of a pain because it often means allocating an array of char *),
but here's a heapsort I've been using.  (I wrote it.)
*/

#if 0
heapsort(ptrs,nels,cmp)
char **ptrs; /* should be `<unknown>**ptrs', but no such type exists */
int nels;
int (*cmp)();
#endif
/*
        Sorts the ptrs vector.  (*cmp)(ptrs[i],ptrs[j]) should return:

                < 0        if *ptrs[i] < *ptrs[j]
                = 0        if *ptrs[i] = *ptrs[j]
                > 0        if *ptrs[i] > *ptrs[j]

        For example, if the ptrs are actually pointers to int, it would
        be perfectly good to write a cmp function as follows (unless the
        integers are so large that overflow can occur in the subtraction):

                int cmp(p1,p2)
                char *p1;
                char *p2;
                {
                 return(*(int *)p1 - *(int *)p2);
                }

        Tip:  If the ptrs are character string pointers, the standard
        strcmp() function is a good cmp function.

        The vector will be in non-decreasing order by this criterion on
        return from heapsort.
*/

static _heapsort_bubble_up(size,ptrs,cmp)
int size;
char **ptrs;
int (*cmp)();
{
 int i;
 int p;
 char *temp;

 i = size;
 while (1)
  { if (i == 0)
     { return;
     }
    p = (i - 1) >> 1;
    if ((*cmp)(ptrs[i],ptrs[p]) > 0)
     { temp = ptrs[i];
       ptrs[i] = ptrs[p];
       ptrs[p] = temp;
       i = p;
     }
    else
     { return;
     }
  }
}

static _heapsort_bubble_down(size,ptrs,cmp)
int size;
char **ptrs;
int (*cmp)();
{
 int i;
 int j;
 int l;
 int r;
 int cl;
 int cr;
 char *temp;

 i = 0;
 while (1)
  { if (i >= size)
     { return;
     }
    l = i + i + 1;
    r = l + 1;
    cl = (l >= size) ? 1 : ((*cmp)(ptrs[i],ptrs[l]) >= 0);
    cr = (r >= size) ? 1 : ((*cmp)(ptrs[i],ptrs[r]) >= 0);
    switch ((cl<<1)|cr)
     { case 0:
          j = ((*cmp)(ptrs[l],ptrs[r]) > 0) ? l : r;
          break;
       case 1:
          j = l;
          break;
       case 2:
          j = r;
          break;
       case 3:
          return;
     }
    temp = ptrs[j];
    ptrs[j] = ptrs[i];
    ptrs[i] = temp;
    i = j;
  }
}

heapsort(ptrs,nels,cmp)
char **ptrs;
int nels;
int (*cmp)();
{
 int size;
 char *temp;

 if (nels <= 1)
  { return;
  }
 size = 1;
 while (size < nels)
  { _heapsort_bubble_up(size,ptrs,cmp);
    size ++;
  }
 while (size > 1)
  { size --;
    temp = ptrs[size];
    ptrs[size] = ptrs[0];
    ptrs[0] = temp;
    _heapsort_bubble_down(size,ptrs,cmp);
  }
}
