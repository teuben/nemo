/*
 * ISMAX: find maximum in an array
 *	Can be optimized for Cray's
 *
 * Based on a Usenet discussion by:
 *   Tim Cullip
 *   cullip@vangogh.cs.unc.edu
 *   (Usenet - comp.unix.cray - October 1990)
 */

#include <stdinc.h>

real find_max1(a, limit) 
real *a; int limit;
{  
   real max;
   int i;

   max = a[0];
   for (i = 1; i < limit; i++) {
      if (a[i] > max) max = a[i];
   }
   return(max);
}

/* 
 *  The above code is straightforward, but doesn't vectorize because the
 *  test on max, and then possible change in the value of max creates a
 *  recurence. 
 *
 *  Using the #pragma ivdep directive doesn't solve the problem, since
 *  this recurence really does screw up the pipeline.
 */


real find_max2(a,limit) 
real *a; 
int limit;
{  
   real max;
   int i;

   max = a[0];
   i = 1;
   while (1) {
      while (i < limit) {
	 if (a[i] > max) break;
	    i++;
      }
      if (i == limit) break;
       max = a[i];
       i++;
   }
   return(max);
}

/*
 * This code seems to run about 10 times faster than the previous
 * method (on the type data that I need it to run on).  The inner
 * loop does vectorize, but everytime you find a new maximum you
 * fall out of the loop.
 */



#if defined(unicos)
real
find_max3(a,limit) real *a; int limit;
{  
    ISMAX(a,&limit);        /* Cray subroutine from LIBSCI - see sr-2081 */
}
#endif

#if 0
============== NOW SOME COMMENTS FROM THE NET ================================
   You have picked a very interesting example to start off with.  The
code you submit is actually nearly optimal (asymptotically, it is
optimal, if compiled correctly) in the case of a randomly-ordered
array.  Note that it can easily be modified to return the location of
the maximum as well as the value by saving it when it occurs.

   The alternative is to use the Cray library routine ISMAX(3SCI).
This has both advantages and disadvantages over your routine.  The big
disadvantage is that it will be a factor of two-plus slower in the
randomly-ordered case (as I discovered the first time I tried to use
it in an application where I cared about speed).

   It uses a search algorithm which is a common paradigm for Cray
programming.  Generate a vector with the first 64 elements of your
array.  Then read the next 64 elements; compare each to the corres-
ponding element of your vector, and choose the larger.  Repeat until
your array is exhausted.  One of the 64 elements of your vector will
be the maximum.  Of course, you have to deal with fragments which are
not multiples of 64, and further it is not trivial to write efficient
high-level code for this algorithm (the actual library routine is
written in assembler).


   It takes only a little knowledge of Cray hardware to realize that
the above algorithm, when perfectly optimized, is about two times
slower than your algorithm, as long as you stay in your inner loop.
Your code just does one-plus chimes of work on each vector, setting
the vector mask and then testing it for nonzero (the memory fetch and
floating subtract should be overlapped or chained); the Cray code does
two-plus chimes, setting the vector mask and then always a vector
merge.  The "plus" is then to go back and find the maximum, which
requires looking again at every 64th element.  If you want to be
clever you can avoid a stride of 64 here either by setting the vector
length to 63 or (more complicated) by using a nonunit stride near n/64
in the first pass.


In the random-long-vector case, therefore, your code will be much
better.  The break out of the vector loop only occurs, on average,
about lg n times out of n/64 vector iterations, which is
asymptotically negligible.  The Cray code has the distinct advantage,
however, that its worst case is the same as its average case.  Your
worst case is incredibly bad; if the array happens to be monotonically
increasing, I would expect your code to be hundreds of times slower
than in the random case.  You may also pay more overhead than you like
when the array is fairly short.  (Solving this problem efficiently on
short arrays is difficult.  But in general you should almost never
attack a problem this way in a time-critical application with short
arrays, but rather find 64 maxima in parallel.)


It might make sense to merge the two algorithms by implementing the
Cray algorithm but with a check for zero-vector-mask.  This would mean
that in the case that none of the 64 array elements is greater than
the corresponding vector element (still "almost always" in the
asymptotic random case), that a single chime would suffice, while the
worst case would not be much worse than the the existing Cray code.
Also, of course, if you want only the value of the maximum and not its
location, a version of the Cray code without that feature would be
slightly faster (this is more significant on a Cray-2 than on an
[XY]-MP, because the stride-64 problem is greater).

   -- David desJardins
=============================================================================
#endif
