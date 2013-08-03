
#undef TEST
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#define RAND_MAX 32767

/*

In alt.folklore.computers, there has been a discussion of the worst
possible sorting algorithm, which has been called "bogosort".

My interpretation of the algorithm is this: Exchange two randomly
chosen elements of the array to be sorted.  Check to see if the array
is in order.  Iterate until done.

From what I've read, theoretical analysis of this algorithm gives it a
performance of O(n!), which means that the time to sort is
proportional to the factorial of the number of elements.  And since
the algorithm is random in nature, it could range from instantaneous
(only two entries out of order, and it happens to exchange them first)
to to infinite (it might never succeed).

So, for kicks I coded up a bogosort routine.

In my testing, I discovered that the mean time to sort an array of ten
integers was 75 seconds (25MHz 486, Unix, gcc 2.1, "optimized").

Extrapolating from this, assuming O(n!), gives 312 days to sort 15
integers and 1,593,378 years to sort 20 integers.  Someone with a much
faster machine than mine will have to verify these figures.

-- 
Richard Krehbiel                                 richk@grebyn.com
OS/2 2.0 will do for me until AmigaDOS for the 386 comes along...


===== cut here ====== bogosort.c =====================================

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#define TRUE 1
#define FALSE 0

#if 0
#define range_rand(range) ((long)rand() * (long)range / (RAND_MAX + 1))

#if ! defined(range_rand)
int range_rand(range) int range;
{
	return (long)rand() * (long)range / (RAND_MAX + 1);
}
#endif

#endif

/* bogo sort is sensitive to rand()'s repeating-lower-digits feature 
   plus the above macro and function gave me different results with
   and without -O, so I replaced them with this slow but serviceable
   function. */

int range_rand(range) 
  int range;
{
  int result;
  extern long random();

  do { 
    result = (int) random(); 
  } 
  while (result < 0);
  
  result = result % range;

  return result;
}


void swap_elem(elem1, elem2, elem_size) 
			   register char *elem1;
			   register char *elem2;
			   register size_t elem_size;

{
	/***
	int temp;

	temp = *(int*) elem1;
	*(int*) elem1 = *(int*) elem2;
	*(int*) elem2 = temp;
	***/

	while(elem_size-- > 0)
	{
		char temp;
		temp = *elem1;
		*elem1++ = *elem2;
		*elem2++ = temp;
	}

}

int in_order(base, n_elem, elem_size, compare)
   register char *base;
   register int n_elem;
   register size_t elem_size;
   int (*compare)();
{
	while(--n_elem > 0)
	{
		if(compare(base, base+elem_size) > 0)
			return FALSE;
		base += elem_size;
	}
	return TRUE;
}

/*
  The bogo() function:

  This function is called with the same arguments as qsort.  When it returns,
  the elements of the given array are in order.

  You may wish to call srand() before using bogosort.
*/

int bogo_sort(base, n_elem, elem_size, compare)
		  char *base;
		  int n_elem;
		  size_t elem_size;
		  int (*compare)();
{
	assert(n_elem <= RAND_MAX);

	while(!in_order(base, n_elem, elem_size, compare))
	{
		register char *elem1, *elem2;

		elem1 = base + (range_rand(n_elem) * elem_size);
		elem2 = base + (range_rand(n_elem) * elem_size);

		swap_elem(elem1, elem2, elem_size);
	}
}

#ifdef TEST

int array[100];		/* Up to 100 elements - no further */

int int_compare(i1, i2)
  char *i1; char *i2;
{
	return *(int *)i1 - *(int *)i2;
}

main(argc, argv)
  int argc; char *argv[];
{
	time_t now;
	int n_elem;
	int i;

	if(argc < 2)
	{
		fprintf(stderr, "useage: %s <number of elements>\n",
				argv[0]);
		exit(EXIT_FAILURE);
	}

	n_elem = atoi(argv[1]);
	if(n_elem > 100)
	{
		fprintf(stderr, 
  "No more than 100 elements, please (as if your life is that long...)\n");
		exit(EXIT_FAILURE);
	}

	now = time((time_t *)NULL);
	srandom(now);

	fputs("Starting array:\n", stdout);

	for(i = 0; i < n_elem; i++)
	{
		array[i] = random();
		printf("%d ", array[i]);
	}
	fputs("\n", stdout);

	bogo_sort((char *)array, n_elem, sizeof(int), int_compare);

	fputs("Ending array:\n", stdout);
	for(i = 0; i < n_elem; i++)
	{
		printf("%d ", array[i]);
	}
	fputs("\n", stdout);
}
#endif

