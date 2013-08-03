
/* 

NAME
  flogger

DESCRIPTION

  Put the sort routines through the paces.  Verify that they actually
  order data properly and that they don't molest the memory locations
  immediately above or below the array.  Tries some common worst-case
  situations like already-ordered data and comparison functions which
  are defective.

  Prints out estimates for each sort algorithm for usage of heap space,
  stack space, and run time and also for the number of calls to the
  comparison function.

  Takes one command line argument which specifies the number of items
  to put in the array.  Larger values take longer to run.

AUTHORSHIP

  Mike Lee, currently mikey@ontek.com

REFERENCES

  Knuth, Art of Computer Programming Vol 3: Searching and Sorting
  Kernighan & Richtie, The C Programming Language, Second Edition

WORK REMAINING

  The main loop is hopeless.
  See the TODO document (which should accompany this program.)

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
#include <signal.h>
#include <setjmp.h>

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

#include "sorting.h"

#define TEST_TYPE double
#define TEST_FUNC double_compare
#define DEFAULT_COUNT 2220 /* product of 4 and several primes */

#define CANARY_LOW -77.0
#define CANARY_HIGH -88.0

#define TEST_RANDOM 1
#define TEST_ASCEND 2
#define TEST_DESCEND 3
#define TEST_FIB_ASC 4
#define TEST_FIB_DESC 5
#define TEST_SURPRISE 6
#define TEST_MOSTLY 7
#define TEST_EQUIV 8

static int compare_count;
static int s_heap;
static char * s_low, * s_high;

#define CHECK_CANARIES \
  { \
    if (*(TEST_TYPE *) a == CANARY_LOW || *(TEST_TYPE *) b == CANARY_LOW)\
    { \
      printf("ran off the bottom of the array!\n"); \
      fflush(stdout); \
      longjmp(next, 1); \
    } \
    if (*(TEST_TYPE *) a == CANARY_HIGH || *(TEST_TYPE *) b == CANARY_HIGH)\
    { \
      printf("ran off the top of the array!\n"); \
      fflush(stdout); \
      longjmp(next, 1); \
    } \
  }

#define UPDATE_S_HEAP \
  { struct mallinfo mi; \
    mi = mallinfo();  \
    if (s_heap < mi.uordblks) s_heap = mi.uordblks;  \
    if (s_low == NULL || where < s_low) s_low = where; \
    if (s_high == NULL || where > s_high) s_high = where; \
    compare_count ++; }

static jmp_buf next;

void catch_quit()
{
  printf("\n");
  fflush(stdout);
  longjmp(next, 1);
}

void test_sort_func();

int int_compare(a, b) int * a, *b; 
{ 
  char foo;
  char * where = &foo;

  CHECK_CANARIES;
  UPDATE_S_HEAP;

  if (*a > *b) return 1;
  if (*a < *b) return -1;
  return 0;
}

int double_compare(a, b) double * a, *b; 
{ 
  char foo;
  char * where = &foo;

  CHECK_CANARIES;
  UPDATE_S_HEAP;

  if (*a > *b) return 1;
  if (*a < *b) return -1;
  return 0;
}

/*ARGSUSED*/
int lie_ascending(a, b) char * a, *b; 
{ 
  char foo;
  char * where = &foo;

  CHECK_CANARIES;
  UPDATE_S_HEAP;

  return -1;
}

/*ARGSUSED*/
int lie_descending(a, b) char * a, *b; 
{ 
  char foo;
  char * where = &foo;

  CHECK_CANARIES;
  UPDATE_S_HEAP;

  return 1;
}

/*ARGSUSED*/
int lie_equal(a, b) char * a, *b; 
{ 
  char foo;
  char * where = &foo;

  CHECK_CANARIES;
  UPDATE_S_HEAP;

  return 0;
}

/*ARGSUSED*/
int surprise(a, b) char * a, *b; 
{ 
  char foo;
  char * where = &foo;

  CHECK_CANARIES;
  UPDATE_S_HEAP;

  foo = rand() >> 23;
  if ((unsigned char) foo < 85) return -1;
  if ((unsigned char) foo < 171) return 0;
  return 1;
}

int main(argc, argv) int argc; char * argv[];
{
  int n;
  int stage = 0;
  int done = 0;

  printf("sort flogger version %d patch level %d.\n",
    FLOGGER_VERSION, FLOGGER_PATCHLEVEL);
  if (argc > 1) 
    n = atoi(argv[1]);
  else
    n = DEFAULT_COUNT;
  if (n < 0) n = -n;

  printf("timer resolution = %1.6f\n", 1.0/(double)HZ);
  printf("element size = %d\n", sizeof(TEST_TYPE));
  printf("number of elements = %d", n);
  fflush(stdout);

  setjmp(next);
  signal(SIGINT, catch_quit);

  /* apologies for the bizarre code layout */

  while (! done) switch (stage)
  {
    case 0: printf("\n*** qsort ***\n");
            printf("data         compares      stack       ");
            printf("heap       user       system\n");
            fflush(stdout);
            stage ++; test_sort_func(qsort, n, 1); break;
    case 1: stage ++; test_sort_func(qsort, n, 2); break;
    case 2: stage ++; test_sort_func(qsort, n, 3); break;
    case 3: stage ++; test_sort_func(qsort, n, 4); break;
    case 4: stage ++; test_sort_func(qsort, n, 5); break;
    case 5: stage ++; test_sort_func(qsort, n, 6); break;
    case 6: stage ++; test_sort_func(qsort, n, 7); break;
    case 7: stage ++; test_sort_func(qsort, n, 8); break;

    case 10: _maybe_sorted = 0;
             printf("\n*** merge_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(merge_sort, n, 1); break;
    case 11: stage ++; test_sort_func(merge_sort, n, 2); break;
    case 12: stage ++; test_sort_func(merge_sort, n, 3); break;
    case 13: stage ++; test_sort_func(merge_sort, n, 4); break;
    case 14: stage ++; test_sort_func(merge_sort, n, 5); break;
    case 15: stage ++; test_sort_func(merge_sort, n, 6); break;
    case 16: stage ++; test_sort_func(merge_sort, n, 7); break;
    case 17: stage ++; test_sort_func(merge_sort, n, 8); break;

    case 20: _maybe_sorted = 1;
             printf("\n*** merge_sort(_maybe_sorted) ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(merge_sort, n, 1); break;
    case 21: stage ++; test_sort_func(merge_sort, n, 2); break;
    case 22: stage ++; test_sort_func(merge_sort, n, 3); break;
    case 23: stage ++; test_sort_func(merge_sort, n, 4); break;
    case 24: stage ++; test_sort_func(merge_sort, n, 5); break;
    case 25: stage ++; test_sort_func(merge_sort, n, 6); break;
    case 26: stage ++; test_sort_func(merge_sort, n, 7); break;
    case 27: stage ++; test_sort_func(merge_sort, n, 8); break;

    case 30: printf("\n*** heap_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(heap_sort, n, 1); break;
    case 31: stage ++; test_sort_func(heap_sort, n, 2); break;
    case 32: stage ++; test_sort_func(heap_sort, n, 3); break;
    case 33: stage ++; test_sort_func(heap_sort, n, 4); break;
    case 34: stage ++; test_sort_func(heap_sort, n, 5); break;
    case 35: stage ++; test_sort_func(heap_sort, n, 6); break;
    case 36: stage ++; test_sort_func(heap_sort, n, 7); break;
    case 37: stage ++; test_sort_func(heap_sort, n, 8); break;

    case 40: printf("\n*** shell_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(shell_sort, n, 1); break;
    case 41: stage ++; test_sort_func(shell_sort, n, 2); break;
    case 42: stage ++; test_sort_func(shell_sort, n, 3); break;
    case 43: stage ++; test_sort_func(shell_sort, n, 4); break;
    case 44: stage ++; test_sort_func(shell_sort, n, 5); break;
    case 45: stage ++; test_sort_func(shell_sort, n, 6); break;
    case 46: stage ++; test_sort_func(shell_sort, n, 7); break;
    case 47: stage ++; test_sort_func(shell_sort, n, 8); break;

    case 50: printf("\n*** quick_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(quick_sort, n, 1); break;
    case 51: stage ++; test_sort_func(quick_sort, n, 2); break;
    case 52: stage ++; test_sort_func(quick_sort, n, 3); break;
    case 53: stage ++; test_sort_func(quick_sort, n, 4); break;
    case 54: stage ++; test_sort_func(quick_sort, n, 5); break;
    case 55: stage ++; test_sort_func(quick_sort, n, 6); break;
    case 56: stage ++; test_sort_func(quick_sort, n, 7); break;
    case 57: stage ++; test_sort_func(quick_sort, n, 8); break;

    case 60: printf("\n*** insertion_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(insertion_sort, n, 1); break;
    case 61: stage ++; test_sort_func(insertion_sort, n, 2); break;
    case 62: stage ++; test_sort_func(insertion_sort, n, 3); break;
    case 63: stage ++; test_sort_func(insertion_sort, n, 4); break;
    case 64: stage ++; test_sort_func(insertion_sort, n, 5); break;
    case 65: stage ++; test_sort_func(insertion_sort, n, 6); break;
    case 66: stage ++; test_sort_func(insertion_sort, n, 7); break;
    case 67: stage ++; test_sort_func(insertion_sort, n, 8); break;

    case 70: printf("\n*** bubble_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(bubble_sort, n, 1); break;
    case 71: stage ++; test_sort_func(bubble_sort, n, 2); break;
    case 72: stage ++; test_sort_func(bubble_sort, n, 3); break;
    case 73: stage ++; test_sort_func(bubble_sort, n, 4); break;
    case 74: stage ++; test_sort_func(bubble_sort, n, 5); break;
    case 75: stage ++; test_sort_func(bubble_sort, n, 6); break;
    case 76: stage ++; test_sort_func(bubble_sort, n, 7); break;
    case 77: stage ++; test_sort_func(bubble_sort, n, 8); break;

    case 80: if (n >= 15)
             {
               stage += 10;   
               break;        /* algorithm is factorial!  give up! */
             }
             printf("\n*** bogo_sort ***\n");
             printf("data         compares      stack       ");
             printf("heap       user       system\n");
             fflush(stdout);
             stage ++; test_sort_func(bogo_sort, n, 1); break;
    case 81: stage ++; test_sort_func(bogo_sort, n, 2); break;
    case 82: stage ++; test_sort_func(bogo_sort, n, 3); break;
    case 83: stage ++; test_sort_func(bogo_sort, n, 4); break;
    case 84: stage ++; test_sort_func(bogo_sort, n, 5); break;
    case 85: stage ++; test_sort_func(bogo_sort, n, 6); break;
    case 86: stage ++; test_sort_func(bogo_sort, n, 7); break;
    case 87: stage ++; test_sort_func(bogo_sort, n, 8); break;

    case 90: done = 1;

    default: stage++;
  }

  return 0;
}

void test_sort_func(func, n, situation)
  int             (*func)();
  int             n;
  int             situation;
{
  unsigned int i;
  static TEST_TYPE * foo = NULL;
  TEST_TYPE temp;
  int (*fp)();
  int old_heap;
  struct tms start, finish;

  /* in case of ctrl-c */
  if (foo != NULL) free((char *) foo);

  foo = (TEST_TYPE *) malloc((n+2) * sizeof(TEST_TYPE));
  if (foo == NULL)
  {
    printf("insufficient memory to conduct test.\n");
    exit(1);
  }

  /* store magic values above and below the array so that we
     can verify the the sort didn't run off either end */

  foo[0] = (TEST_TYPE) CANARY_LOW;
  foo[n+1] = (TEST_TYPE) CANARY_HIGH;

  srand((int) n);

  /* initialize the contents of the array appropriately for
     the situation we are presenting to the sort function */

  for (i = 1; i <= n; i ++) 
  {
    switch(situation)
    {
      case TEST_RANDOM:
      case TEST_FIB_ASC:
      case TEST_FIB_DESC:
      case TEST_SURPRISE:
        do {
          foo[i] = rand();
        } while (foo[i] == CANARY_LOW || foo[i] == CANARY_HIGH);
      break;

      case TEST_ASCEND:
      case TEST_MOSTLY: 
      case TEST_EQUIV: 
        foo[i] = i;
      break;

      case TEST_DESCEND: 
        foo[i] = n - i;
      break;
    }
  }

  if (situation == TEST_MOSTLY && n > 0)
  {
    temp = foo[1];
    foo[1] = foo[n];
    foo[n] = temp;
  }

  compare_count = 0;
  s_low = NULL;
  s_high = NULL;
  old_heap = mallinfo().uordblks;
  s_heap = old_heap;

  switch(situation)
  {
    case TEST_RANDOM:   printf("random:   "); fp = TEST_FUNC; break;
    case TEST_ASCEND:   printf("ascend:   "); fp = TEST_FUNC; break;
    case TEST_DESCEND:  printf("descend:  "); fp = TEST_FUNC; break;
    case TEST_FIB_ASC:  printf("fib asc:  "); fp = lie_ascending; break;
    case TEST_FIB_DESC: printf("fib desc: "); fp = lie_descending; break;
    case TEST_SURPRISE: printf("surprise: "); fp = surprise;  break;
    case TEST_MOSTLY:   printf("mostly:   "); fp = TEST_FUNC;  break;
    case TEST_EQUIV:    printf("equiv:    "); fp = lie_equal;  break;
    default: fp = NULL; break;
  }

  fflush(stdout);

  times(&start);

  /* actually call the sort function */
  (*func)((char *) &foo[1], (int) n, sizeof(TEST_TYPE), fp);

  times(&finish);

/*
  printf("\n");
  for (i = 0; i < n; i ++) printf(i % 5 == 4 ? "%11d\n" : "%11d ", foo[i]);
  printf("\n");
*/

  /* print out one row of data */

  printf("%11d", compare_count);
  printf("%11ld", (long) (s_high - s_low));
  printf("%11d", s_heap - old_heap);
  printf("%11.2f", 
    (double) (finish.tms_utime - start.tms_utime) / (double) HZ);
  printf("%11.2f", 
    (double) (finish.tms_stime - start.tms_stime) / (double) HZ);
  fflush(stdout);

  if (mallinfo().uordblks - old_heap != 0)
    printf(" (%d leak!)", mallinfo().uordblks - old_heap);
  printf("\n");

  /* check the areas immediately above and immediately below
     the array for contamination */

  if (foo[0] != CANARY_LOW)
    printf("wrote before beginning of array!\n");

  if (foo[n+1] != CANARY_HIGH)
    printf("wrote past end of array!\n");

  /* if appropriate, make sure the data was actually sorted */

  if (situation == TEST_RANDOM || situation == TEST_ASCEND ||
      situation == TEST_DESCEND || situation == TEST_MOSTLY)
    for (i = 2; i <= n; i ++)
    {
      if ((*TEST_FUNC)(&foo[i], &foo[i-1]) < 0) 
      {
        printf("out of order at position %d\n", i);
        break;
      }
    }
  
  /* check for sortedness, but for this case it's not an error
     just the normal behavior of the algorithm */

  if (situation == TEST_EQUIV)
  {
    int stable = 1;

    for (i = 2; i <= n && stable; i ++)
    {
      if ((*TEST_FUNC)(&foo[i], &foo[i-1]) < 0) stable = 0;
    }

    if (stable)
      printf("algorithm appears to be stable.\n");
    else
      printf("algorithm is not stable.\n");
  }
  

  if (foo != NULL) 
  {
    free((char *) foo);
    foo = NULL;
  }

  return;
}

