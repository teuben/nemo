
/* use the next line if your C library does not have memmove */

#define memmove(a, b, n) bcopy(b, a, n)

/* use the next line if your C library does not have memcpy */

/* #define memcpy(to, from, n) bcopy(from, to, n)  */


#define MAYBE_THRESHOLD 8
extern int _maybe_sorted;

extern void exit();
extern int qsort();

extern int merge_sort();
extern int heap_sort();
extern int bubble_sort();
extern int quick_sort();
extern int shell_sort();
extern int insertion_sort();
extern int bogo_sort();

/* It's more important that these represent common sizes on your
   machine than that they represent a certain number of bytes. */

typedef short               chunk2;
typedef long                chunk4;
typedef double              chunk8;

#define FLOGGER_VERSION 0
#define FLOGGER_PATCHLEVEL 0

