/*
 * ALLOCATE: memory (re)allocation with fatal error checking.
 *
 *	<dark ages>	created			Josh?
 *	 4-oct-90	added realloc		Peter
 *      16-nov-90       made realloc understand NULL pointers   PJT
 *	25-feb-92       somebody commented out malloc()?
 *      15-apr-92       AIX			
 *      26-jan-94       extra error() when nb < 1       pjt
 *	 3-jul-94       hmmm, allow 0 bytes, but alloc 1!   pjt
 *	22-jan-95	proto
 *	 3-may-95	added dprintf() 
 *	 6-apr-01	changed malloc -> calloc		pjt
 *         jan-02       experimenting with exception handling	pjt/nas
 *       7-sep-05       TOOLBOX benchmark pjt/chapman
 *      31-may-07       use size_t to better match malloc() 	pjt/Pierre Fortin <pierre.fortin@oamp.fr>
 *      12-jun-08       allocate_FL etc, see stdinc.h           WD
 *      12-jun-08       removed tests for size_t < 0            WD
 *      03-oct-08       debugged error in debug_info reporting  WD
 *       4-jan-11       add local/static arrays to show where they go in the TESTBED version
 */

#include <stdinc.h>
#include <errno.h>

void *allocate_FL(size_t nb, const_string file, int line)
{
    void *mem;

    if (sizeof(size_t) == 4 && nb > 2147483647)
      warning("allocate: 32bit machine allocate");

    if (nb==0) nb++;       /* never allocate 0 bytes */

    mem = (void *) calloc(nb, 1);

    if (mem == NULL)  {
	nemo_dprintf(0,"solaris csh: limit datasize unlimited\n");
        nemo_dprintf(0,"solaris ksh: ulimit -d unlimited\n");
	if(file) error("[%s:%d]: cannot allocate %lu bytes",file,line,nb);
	else     error("cannot allocate %lu bytes",nb);
    }

    if(file)
	nemo_dprintfN(8,"[%s:%d]: allocated %lu bytes @ %p\n",file,line,nb,mem);
    else
	nemo_dprintfN(8,"allocated %lu bytes @ %p\n",nb,mem);

    return mem;
}

void *reallocate_FL(void *bp, size_t nb, const_string file, int line)
{
    void *mem;

    if (nb == 0) nb++;

    if(bp==NULL)
        mem = (void *) calloc(nb, 1);
    else
        mem = (void *) realloc((void *)bp,nb);
    if (mem == NULL)  {
	if(file) error("[%s:%d]: cannot reallocate %lu bytes",file,line,nb);
	else     error("cannot reallocate %lu bytes",nb);
    }

    if(file)
	nemo_dprintfN(8,"[%s:%d]: reallocated %lu bytes @ %p\n",file,line,nb,mem);
    else
	nemo_dprintfN(8,"reallocated %lu bytes @ %p\n",nb,mem);

    return mem;
}

void *
my_calloc(size_t nmemb, size_t size)
{
  static int my_counter = 0;

  my_counter++;
  if (my_counter == 5) return NULL;
  return calloc( nmemb, size);
}

#if defined(TESTBED)

#include <nemo.h>

string defv[] = {
  "size=16\n       Size of a single bucket (in kB) to allocate",
  "nalloc=0\n      Number of extra times to allocate <size>",
  "incr=16\n       Increment size (in kB) to reallocate the size\n",
  "nrealloc=0\n    Number of times to increment and reallocate\n",
  "repeat=1\n      How often to repeat the whole test",
  "big1=100\n      Allocate the product of these two",
  "big2=100\n      Allocate the product of these two",
  "doubling=f\n    Doubling the big1*big2 allocation until failure",
  "VERSION=2.1\n   4-jan-2011 PJT",
  NULL,
};

string usage = "(re)allocate benchmark";

#define MAXTEST 128

static int test1[MAXTEST];

void nemo_main(void) {
  int size, size0 = getiparam("size")*1024;
  int nalloc = getiparam("nalloc");
  int incr = getiparam("incr")*1024;
  int nrealloc = getiparam("nrealloc");
  int repeat = getiparam("repeat");
  int big1 = getiparam("big1");
  int big2 = getiparam("big2");
  int big = big1*big2;
  size_t big64 = (size_t)big1*(size_t)big2;    /* this is crucial to cast */
  size_t i64;
  int i;
  int test2[MAXTEST];
  static int test3[MAXTEST];
  bool Qdouble = getbparam("doubling");
  char *data;

  nemo_dprintf(0,"static test1 @ %p\n",test1);
  nemo_dprintf(0,"       test2 @ %p\n",test2);
  nemo_dprintf(0,"static test3 @ %p\n",test3);


  nemo_dprintf(0,"  Alloc:  %d * %d bytes\n",nalloc,size0);
  nemo_dprintf(0,"ReAlloc:  %d * %d bytes\n",nrealloc,incr);

  if (big   < 0) warning("big   < 0: %d overflow? (%d x %d) %ld",big,big1,big2,big64);
  if (big64 < 0) warning("big64 < 0: %d overflow? (%d x %d)"    ,big,big1,big2);
  dprintf(0,"sizeof(size_t) = %d\n",sizeof(size_t));
  data = allocate(big64);
  free(data);
  dprintf(0,"Passed big64 allocating %ld\n",big64);
  while (Qdouble) {
    big64 *= 2;
    data = allocate(big64);
    for (i64=0; i64<big64; i64++)
      data[i64] = 0;
    free(data);
    dprintf(0,"Passed big64 allocating %ld\n",big64);
  }


  while (repeat-- > 0) {              /* repeat loop */

    size = size0;                     /* allocate loop */
    data = allocate(size);
    if (nalloc > 1) {
      for (i=0; i<nalloc; i++) {
	free(data);
	size += size0;
	data = allocate(size);
      }
      free(data);
    }
    
    size = size0;                      /* reallocate loop */
    data = reallocate(0,size);
    nemo_dprintf(1,"%d %d %d\n",repeat,size,incr);
    if (nrealloc > 1) {
      for (i=0; i<nrealloc; i++) {
	size += incr;
	data = reallocate(data,size);
      }
    }
    free(data);
  }
}

#endif




