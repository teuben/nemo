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
 *      25-may-07       use size_t to better match malloc() 	pjt/Pierre Fortin <pierre.fortin@oamp.fr>
 */

#include <stdinc.h>
#include <errno.h>

void *allocate(size_t nb)
{
    void *mem;

    /* how should this kind of error be processed ? */
    if (nb < 0) error("allocate < 0: cannot allocate %d bytes",nb);
    if (nb==0) nb++;       /* never allocate 0 bytes */
    mem = (void *) calloc(nb, 1);
    if (mem == NULL)  {
	nemo_dprintf(0,"solaris csh: limit datasize unlimited\n");
        nemo_dprintf(0,"solaris ksh: ulimit -d unlimited\n");
	error("allocate: not enough memory for %d bytes", nb);
    }
    nemo_dprintf(8,"allocate: %d bytes @ %d (0x%x)\n",nb, mem, mem);
    return mem;
}

void *reallocate(void *bp, size_t nb)
{
    void *mem;

    /* how should this kind of error be processed ? */
    if (nb < 0) error("reallocate: cannot allocate %d bytes",nb);
    if (nb == 0) nb++;
    if(bp==NULL)
        mem = (void *) calloc(nb, 1);
    else
        mem = (void *) realloc((void *)bp,nb);
    if (mem == NULL)  {
	error("reallocate: not enough memory for %d bytes", nb);
    }
    nemo_dprintf(8,"reallocate: %d bytes @ %d \n",nb, mem);
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
  "VERSION=1.0\n   7-sep-2005 PJT",
  NULL,
};

string usage = "(re)allocate benchmark";

void nemo_main(void) {
  int size, size0 = getiparam("size")*1024;
  int nalloc = getiparam("nalloc");
  int incr = getiparam("incr")*1024;
  int nrealloc = getiparam("nrealloc");
  int repeat = getiparam("repeat");
  int i;
  char *data;

  nemo_dprintf(0,"  Alloc:  %d * %d bytes\n",nalloc,size0);
  nemo_dprintf(0,"ReAlloc:  %d * %d bytes\n",nrealloc,incr);

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




