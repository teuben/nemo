/*
 * ALLOCATE: memory (re)allocation with error checking.
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
 */

#include <stdinc.h>
#include <errno.h>

#include "exception.h"

// #define calloc my_calloc
// #define NEW_EH

void *allocate(int nb)
{
    void *mem;
#ifdef NEW_EH
    if (BeginBlock() == -1)  {
	RestoreUserContext();
	return NULL;
    }
#endif
    // how should this kind of error be processed ?
    if (nb < 0) error("allocate: cannot allocate %d bytes",nb);
    if (nb==0) nb++;
    mem = (void *) calloc((size_t)nb, 1);
    if (mem == NULL)  {
#ifdef NEW_EH
	errno = ENOMEM;
	RaiseException (errno);
#else
	error("allocate: not enough memory for %d bytes", nb);
#endif
    }
    dprintf(8,"allocate: %d bytes @ %d \n",nb, mem);
#ifdef NEW_EH
    EndBlock();
#endif
    return mem;
}

void *reallocate(void *bp, int nb)
{
    void *mem;
#ifdef NEW_EH
    if (BeginBlock() == -1)  {
	RestoreUserContext();
	return NULL;
    }
#endif
    // how should this kind of error be processed ?
    if (nb < 0) error("reallocate: cannot allocate %d bytes",nb);
    if (nb == 0) nb++;
    if(bp==NULL)
        mem = (void *) calloc((size_t)nb, 1);
    else
        mem = (void *) realloc((void *)bp,(size_t)nb);
    if (mem == NULL)  {
#ifdef NEW_EH
	RaiseException (errno);
#else
	error("allocate: not enough memory for %d bytes", nb);
#endif
    }
    dprintf(8,"reallocate: %d bytes @ %d \n",nb, mem);
#ifdef NEW_EH
    EndBlock();
#endif
    return mem;
}

#undef calloc


void *
my_calloc(size_t nmemb, size_t size)
{
  static int my_counter = 0;

  my_counter++;
  if (my_counter == 5) return NULL;
  return calloc( nmemb, size);
}




