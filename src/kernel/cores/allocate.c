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
 */

#include <stdinc.h>

void *allocate(int nb)
{
    char *mem;

    dprintf(8,"allocate: %d bytes\n",nb);
    if (nb < 0) error("allocate: cannot allocate %d bytes",nb);
    if (nb==0) nb++;
    mem = (char *) malloc((size_t)nb);
    if (mem == NULL)
	error("allocate: not enough memory for %d bytes", nb);
    return mem;
}

void *reallocate(void *bp, int nb)
{
    char *mem;

    dprintf(8,"reallocate: %d bytes\n",nb);
    if (nb < 0) error("reallocate: cannot allocate %d bytes",nb);
    if (nb == 0) nb++;
    if(bp==NULL)
        mem = (char *) malloc((size_t)nb);
    else
        mem = (char *) realloc((char *)bp,(size_t)nb);
    if (mem == NULL)
	error("reallocate: not enuf memory (%d bytes)", nb);
    return mem;
}

