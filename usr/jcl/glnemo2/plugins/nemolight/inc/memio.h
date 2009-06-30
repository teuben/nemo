/*
 * MEMIO.H: safe parsing of misaligned data
 */

#ifndef _memio_h
#define _memio_h

typedef struct mem {
    char *buf;      /* pointer to memory buffer */
    int buflen;     /* length of memory buffer */
    char *bp;       /* current pointer inside buffer */
} mem;

extern mem *memopen(char *buf, int buflen);
extern void memclose(mem *mp);
extern float memfread(mem *mp);
extern double memdread(mem *mp);
extern int memiread(mem *mp);
extern void memseed(mem *mp, int *loc, int mode);

#endif
