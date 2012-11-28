/*
 * EXTSTRING.C: an extension of the standard C concept of character
 * strings to strings of n-byte (nonzero) values, terminated by a
 * marker of n zero bytes.
 *	18-nov-91  malloc() into header declaration		PJT
 *	22-nov-91  allocate()
 *	25-feb-92  happy gcc2.0
 *       7-feb-95  prototyped
 */

#include <stdinc.h>
#include <extstring.h>

#define MAXLEN  1024

void *getxstr(                  /* returns ptr to ext. string */
    stream inpt,                /* stdio input stream */
    int nbyt)                   /* bytes per value */
{
    char buf[MAXLEN], *bp;
    bool lpflg;
    int i, ch;

    bp = &buf[0];                               /* point into buffer */
    do {                                        /* loop to read in */
        lpflg = FALSE;                          /*   init loop flag */
        for (i = 0; i < nbyt; i++) {            /*   loop over bytes */
            ch = getc(inpt);                    /*     input one byte */
            if (bp > &buf[MAXLEN-1])            /*     no space left? */
                error("getxstr: buffer overflow");
            *bp =(ch != EOF ? (char) ch : 0);   /*     map EOF to NULL */
            if (*bp++ != 0)                     /*     a byte of data? */
                lpflg = TRUE;                   /*       set loop flag */
        }
    } while (lpflg);                            /* until a NULL value */
    return (copxstr(&buf[0], nbyt));            /* alloc, return copy */
}

bool putxstr(                       /* returns TRUE on success */
    stream outp,                    /* stdio stream to write to */
    void *xspt,                     /* ptr to ext. string to copy */
    int nbyt)                       /* bytes per value */
{
    int n;
    char c, *cp = (char *) xspt;

    n = nbyt * xstrlen(xspt, nbyt);             /* get length in bytes */
    while (--n >= 0) {                          /* loop over bytes */
        c = *cp++;				/*   get byte to output */
        putc(0xff & c, outp);			/*   and write it out */
        if (ferror(outp))			/*   did output fail? */
            return (FALSE);                     /*     then so did we */
    }
    return (TRUE);                              /* return sign of success */
}

void *copxstr(                 /* returns ptr to new copy */
    void *xspt,                /* ptr to ext. string to copy */
    int nbyt)                  /* bytes per value */
{
    int n;
    char *dest, *dp, *cp = (char *) xspt;

    n = nbyt * xstrlen(xspt, nbyt);             /* get length in bytes */
    dp = dest = (char *) allocate(n);           /* allocate new storage */
    while (--n >= 0)                            /* loop over bytes */
        *dp++ = *cp++;                          /*   copy each in turn */
    return (dest);                              /* return copy string */
}

int xstrlen(                 /* returns count of values (w. NULL) */
    void *xspt,              /* ptr to ext. string to copy */
    int nbyt)                /* bytes per value */
{
    int nval, i;
    bool lpflg;
    char *cp = (char *) xspt;

    nval = 0;                                   /* init count of values */
    do {                                        /* loop over values */
        nval++;                                 /*   count one more */
        lpflg = FALSE;                          /*   init loop flag */
        for (i = 0; i < nbyt; i++)              /*   loop over bytes */
            if (*cp++ != 0)                     /*     a byte of data? */
                lpflg = TRUE;                   /*       set loop flag */
    } while (lpflg);                            /* until a NULL value */
    return (nval);                              /* return total count */
}

bool xstreq(             /* returns TRUE if equal */
    void *xp1,
    void *xp2,           /* ptrs to ext. strings to test */
    int nbyt)            /* bytes per value */
{
    int n;
    char *cp1 = (char *) xp1, *cp2 = (char *)xp2;

    n = nbyt * xstrlen(xp1, nbyt);              /* get length in bytes */
    while (--n >= 0)                            /* loop over bytes */
        if (*cp1++ != *cp2++)                   /*   bytes not equal? */
            return FALSE;                       /*     then strs unequal */
    return TRUE;                                /* indicate equality */
}

#ifdef TESTBED

main(int argc, char **argv)
{
    register int i;
    long lstr[32], *lcop;
    FILE *opt, *ipt;

    for (i = 0; i < 32; i++)
        lstr[i] = (i < 21 ? 12345 + 512 * i : 0);
    printf("xstrlen(lstr, %d) == %d\n",
           sizeof(long), xstrlen(lstr, sizeof(long)));
    lcop = (long *) copxstr(lstr, sizeof(long));
    printf("xstrlen(lcop, %d) == %d\n",
           sizeof(long), xstrlen(lcop, sizeof(long)));
    printf("xstreq(lstr, lcop, %d) == %d\n",
           sizeof(long), xstreq(lstr, lcop, sizeof(long)));
    printf("changing *lcop\n");
    *lcop = -1;
    printf("xstreq(lstr, lcop, %d) == %d\n",
           sizeof(long), xstreq(lstr, lcop, sizeof(long)));
    opt = stropen("foobar.dat", "w!");
    printf("putxstr(opt, lstr, %d) == %d\n",
           sizeof(long), putxstr(opt, lstr, sizeof(long)));
    printf("putxstr(opt, lcop, %d) == %d\n",
           sizeof(long), putxstr(opt, lcop, sizeof(long)));
    fclose(opt);
    ipt = stropen("foobar.dat", "r");
    lcop = (long *) getxstr(ipt, sizeof(long));
    if (lcop == NULL)
        printf("getxstr(ipt, %d) failed\n", sizeof(long));
    printf("xstreq(lstr, lcop, %d) == %d\n",
           sizeof(long), xstreq(lstr, lcop, sizeof(long)));
}

#endif
