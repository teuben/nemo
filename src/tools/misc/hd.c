/***    hd - hex dump
 *
 *      Revision history:
 *
 *          Version 2.0:  Modified by Donald L. Nash to run on IBM PCs,
 *                        generic MS-DOS machines, generic UN!X machines
 *                        with big-endian architectures, and non-ASCII
 *                        machines which support the UN!X standard I/O
 *                        library.  Also modified to make dumps fit on
 *                        an 80 column screen and printer for machines
 *                        like an IBM PC which have 80 column printers.
 *
 *          Version 1.0:  Written by Kenneth J. Montgomery for the DEC
 *                        VAX running UN!X 4.2 BSD.  Originally made
 *                        dumps wider than 80 columns, indended for
 *                        line printers only.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "hd.h"

CHAR *titlep;                       /* pointer to name of current file */
long addr = 0l;                     /* current address in file */
int linel;                          /* current line on page */
CHAR encode[]="0123456789ABCDEF";   /* to convert int to hex number */

main(argc, argv)
int argc;
CHAR *argv[];

{

    FILE *fp, *nextfile();

    while ( (fp = nextfile(argc, argv)) != NULL )   /* while more files */
        while (dumpline(fp))        /* while more data */
            if (!--linel) page();   /* if at end of page, page eject */
    exit(0);

} /* main */


address()

/*
    This function prints out the current offset into the file.  If the
    offset is a multiple of 512 (0x200), then a '+' precedes the address.
    Otherwise a ' ' precedes it.
*/

{

    CHAR buf[10];                   /* to put address string in */
    register CHAR *bp;              /* to step through buf */
    register long a;                /* to manipulate current address */

    a = addr;                       /* init a to current address */
    bp = buf + 9;                   /* init bp to end of buff */
    *bp-- = '\0';                   /* null terminate buf */

    if (a % 0x200 == 0)             /* if on a 512 byte boundary */
         buf[0] = '+';              /* precede address by a '+' */
    else
         buf[0] = ' ';              /* else by a ' ' */

    while (bp >= buf + 1) {         /* while not at buf[1] */
        *bp-- = encode[a & 0xf];    /* convert first nibble in a to hex */
        a >>= 4;                    /* shift to get next nibble */
    } /* while */

    printf(" %s", buf);             /* print address */

} /* address */


dumpline(fp)
FILE *fp;

/*
    This function dumps a CBSIZE byte line to stdout.  Each line contains
    the current offset in the file for the start of the line, the hex digits
    for the CBSIZE bytes, and the ASCII representation of the CBSIZE bytes.
*/

{

    CHAR cbuf[CBSIZE];              /* input buffer (data from file) */
    CHAR fbuf[FBSIZE];              /* formatted data buffer */
    register int ccount;

    ccount = getbuf(cbuf, fp);      /* get CBSIZE bytes from fp to cbuf */
    if (!ccount) return(0);         /* if ccount == 0 no bytes read, so done */
    address();                      /* print address */
    format(cbuf, ccount, fbuf);     /* format ccount bytes from cbuf to fbuf */
    printf("%s    ", fbuf);         /* print formatted output */
    xlate(cbuf, ccount, fbuf);      /* xlate ccount bytes from cbuf to fbuf */
    printf("%s\n", fbuf);           /* print translated output */
    addr += ccount;                 /* update current address in file */
    return(ccount == CBSIZE);       /* if we read < CBSIZE bytes, return 0 */

} /* dumpline */


format(cbuf, ccount, fbuf)
CHAR *cbuf;
int ccount;
CHAR *fbuf;

/*
    This function reads ccount bytes from cbuf and turns them into
    hexadecimal digit pairs in fbuf.
*/

{

    register CHAR c, *src, *dst, *lmt;

    for (dst = fbuf, lmt = dst + FBSIZE; dst < lmt;)
        *dst++ = ' ';               /* blank fill output string */
    lmt = cbuf + ccount;

#if BIGENDIAN

/*
    If this is a big-endian machine, then dump from left to right, so start
    at the beginning of the output buffer.
*/

    dst = fbuf;
    dst[FBSIZE - 1] = '\0';         /* terminate string */
#else

/*
    If this is a little-endian machine, then dump from right to left, so
    start at the end of the output buffer.
*/

    dst = fbuf + FBSIZE - 1;
    *dst-- = '\0';                  /* terminate string */
#endif

    for (src = cbuf; src < lmt; src++) {    /* step through source buffer */
        c = *src;

#if BIGENDIAN
        *dst++ = ' ';                       /* put a space before hex number */
        *dst++ = encode[(c >> 4) & 0xf];    /* put in upper nibble */
        *dst++ = encode[c & 0xf];           /* put in lower nibble */
#else
        *dst-- = encode[c & 0xf];           /* put in lower nibble */
        *dst-- = encode[(c >> 4) & 0xf];    /* put in upper nibble */
        *dst-- = ' ';                       /* put a space after hex number */
#endif

    } /* for */


} /* format */


int getbuf(bp, fp)
CHAR *bp;
FILE *fp;

/*
    This function gets the next CBSIZE bytes from fp.  The bytes read are
    put in the buffer pointed to by bp.  The actual number of bytes read
    is returned.
*/

{

    register CHAR *abp;
    register int i, x;

    for (i = 0, abp = bp; i < CBSIZE; i++) {    /* read at most CBSIZE bytes */
        if ((x = getc(fp)) == EOF) break;       /* if EOF then done */
        *abp++ = x;                             /* put byte in buffer */
    } /* for */

    return(i);                      /* i counts # of bytes read */

} /* getbuf */


FILE *nextfile(argc, argv)
int argc;
CHAR *argv[];

/*
    This function steps through the files on the command line opening
    each file in turn.  It returns the FILE * for the file it opens.
!!!!	this function appeared a day later as a bug fix to original source
*/

{

    static int fn = 0;              /* keeps track of which argv[] is in use */
    static FILE *filep = (FILE *)NULL;

    addr = 0l;                      /* reset addr to zero */
    if (argc < 2) {                 /* if no args, use stdin */
        if (++fn == 1) {            /* only do this part once */
            titlep = (CHAR *)NULL;  /* set titlep */
            page();                 /* start a new page */
            return(stdin);          /* return stdin */
        } /* if */
        else
            return(NULL);           /* no more, so return NULL */
    } /* if */

loop:

    if (++fn < argc) {              /* if still more args to go */
        if (filep != NULL)          /* close any open file */
            fclose(filep);

#if MSDOS
        filep = fopen(argv[fn], "rb");  /* open next file */
#else
        filep = fopen(argv[fn], "r");   /* open next file */
#endif

        titlep = argv[fn];          /* set titlep */
        page();                     /* start a new page */

        if (filep == NULL) {        /* check for bad filename */
            printf("Cannot open file %s\n", argv[fn]);
            goto loop;
        } /* if */

        return(filep);              /* return file pointer */
    } /* if */

    return(NULL);                   /* no more files, return NULL */

} /* nextfile */


page()

/*
    This function starts a new page, prints the name of the input file or
    "Standard Input" if there are no files on the command line, and prints
    out the byte indecies for the line.
*/

{

    putc('\014', stdout);           /* page eject */

    if (titlep)
        printf("File:  %s", titlep); /* print filename */
    else
        printf("Standard Input");

/*
    Print index for hex part of dump.  Go left to right for big-endian,
    right to left for little-endian.
*/

#if BIGENDIAN
    printf("\n\n           +0 +1 +2 +3 +4 +5 +6 +7 +8 +9 +A +B +C +D +E +F");
#else
    printf("\n\n           +F +E +D +C +B +A +9 +8 +7 +6 +5 +4 +3 +2 +1 +0");
#endif

/*
    Print index for ASCII part of dump.  It is always left to right.
*/

    printf("    0123456789ABCDEF\n\n");
    linel = 55;                     /* set # of lines left on this page */

} /* page */


xlate(cbuf, ccount, fbuf)
CHAR *cbuf;
int ccount;
CHAR *fbuf;

/*
    This function reads ccount bytes from cbuf and puts them in fbuf.  A
    byte is copied into fbuf only if it is a printable character.  If it
    is not a printable character, a '.' is put in its place in fbuf.
*/

{

    register CHAR c, *src, *dst, *lmt;

    src = cbuf;                     /* pointer to source buffer */
    dst = fbuf;                     /* pointer to destination buffer */
    lmt = src + ccount;             /* limit of source buffer */

    while (src < lmt) {             /* while more bytes in source buffer */
        c = *src++;                 /* get next byte */
        if ( isprint(c) )           /* if it is printable */
            *dst++ = c;             /* then put it in destination buffer */
        else
            *dst++ = '.';           /* else put a '.' in destination buffer */
    } /* while */

    *dst++ = '\0';                  /* terminate the string */

} /* xlate */
