/*
 * EXTSTRING.H: an extension of the standard C concept of character
 * strings to strings of n-byte (nonzero) values, terminated by a
 * marker of n zero bytes.
 *
 * These should not be confused with unicode type extensions to strings,
 * but more general.
 *
 *  20-sep-93       ANSI header file created                pjt
 *
 */

#if defined(__cplusplus)
extern "C" {
#endif

extern void *getxstr (stream inpt, int nbyt);
extern bool putxstr  (stream outp, void *xspt, int nbyt);
extern void *copxstr (void *xspt, int nbyt);
extern int  xstrlen  (void *xspt, int nbyt);
extern bool xstreq   (void *xp1, void *xp2, int nbyt);

#if defined(__cplusplus)
}
#endif

