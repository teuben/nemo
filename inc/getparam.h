/*
 * GETPARAM.H: command-line processing definitions
 *
 * Mar. 1987  Joshua Edward Barnes, Princeton NJ.
 * Sep. 1990  -- added finiparam() to list              PJT
 * Nov. 1991  -- added hasvalue()                       PJT
 * July 1993 -- C++ entrant
 * Mar  1994 -- added nemoinpf/r
 * Feb  1995 -- added updparam, no more ARGS
 * Jun  2001 -- added some convenient ZENO macro	PJT
 */

#ifndef _getparam_h
#define _getparam_h

#ifdef __cplusplus
extern "C" {
#endif

extern void    initparam (string *, string *);
extern void    finiparam (void);
extern void    stop (int);
extern void    putparam (string, string);
extern void    promptparam (string, string);
extern int     cntparam (void);
extern bool    isaparam (string);
extern bool    hasvalue (string);
extern bool    updparam (string);
extern int     indexparam(string, int);

extern string  getparam  (string);
extern int     getiparam (string);
extern long    getlparam (string);
extern bool    getbparam (string);
extern double  getdparam (string);

extern string  getparam_idx  (string, int);
extern int     getiparam_idx (string, int);
extern long    getlparam_idx (string, int);
extern bool    getbparam_idx (string, int);
extern double  getdparam_idx (string, int);

extern int nemoinpi (string, int *, int);
extern int nemoinpd (string, double *, int);
extern int nemoinpf (string, float *, int);
extern int nemoinpl (string, long *, int);
extern int nemoinpb (string, bool *, int);

extern int    natoi(string);
extern double natof(string);

#if !defined(SINGLEPREC)
# define nemoinpr  nemoinpd
# define getrparam getdparam
#else
# define nemoinpr  nemoinpf
# define getrparam getfparam
#endif

#ifdef __cplusplus
}
#endif

/*
 * Macro used to obtain name and version of program.
 */

#define getargv0()      (getparam("argv0"))
#define getversion()    (getparam("VERSION"))


#endif /* _getparam_h */
