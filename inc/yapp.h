/*
 *  YAPP: Yet Another Plotting Package:
 *
 *	definitions of accessible functions in most YAPP's
 */

#if !defined(_yapp_inc)
#define _yapp_inc

#ifdef __cplusplus
extern "C" {
#endif
                                                /* standard functions: */
extern int plinit        (string,real,real,real,real);
extern int plltype       (int,int);
extern int plline        (real,real);
extern int plmove        (real,real);
extern int plpoint       (real,real);
extern int plcircle      (real,real,real);
extern int plcross       (real,real,real);
extern int plbox         (real,real,real);
extern int pltext        (string,real,real,real,real);
extern int pljust        (int);
extern int plflush       (void);
extern int plframe       (void);
extern int plstop        (void);

                                                /* hardly ever used ? */
extern int plswap        (void);
extern real plxscale     (real,real);
extern real plyscale     (real,real);

                                        	/* -DCOLOR or perhaps stubbed */
extern void plcolor      (int);
extern int  plncolor     (void);
extern void plpalette    (real*,real*,real*,int);
                                                /* non_standard */ 
extern int pl_interp     (string);
extern int pl_screendump (string);
extern int pl_frame      (void);
extern int pl_getpoly    (float *,float *,int);
extern int pl_cursor     (real *, real *, char *);
extern int pl_matrix     (real *,int,int,real,real,real,real,real,real,real);

#ifdef __cplusplus
}
#endif


#endif
