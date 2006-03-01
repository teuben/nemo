
/*
 * Fortran and C I/O interface
 *	some of this file depends on the Fortran_to_C interface
 */



/* f2c */
extern void nb4open_  (string, int);
extern void nb4put_   (int*, float*, float*, float*);
extern void nb4get_   (int*, float*, float*, float*);
extern void nb4close_ (void);

/* f2c */
extern void nb3open_  (string, int);
extern void nb3header_(int*, int*, int*, int*);
extern void nb3data_  (int*, int*, float*, float*, float*, float*, int*);
extern void nb3close_ (void);

/* C */
extern void nb3open_c  (string, int, bool);
extern void nb3header_c(int*, int*, int*, int*);
extern void nb3data_c  (int*, int*, int *, float*, float*, float*, float*, float *, float *, int*, int*);
extern void nb3close_c (void);

