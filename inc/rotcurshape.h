#include <nemo.h>

typedef real (*rcproc)(real, int, real *, real *);

/* see: $NEMO/src/image/rotcur/rotcurs.c */

extern real rotcur_flat    (real r, int n, real *p, real *d);
extern real rotcur_linear  (real r, int n, real *p, real *d);
extern real rotcur_poly    (real r, int n, real *p, real *d);
extern real rotcur_core    (real r, int n, real *p, real *d);
extern real rotcur_core1   (real r, int n, real *p, real *d);
extern real rotcur_core2   (real r, int n, real *p, real *d);
extern real rotcur_plummer (real r, int n, real *p, real *d);
extern real rotcur_tanh    (real r, int n, real *p, real *d);
extern real rotcur_iso     (real r, int n, real *p, real *d);
extern real rotcur_exp     (real r, int n, real *p, real *d);
extern real rotcur_nfw     (real r, int n, real *p, real *d);
extern real rotcur_moore   (real r, int n, real *p, real *d);
extern real rotcur_brandt  (real r, int n, real *p, real *d);
extern real rotcur_power   (real r, int n, real *p, real *d);
extern real rotcur_disk1   (real r, int n, real *p, real *d);
