
/*
 *   src/misc/lsq.c
 */

void lsq_zero(int n, real *mat, real *vec);
void lsq_accum(int n, real *mat, real *vec, real *a, real w);
void lsq_solve(int n, real *mat, real *vec, real *sol);
void lsq_cfill(int n, real * mat, int c, real *vec);
