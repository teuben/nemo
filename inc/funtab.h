/*
 * FunctionTable:
 */

/* lookup modes */

#define FUNTAB_LINEAR  0x01
#define FUNTAB_SPLINE  0x02
#define FUNTAB_NEAR    0x03

typedef struct FunctionTable {
    string name;        /* ID or filename */
    int mode;           /* lookup mode (one of the above FUNTAB_xxx */
    int n;              /* Number of points in table */
    real *x;		/* pointer to array of X values */
    real *y;		/* pointer to array of Y values */
    real *coeff;        /* (spline) coefficients, if used */
    int errors;         /* cumulative errors in interpolation */
} FunctionTable;


/* funtab.c */
FunctionTable *ft_open   (string, int, int, int);
real           ft_spline (FunctionTable *, real);
real           ft_linear (FunctionTable *, real);
int            ft_close  (FunctionTable *);

