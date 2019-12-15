
/*
 * simple multi-dimensional arrays  - sorry, no template programming
 * (real = double or float, depending on the setting in NEMO)
 *
 * MDMAXDIM: maximum number of dimensions we support here;
 * it's hardcoded.
 *
 * See: NEMO/src/kernel/misc/mdarray.c
 *
 */

#define MDMAXDIM   8

typedef real        *mdarray1;    /* v1[n1]                              */
typedef mdarray1    *mdarray2;    /* v2[n2][n1]                          */
typedef mdarray2    *mdarray3;    /* v3[n3][n2][n1]                      */
typedef mdarray3    *mdarray4;    /* v4[n4][n3][n2][n1]                  */
typedef mdarray4    *mdarray5;    /* v5[n5][n4][n3][n2][n1]              */
typedef mdarray5    *mdarray6;    /* v6[n6][n5][n4][n3][n2][n1]          */
typedef mdarray6    *mdarray7;    /* v7[n7][n6][n5][n4][n3][n2][n1]      */
typedef mdarray7    *mdarray8;    /* v8[n8][n7][n6][n5][n4][n3][n2][n1]  */

mdarray1 allocate_mdarray1(int n1);
mdarray2 allocate_mdarray2(int n2, int n1);
mdarray3 allocate_mdarray3(int n3, int n2, int n1);
mdarray4 allocate_mdarray4(int n4, int n3, int n2, int n1);
mdarray5 allocate_mdarray5(int n5, int n4, int n3, int n2, int n1);
mdarray6 allocate_mdarray6(int n6, int n5, int n4, int n3, int n2, int n1);
mdarray7 allocate_mdarray7(int n7, int n6, int n5, int n4, int n3, int n2, int n1);
mdarray8 allocate_mdarray8(int n8, int n7, int n6, int n5, int n4, int n3, int n2, int n1);

void free_mdarray1(mdarray1 x, int n1);
void free_mdarray2(mdarray2 x, int n2, int n1);
void free_mdarray3(mdarray3 x, int n3, int n2, int n1);
void free_mdarray4(mdarray4 x, int n4, int n3, int n2, int n1);
void free_mdarray5(mdarray5 x, int n5, int n4, int n3, int n2, int n1);
void free_mdarray6(mdarray6 x, int n6, int n5, int n4, int n3, int n2, int n1);
void free_mdarray7(mdarray7 x, int n7, int n6, int n5, int n4, int n3, int n2, int n1);
void free_mdarray8(mdarray8 x, int n8, int n7, int n6, int n5, int n4, int n3, int n2, int n1);

