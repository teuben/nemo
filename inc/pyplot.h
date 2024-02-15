/*
 *   simple python script template generator
 */

stream pyplot_init(string fname);
void  pyplot_close(stream str);
void  pyplot_plot(stream str, string tabname, int *xcol, int *ycol, int *dycol, real *xrange, real *yrange);
void  pyplot_plot2(stream str, string tabname1, string tabname2, int *xcol, int *ycol, real *xrange, real *yrange);
void  pyplot_hist(stream str, string tabname, int *xcol, real *xrange, int bins);
