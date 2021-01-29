//#define PG6

//            a few "long" (>6) names have a short equivalent in the library
#ifdef PG6
#define  pgbegin_  pgbeg_
#define  pgcurse_  pgcurs_
#define  pglabel_  pglab_
#define  pgmtext_  pgmtxt_
#define  pgncurse_ pgncur_
#define  pgpaper_  pgpap_
#define  pgpoint_  pgpt_
#define  pgptext_  pgptxt_
#define  pgvport_  pgsvp_
#define  pgvsize_  pgvsiz_
#define  pgvstand_ pgvstd_
#define  pgwindow_ pgswin_
#endif

//  yapp_pglot and yapp_giza only use the following pgplot functions
//  these are the BSD style fortran names, with underscored and int 

int  pgbegin_(int *unit, char *file, int *nxsub, int *nysub, int filelen);
void pgask_(int *flag);
void pgqvsz_(int *units, float *x1, float *x2, float *y1, float *y2);
void pgsvp_(float *xleft, float *xright, float *ybot, float *ytop);
void pgpaper_(float *width, float *aspect);
void pgwindow_(float *x1, float *x2, float *y1, float *y2);
void pgslw_(int *lw);
void pgsls_(int *ls);
int  pgcurse_(float *x, float *y, char *ch_scalar, int chlen);
void pgpoint_(int *n, float *xpts, float *ypts, int *symbol);
void pgptext_(float *x, float *y, float *angle, float *fjust, char *text, int text_len);
void pgswindow_(float *x1, float *x2, float *y1, float *y2);
void pgdraw_(float *x, float *y);
void pgmove_(float *x, float *y);
void pgsch_(float *size);
void pgqinf_(char *item, char *value, int *value_length, int itemlen, int valuelen);
void pgupdt_();
void pgpage_();
int  pgiden_();   
int  pgend_();
void pggray_(float *a, int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2, float *fg, float *bg, float *tr);
void pgcont_(float *a, int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2, float *c, int *nc, float *tr);
void pgsci_(int *ci);
void pgscr_(int *ci, float *red, float *green, float *blue);


//   see cpgplot.h for any new ones that may be used later
