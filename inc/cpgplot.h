#ifndef cpgplot_h
#define cpgplot_h

#ifdef __cplusplus
extern "C" {
#endif

typedef int Logical;

void cpgarro(float x1, float y1, float x2, float y2);
void cpgask(Logical flag);
void cpgaxis(const char *opt, float x1, float y1, float x2, float y2, float v1, float v2, float step, int nsub, float dmajl, float dmajr, float fmin, float disp, float orient);
int cpgband(int mode, int posn, float xref, float yref, float *x, float *y, char *ch_scalar);
void cpgbbuf(void);
int cpgbeg(int unit, const char *file, int nxsub, int nysub);
void cpgbin(int nbin, const float *x, const float *data, Logical center);
void cpgbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub);
void cpgcirc(float xcent, float ycent, float radius);
void cpgclos(void);
void cpgconb(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr, float blank);
void cpgconf(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float c1, float c2, const float *tr);
void cpgconl(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float c, const float *tr, const char *label, int intval, int minint);
void cpgcons(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr);
void cpgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr);
void cpgctab(const float *l, const float *r, const float *g, const float *b, int nc, float contra, float bright);
int cpgcurs(float *x, float *y, char *ch_scalar);
void cpgdraw(float x, float y);
void cpgebuf(void);
void cpgend(void);
void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
void cpgeras(void);
void cpgerr1(int dir, float x, float y, float e, float t);
void cpgerrb(int dir, int n, const float *x, const float *y, const float *e, float t);
void cpgerrx(int n, const float *x1, const float *x2, const float *y, float t);
void cpgerry(int n, const float *x, const float *y1, const float *y2, float t);
void cpgetxt(void);
void cpggray(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float fg, float bg, const float *tr);
void cpghi2d(const float *data, int nxv, int nyv, int ix1, int ix2, int iy1, int iy2, const float *x, int ioff, float bias, Logical center, float *ylims);
void cpghist(int n, const float *data, float datmin, float datmax, int nbin, int pgflag);
void cpgiden(void);
void cpgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
void cpglcur(int maxpt, int *npt, float *x, float *y);
void cpgldev(void);
void cpglen(int units, const char *string, float *xl, float *yl);
void cpgline(int n, const float *xpts, const float *ypts);
void cpgmove(float x, float y);
void cpgmtxt(const char *side, float disp, float coord, float fjust, const char *text);
void cpgncur(int maxpt, int *npt, float *x, float *y, int symbol);
void cpgnumb(int mm, int pp, int form, char *string, int *string_length);
void cpgolin(int maxpt, int *npt, float *x, float *y, int symbol);
int cpgopen(const char *device);
void cpgpage(void);
void cpgpanl(int nxc, int nyc);
void cpgpap(float width, float aspect);
void cpgpixl(const int *ia, int idim, int jdim, int i1, int i2, int j1, int j2, float x1, float x2, float y1, float y2);
void cpgpnts(int n, const float *x, const float *y, const int *symbol, int ns);
void cpgpoly(int n, const float *xpts, const float *ypts);
void cpgpt(int n, const float *xpts, const float *ypts, int symbol);
void cpgpt1(float xpt, float ypt, int symbol);
void cpgptxt(float x, float y, float angle, float fjust, const char *text);
void cpgqah(int *fs, float *angle, float *barb);
void cpgqcf(int *font);
void cpgqch(float *size);
void cpgqci(int *ci);
void cpgqcir(int *icilo, int *icihi);
void cpgqclp(int *state);
void cpgqcol(int *ci1, int *ci2);
void cpgqcr(int ci, float *cr, float *cg, float *cb);
void cpgqcs(int units, float *xch, float *ych);
void cpgqdt(int n, char *type, int *type_length, char *descr, int *descr_length, int *inter);
void cpgqfs(int *fs);
void cpgqhs(float *angle, float *sepn, float *phase);
void cpgqid(int *id);
void cpgqinf(const char *item, char *value, int *value_length);
void cpgqitf(int *itf);
void cpgqls(int *ls);
void cpgqlw(int *lw);
void cpgqndt(int *n);
void cpgqpos(float *x, float *y);
void cpgqtbg(int *tbci);
void cpgqtxt(float x, float y, float angle, float fjust, const char *text, float *xbox, float *ybox);
void cpgqvp(int units, float *x1, float *x2, float *y1, float *y2);
void cpgqvsz(int units, float *x1, float *x2, float *y1, float *y2);
void cpgqwin(float *x1, float *x2, float *y1, float *y2);
void cpgrect(float x1, float x2, float y1, float y2);
float cpgrnd(float x, int *nsub);
void cpgrnge(float x1, float x2, float *xlo, float *xhi);
void cpgsah(int fs, float angle, float barb);
void cpgsave(void);
void cpgunsa(void);
void cpgscf(int font);
void cpgsch(float size);
void cpgsci(int ci);
void cpgscir(int icilo, int icihi);
void cpgsclp(int state);
void cpgscr(int ci, float cr, float cg, float cb);
void cpgscrl(float dx, float dy);
void cpgscrn(int ci, const char *name, int *ier);
void cpgsfs(int fs);
void cpgshls(int ci, float ch, float cl, float cs);
void cpgshs(float angle, float sepn, float phase);
void cpgsitf(int itf);
void cpgslct(int id);
void cpgsls(int ls);
void cpgslw(int lw);
void cpgstbg(int tbci);
void cpgsubp(int nxsub, int nysub);
void cpgsvp(float xleft, float xright, float ybot, float ytop);
void cpgswin(float x1, float x2, float y1, float y2);
void cpgtbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub);
void cpgtext(float x, float y, const char *text);
void cpgtick(float x1, float y1, float x2, float y2, float v, float tikl, float tikr, float disp, float orient, const char *str);
void cpgupdt(void);
void cpgvect(const float *a, const float *b, int idim, int jdim, int i1, int i2, int j1, int j2, float c, int nc, const float *tr, float blank);
void cpgvsiz(float xleft, float xright, float ybot, float ytop);
void cpgvstd(void);
void cpgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
void cpgwnad(float x1, float x2, float y1, float y2);

#ifdef __cplusplus
}
#endif

#endif
