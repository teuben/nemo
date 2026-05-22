/* TABCOG:   integrate a sorted table using the Curve of Growth method
 *
 *
 *      19-may-2026  Q&D first version, cloned off tabint
 */

#include <stdinc.h> 
#include <getparam.h>
#include <spline.h>
#include <extstring.h>
#include <table.h>
#include <mdarray.h>

string defv[] = {
  "in=???\n         Input table file",
  "xcol=1\n         Column with X coordinate (must be sorted in that column)",
  "ycol=2\n         Column with Y coordinate of function",
  "xc=\n            Override central point",
  "xmin=\n          value if data below xmin to be discarded",
  "xmax=\n          value if data above xmax to be discarded",
  "tol=0.1\n        flat tolerance",
  "bench=1\n        How often to call curve_of_growth",
  "VERSION=0.2\n    20-may-2026 PJT",
  NULL,

};

string usage="integrate a function using curve of growth method";

/* ------------------------------------------------------------------ */
/* Result struct (mirrors the dict returned by curve_of_growth)       */
/* ------------------------------------------------------------------ */

#define MAX_WIDTH_FRAC 10

typedef struct {
    double flux,      flux_std;
    double flux_r,    flux_r_std;
    double flux_b,    flux_b_std;
    double A_F, A_C, C_V;
    double rms;
    int    bchan, echan;
    double vel,       vel_std;
    /* parallel arrays for each width fraction */
    int    n_wf;
    double wf[MAX_WIDTH_FRAC];
    double width[MAX_WIDTH_FRAC];
    double width_std[MAX_WIDTH_FRAC];
} CogResult;


extern int minmax(int, real *, real *, real *);
local void reverse(int n, real *x, real *y);
local void sortxy(int n, real *x, real *y);

local void curve_of_growth(const double *x, const double *y, int n,
			   double vc,
			   const double *wf_in, int n_wf_in,
			   int bchan, int echan,
			   double flat_tol, double fw,
			   CogResult *result);


void nemo_main()
{
  int colnr[2], i, n, nmax, nsteps;
  real *xdat, *ydat, xmin, xmax, ymin, ymax, zmin, zmax;
  real x, y, z, s, xold, yold, dx, dz, sum, sum0, *sdat;
  string spectrum = getparam("in");
  bool Qmin = hasvalue("xmin");
  bool Qmax = hasvalue("xmax");
  int nbench = getiparam("bench");
  int mom = 0;
  double xc = NAN;
  double flat_tol = getdparam("tol");
  
  /* read the data */

  tableptr t = table_open(stropen(spectrum,"r"),0);
  n = nmax = table_nrows(t);
  int ncols = table_ncols(t);
  dprintf(1,"%s has %d x %d table\n", spectrum, n, ncols);
  colnr[0] = getiparam("xcol");
  colnr[1] = getiparam("ycol");
  mdarray2 d = table_md2cr(t,2,colnr,0,0);
  xdat = &d[0][0];
  ydat = &d[1][0];

  if (hasvalue("xc"))  xc = getdparam("xc");

  /* reverse arrays if not sorted properly */
  if (xdat[0] > xdat[1])
    reverse(n, xdat, ydat);

  /* check if array is now properly sorted */
  int nbad = 0;
  for (i=1; i<n; i++)
    if (xdat[i] < xdat[i-1]) nbad++;
  if (nbad > 0) sortxy(n, xdat, ydat);
  
  /* figure out an appropriate min/max */
  minmax(n,xdat,&xmin,&xmax);
  minmax(n,ydat,&ymin,&ymax);
  dprintf(1,"X range: %g : %g\n",xmin,xmax);
  dprintf(1,"Y range: %g : %g\n",ymin,ymax);
  if (Qmin || Qmax) {
    if (Qmin) xmin = getrparam("xmin");
    if (Qmax) xmax = getrparam("xmax");
    dprintf(1,"X range: %g : %g   (reset)\n",xmin,xmax) ;
    int imin = n;
    int imax = 0;
    for (i=0; i<n; i++) {
      if (Qmin && xdat[i] < xmin) imin = i;
      if (Qmax && xdat[i] < xmax) imax = i;
    }
    dprintf(1,"New range %d - %d   %g - %g\n", imin, imax, xdat[imin], xdat[imax]);
    /* new data */
    n = imax - imin + 1;
    xdat = &xdat[imin];
    ydat = &ydat[imin];
  }

  CogResult res;
  while (nbench-- > 0) {
    //                            vc   wf_in n_wf_in, bchan,echan, flat_tol, fw
    curve_of_growth(xdat, ydat, n, xc, NULL, 0, -1, -1, flat_tol, 1.0, &res);
  }

  printf("flux     = %g\n", res.flux);
  printf("flux_std = %g\n", res.flux_std);
  printf("vel      = %g\n", res.vel);
  printf("vel_std  = %g\n", res.vel_std);
  printf("A_F      = %g\n", res.A_F);
  printf("A_C      = %g\n", res.A_C);
  printf("C_V      = %g\n", res.C_V);
  printf("rms      = %g\n", res.rms);
  printf("bchan    = %d\n", res.bchan);
  printf("echan    = %d\n", res.echan);
  for (int k = 0; k < res.n_wf; k++)
    printf("width[%.2f] = %g +/- %g\n",
	   res.wf[k], res.width[k], res.width_std[k]);

}  

local void reverse(int n, real *x, real *y)
{
  real t;
  for (int i = 0; i<n/2; i++) {
    t = x[i];  x[i] = x[n-i-1];  x[n-i-1] = t;
    t = y[i];  y[i] = y[n-i-1];  y[n-i-1] = t;
  }
}

local void sortxy(int n, real *x, real *y)
{
  error("Your data were not sorted, and sortxy is not implemented yet");
}


/*
 * Curve of Growth - C translation of cog.py
 *
 * Assumes clean input arrays (no NaN / masked values).
 */

/* ------------------------------------------------------------------ */
/* Utility functions                                                   */
/* ------------------------------------------------------------------ */

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

/* np.diff: out[i] = a[i+1]-a[i], length n-1 */
static void array_diff(const double *a, int n, double *out) {
    for (int i = 0; i < n - 1; i++)
        out[i] = a[i + 1] - a[i];
}

/* np.cumsum: out[i] = sum(a[0..i]) */
static void array_cumsum(const double *a, int n, double *out) {
    double s = 0.0;
    for (int i = 0; i < n; i++) { s += a[i]; out[i] = s; }
}

/* np.std (population std) of a[0..n-1] */
static double array_std(const double *a, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a[i];
    double mean = sum / n;
    double var = 0.0;
    for (int i = 0; i < n; i++) var += (a[i] - mean) * (a[i] - mean);
    return sqrt(var / n);
}

/* np.median of a[start..end-1] (sorts a copy) */
static double array_median(const double *a, int start, int end) {
    int n = end - start;
    if (n <= 0) return 0.0;
    double *tmp = (double *)malloc(n * sizeof(double));
    memcpy(tmp, a + start, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_double);
    double med = (n % 2 == 0) ? (tmp[n/2-1] + tmp[n/2]) / 2.0 : tmp[n/2];
    free(tmp);
    return med;
}

/* np.argmin(abs(a - val)) over a[0..n-1] */
static int argmin_abs(const double *a, int n, double val) {
    int idx = 0;
    double best = fabs(a[0] - val);
    for (int i = 1; i < n; i++) {
        double d = fabs(a[i] - val);
        if (d < best) { best = d; idx = i; }
    }
    return idx;
}

/* ------------------------------------------------------------------ */
/* cog_slope                                                           */
/* ------------------------------------------------------------------ */
/*
 * slope     : output, length n-1 (caller allocates)
 * slope_rms : output
 * flat_idx0 : output; first index where |slope| < flat_tol*rms, or -1
 */
static void cog_slope(const double *c, int n, double flat_tol,
                      double *slope, double *slope_rms, int *flat_idx0) {
    int m = n - 1;
    array_diff(c, n, slope);
    *slope_rms = array_std(slope, m);
    *flat_idx0 = -1;
    for (int i = 0; i < m; i++) {
        if (fabs(slope[i]) < flat_tol * (*slope_rms)) {
            *flat_idx0 = i;
            break;
        }
    }
}

/* ------------------------------------------------------------------ */
/* cog_flux                                                            */
/* ------------------------------------------------------------------ */
/*
 * flux      : np.median(c[flat_idx0:])
 * flux_std  : np.std(c[flat_idx0:])
 * slope_out : np.median(slope[:flat_idx0])
 *
 * Python negative-index behaviour when flat_idx0 == -1:
 *   c[-1:]      -> last element only  (flux_start = n-1)
 *   slope[:-1]  -> first n-2 elements (slope_end  = n-2)
 */
static void cog_flux(const double *c, int n, double flat_tol,
                     double *flux, double *flux_std, double *slope_out) {
    double *slope = (double *)malloc((n - 1) * sizeof(double));
    double slope_rms;
    int flat_idx0;
    cog_slope(c, n, flat_tol, slope, &slope_rms, &flat_idx0);

    int flux_start = (flat_idx0 < 0) ? n - 1  : flat_idx0;
    int slope_end  = (flat_idx0 < 0) ? n - 2  : flat_idx0;

    *flux      = array_median(c, flux_start, n);
    *flux_std  = array_std(c + flux_start, n - flux_start);
    *slope_out = array_median(slope, 0, slope_end);

    free(slope);
}

/* ------------------------------------------------------------------ */
/* curve_of_growth                                                     */
/* ------------------------------------------------------------------ */
/*
 * vc       : central velocity; pass NAN to auto-estimate (moment-1)
 * wf_in    : width fractions array; pass NULL for default
 *            {0.25, 0.65, 0.75, 0.85, 0.95}
 * n_wf_in  : length of wf_in (ignored when wf_in == NULL)
 * bchan    : first signal channel; pass -1 to auto-estimate
 * echan    : last  signal channel; pass -1 to auto-estimate
 * flat_tol : tolerance for flat region of the CoG
 * fw       : multiplier on max width when estimating line-free range
 */
void curve_of_growth(const double *x, const double *y, int n,
                     double vc,
                     const double *wf_in, int n_wf_in,
                     int bchan, int echan,
                     double flat_tol, double fw,
                     CogResult *result)
{
    /* ---- width fractions ---------------------------------------- */
    double wf[MAX_WIDTH_FRAC];
    int nwf;
    if (wf_in == NULL) {
        double def[] = {0.25, 0.65, 0.75, 0.85, 0.95};
        nwf = 5;
        memcpy(wf, def, nwf * sizeof(double));
    } else {
        nwf = (n_wf_in < MAX_WIDTH_FRAC) ? n_wf_in : MAX_WIDTH_FRAC;
        memcpy(wf, wf_in, nwf * sizeof(double));
    }
    /* insert 0.25 at front if missing, append 0.85 if missing */
    int has025 = 0, has085 = 0;
    for (int k = 0; k < nwf; k++) {
        if (wf[k] == 0.25) has025 = 1;
        if (wf[k] == 0.85) has085 = 1;
    }
    if (!has025 && nwf < MAX_WIDTH_FRAC) {
        memmove(wf + 1, wf, nwf * sizeof(double));
        wf[0] = 0.25; nwf++;
        has085 = 0;
        for (int k = 0; k < nwf; k++) if (wf[k] == 0.85) has085 = 1;
    }
    if (!has085 && nwf < MAX_WIDTH_FRAC)
        wf[nwf++] = 0.85;

    /* ---- sort x ascending: _x = x[::p], _y = y[::p] ----------- */
    int p = (x[0] <= x[1]) ? 1 : -1;
    double *_x = (double *)malloc(n * sizeof(double));
    double *_y = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        int src = (p == 1) ? i : (n - 1 - i);
        _x[i] = x[src];
        _y[i] = y[src];
    }

    /* ---- channel crop ------------------------------------------- */
    int i0 = 0, i1 = n;
    if (echan >= 0 && echan <= n) i1 = echan;
    if (bchan >= 0)               i0 = bchan;
    int nc = i1 - i0;
    double *cx = _x + i0;   /* views into sorted arrays */
    double *cy = _y + i0;

    /* dx[i] = cx[i+1]-cx[i],  ydx[i] = cy[i]*dx[i]  (length nc-1) */
    double *dx  = (double *)malloc((nc - 1) * sizeof(double));
    double *ydx = (double *)malloc((nc - 1) * sizeof(double));
    array_diff(cx, nc, dx);
    for (int i = 0; i < nc - 1; i++) ydx[i] = cy[i] * dx[i];

    /* ---- central velocity --------------------------------------- */
    double _vc;
    if (!isnan(vc)) {
        _vc = vc;
    } else {
        double num = 0.0, den = 0.0;
        for (int i = 0; i < nc; i++) { num += cx[i] * cy[i]; den += cy[i]; }
        _vc = num / den;
    }
    int vc_idx = argmin_abs(cx, nc, _vc);

    /* ---- blue CoG  b  = cumsum(ydx[:vc_idx][::-1]) ------------ */
    int lb = vc_idx;
    double *b = (double *)malloc((lb > 0 ? lb : 1) * sizeof(double));
    if (lb > 0) {
        double *rev = (double *)malloc(lb * sizeof(double));
        for (int i = 0; i < lb; i++) rev[i] = ydx[lb - 1 - i];
        array_cumsum(rev, lb, b);
        free(rev);
    }

    /* bp = cumsum(ydx[:vc_idx+1][::-1]) */
    int lbp = vc_idx + 1;
    double *bp = (double *)malloc(lbp * sizeof(double));
    {
        double *rev = (double *)malloc(lbp * sizeof(double));
        for (int i = 0; i < lbp; i++) rev[i] = ydx[lbp - 1 - i];
        array_cumsum(rev, lbp, bp);
        free(rev);
    }

    /* ---- red  CoG  r  = cumsum(ydx[vc_idx:]) ------------------ */
    int lr = nc - 1 - vc_idx;
    double *r = (double *)malloc((lr > 0 ? lr : 1) * sizeof(double));
    if (lr > 0) array_cumsum(ydx + vc_idx, lr, r);

    /* rp = cumsum(ydx[vc_idx+1:]) */
    int lrp = nc - 2 - vc_idx;
    double *rp = (double *)malloc((lrp > 0 ? lrp : 1) * sizeof(double));
    if (lrp > 0) array_cumsum(ydx + vc_idx + 1, lrp, rp);

    /* ---- combined  t = b[:s] + r[:s],  dx = dx[:s] ------------ */
    int s = (lb < lr) ? lb : lr;
    double *t   = (double *)malloc((s > 0 ? s : 1) * sizeof(double));
    double *dxt = (double *)malloc((s > 0 ? s : 1) * sizeof(double));
    for (int i = 0; i < s; i++) t[i] = b[i] + r[i];
    memcpy(dxt, dx, s * sizeof(double));

    /* ---- xt = xr[:s] + xb[:s] ---------------------------------- */
    /* xr[i] = cx[vc_idx+i] - _vc,  xb[i] = |cx[vc_idx-1-i] - _vc| */
    int lxr = nc - vc_idx;
    int lxb = vc_idx;
    int sx  = (lxr < lxb) ? lxr : lxb;
    if (sx > s) sx = s;
    double *xt = (double *)malloc((sx > 0 ? sx : 1) * sizeof(double));
    for (int i = 0; i < sx; i++)
        xt[i] = (cx[vc_idx + i] - _vc) + fabs(cx[vc_idx - 1 - i] - _vc);

    /* ---- flux --------------------------------------------------- */
    double flux, flux_std, _slope_dummy;
    cog_flux(t, s, flat_tol, &flux, &flux_std, &_slope_dummy);
    flux_std = sqrt(flux_std*flux_std + (1.03*flux_std)*(1.03*flux_std));

    double flux_b, flux_b_std, slope_b;
    cog_flux(b, lb > 0 ? lb : 1, flat_tol, &flux_b, &flux_b_std, &slope_b);
    double flux_bp, _d1, _d2;
    cog_flux(bp, lbp, flat_tol, &flux_bp, &_d1, &_d2);
    double dflux_b = flux_b - flux_bp;
    flux_b_std = sqrt(flux_b_std*flux_b_std
                    + (1.03*flux_b_std)*(1.03*flux_b_std)
                    + dflux_b*dflux_b);

    double flux_r, flux_r_std, slope_r;
    cog_flux(r, lr > 0 ? lr : 1, flat_tol, &flux_r, &flux_r_std, &slope_r);
    double flux_rp, _d3, _d4;
    if (lrp > 0)
        cog_flux(rp, lrp, flat_tol, &flux_rp, &_d3, &_d4);
    else
        flux_rp = flux_r;
    double dflux_r = flux_r - flux_rp;
    flux_r_std = sqrt(flux_r_std*flux_r_std
                    + (1.03*flux_r_std)*(1.03*flux_r_std)
                    + dflux_r*dflux_r);

    /* ---- line widths -------------------------------------------- */
    double *slope_t = (double *)malloc((s > 1 ? s - 1 : 1) * sizeof(double));
    double slope_rms_t;
    int flat_idx0_t;
    cog_slope(t, s, flat_tol, slope_t, &slope_rms_t, &flat_idx0_t);
    free(slope_t);

    /* nt = t[:flat_idx0] / flux  (same negative-index rule as cog_flux) */
    int fi = (flat_idx0_t < 0) ? (s - 1) : flat_idx0_t;
    if (fi < 1) fi = 1;

    double *nt = (double *)malloc(fi * sizeof(double));
    for (int i = 0; i < fi; i++) nt[i] = t[i] / flux;

    double widths[MAX_WIDTH_FRAC];
    for (int k = 0; k < nwf; k++) {
        int idx = argmin_abs(nt, fi, wf[k]);
        widths[k] = (idx < sx) ? xt[idx] : xt[sx - 1];
    }

    /* ---- rms from line-free channels ---------------------------- */
    double max_wf = wf[0];
    int    max_wf_k = 0;
    for (int k = 1; k < nwf; k++)
        if (wf[k] > max_wf) { max_wf = wf[k]; max_wf_k = k; }

    int _bchan, _echan;
    if (bchan < 0) {
        _bchan = argmin_abs(x, n, _vc - fw * widths[max_wf_k]);
        if (_bchan <= 0) _bchan = 0;
    } else {
        _bchan = bchan;
    }
    if (echan < 0) {
        _echan = argmin_abs(x, n, _vc + fw * widths[max_wf_k]);
        if (_echan >= n) _echan = n;
    } else {
        _echan = echan;
    }
    if (_bchan > _echan) { int tmp = _bchan; _bchan = _echan; _echan = tmp; }

    /* rms = std(hstack(y[:_bchan], y[_echan:])) on original y */
    int rms_n = _bchan + (n - _echan);
    double *rms_buf = (double *)malloc((rms_n > 0 ? rms_n : 1) * sizeof(double));
    memcpy(rms_buf,            y,           _bchan * sizeof(double));
    memcpy(rms_buf + _bchan,   y + _echan,  (n - _echan) * sizeof(double));
    double rms = (rms_n > 0) ? array_std(rms_buf, rms_n) : 0.0;
    free(rms_buf);

    /* ---- error on centroid ------------------------------------- */
    double vc_std = 0.0;
    if (isnan(vc)) {
        double sum_xy = 0.0, sum_y = 0.0, sum_x2 = 0.0;
        for (int i = 0; i < nc; i++) {
            sum_xy += cy[i] * cx[i];
            sum_y  += cy[i];
            sum_x2 += cx[i] * cx[i];
        }
        double fac1 = _vc / sum_xy * rms * sqrt(sum_x2);
        double fac2 = sqrt((double)nc) * rms * _vc / sum_y;
        vc_std = sqrt(fac1*fac1 + fac2*fac2);
    }

    /* ---- errors on widths -------------------------------------- */
    /* nt_std[i] = sqrt((nt[i]/t[i]*rms*dx[i])^2 + (nt[i]/flux*flux_std)^2) */
    double *nt_std = (double *)malloc(fi * sizeof(double));
    for (int i = 0; i < fi; i++) {
        double term1 = (fabs(t[i]) > 0.0) ? nt[i] / t[i] * rms * dxt[i] : 0.0;
        double term2 = nt[i] / flux * flux_std;
        nt_std[i] = sqrt(term1*term1 + term2*term2);
    }

    double widths_std[MAX_WIDTH_FRAC];
    for (int k = 0; k < nwf; k++) {
        int idx   = argmin_abs(nt, fi, wf[k]);
        double sv = nt_std[idx];
        /* idx_m: argmin|nt - sv - wf[k]|  -> nt nearest wf[k]+sv */
        int idx_m = argmin_abs(nt, fi, wf[k] + sv);
        /* idx_p: argmin|nt + sv - wf[k]|  -> nt nearest wf[k]-sv */
        int idx_p = argmin_abs(nt, fi, wf[k] - sv);
        double std_m = (idx_m < sx && idx < sx) ? xt[idx_m] - xt[idx] : 0.0;
        double std_p = (idx_p < sx && idx < sx) ? xt[idx]   - xt[idx_p] : 0.0;
        double dx_idx = (idx < s) ? dxt[idx] : 0.0;
        double ws = fmax(fmax(std_m, std_p), dx_idx);
        widths_std[k] = sqrt(ws*ws + (widths[k]*0.01)*(widths[k]*0.01));
    }

    /* ---- asymmetry & concentration ----------------------------- */
    double a_f = flux_b / flux_r;
    if (a_f < 1.0) a_f = 1.0 / a_f;

    double a_c = slope_b / slope_r;
    if (a_c < 1.0) a_c = 1.0 / a_c;

    double w025 = 0.0, w085 = 0.0;
    for (int k = 0; k < nwf; k++) {
        if (wf[k] == 0.25) w025 = widths[k];
        if (wf[k] == 0.85) w085 = widths[k];
    }
    double c_v = w085 / w025;

    /* ---- fill result ------------------------------------------- */
    result->flux       = flux;       result->flux_std   = flux_std;
    result->flux_r     = flux_r;     result->flux_r_std = flux_r_std;
    result->flux_b     = flux_b;     result->flux_b_std = flux_b_std;
    result->A_F        = a_f;
    result->A_C        = a_c;
    result->C_V        = c_v;
    result->rms        = rms;
    result->bchan      = _bchan;     result->echan      = _echan;
    result->vel        = _vc;        result->vel_std    = vc_std;
    result->n_wf       = nwf;
    memcpy(result->wf,        wf,         nwf * sizeof(double));
    memcpy(result->width,     widths,     nwf * sizeof(double));
    memcpy(result->width_std, widths_std, nwf * sizeof(double));

    /* ---- cleanup ----------------------------------------------- */
    free(_x); free(_y);
    free(dx); free(ydx);
    free(b); free(bp); free(r); free(rp);
    free(t); free(dxt); free(xt);
    free(nt); free(nt_std);
}
