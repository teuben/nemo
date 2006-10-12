/*
 * SNAPKINPLOT.C: plot various diagnostics from snapkinem.

 15-09-2006   WD  deleted redundant and incorrect declartion of sqrt()
 12-10-2006  PJT  main -> nemo_main, and some prototyping

 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in=???\n       Input file (snapshot)",
    "VERSION=1.0\n  12-oct-06 PJT",
    NULL,
};

string usage = "plot various diagnostics from snapkinem";

string cvsid="$Id$";

real xrange[] = {  0.0,  1.0 };
real jrange[] = { -1.0,  1.0 };
real qrange[] = {  0.0,  1.0 };
real arange[] = {  0.0,  1.0 };

real xticks[] = { 0.0, 0.2, 0.4, 0.6, 0.8 };
real jticks[] = { -1.0, -0.5, 0.0, 0.5, 1.0 };
real qticks[] = { 0.0, 0.5, 1.0 };
real aticks[] = { 0.0, 0.5, 1.0 };

bool read_data(stream instr);
void plot_data(void);
real xtrans(real x);
real jtrans(real j);
real qtrans(real q);
real atrans(real a);
int plvect(real x1, real y1, real x2, real y2);

nemo_main()
{
    stream instr;
    bool read_data();

    instr = stropen(getparam("in"), "r");
    plinit("", 0.0, 20.0, 0.0, 20.0);
/*  pltext(getparam("in"), 2.0, 18.6, 0.32, 0.0);  */
    xaxis( 2.0, 16.0, 16.0, &xticks[1], 4, xtrans, "");
    yaxis( 2.0, 14.0,  4.0,  jticks,    5, jtrans, "");
    xaxis( 2.0,  8.0, 16.0,  xticks,    5, xtrans, "");
    yaxis( 2.0,  8.0,  4.0,  qticks,    3, qtrans, "");
    xaxis( 2.0,  2.0, 16.0,  xticks,    5, xtrans, "");
    yaxis( 2.0,  2.0,  4.0,  aticks,    3, atrans, "");
    plltype(0, 2);
    plvect(2.0, 18.0, 18.0, 18.0);
    plvect(2.0, 14.0, 18.0, 14.0);
    plvect(2.0, 12.0, 18.0, 12.0);
    plvect(2.0,  6.0, 18.0,  6.0);
    plltype(2, 1);
    while (read_data(instr))
	plot_data();
    plstop();
}

real w_sum;
real jz_norm;
real q_maj, q_int, q_min;
real qx_pro, qy_pro, qz_pro;

bool read_data(stream instr)
{
    real tkin, j_norm;

    if (fscanf(instr,
	       "%*[ \t\n\f] time: %lf n_tot: %*d n_sum: %*d w_sum: %lf",
	       &tkin, &w_sum) != 2)
	return (FALSE);
    if (fscanf(instr, " length x y z") != 0)
	return (FALSE);
    if (fscanf(instr, " pos: %*f %*f %*f %*f") != 0)
	return (FALSE);
    if (fscanf(instr, " vel: %*f %*f %*f %*f") != 0)
	return (FALSE);
    if (fscanf(instr, " jvec: %lf %*f %*f %lf", &j_norm, &jz_norm) != 2)
	return (FALSE);
    jz_norm = jz_norm / j_norm;
    if (fscanf(instr, " eigval x y z") != 0)
	return (FALSE);
    if (fscanf(instr, " qpole: %lf %lf %*f %*f", &q_maj, &qx_pro) != 2)
	return (FALSE);
    if (fscanf(instr, "        %lf %*f %lf %*f", &q_int, &qy_pro) != 2)
	return (FALSE);
    if (fscanf(instr, "        %lf %*f %*f %lf", &q_min, &qz_pro) != 2)
	return (FALSE);
    if (fscanf(instr, " keten: %*f %*f %*f %*f") != 0)
	return (FALSE);
    if (fscanf(instr, "        %*f %*f %*f %*f") != 0)
	return (FALSE);
    if (fscanf(instr, "        %*f %*f %*f %*f") != 0)
	return (FALSE);
    return (TRUE);
}

bool first_time = TRUE;
real w_tot = 0.0;

real x0, jz0, qx0, qy0, qz0, ba0, ca0;

void plot_data(void)
{
    real x, jz, qx, qy, qz, ba, ca;

    w_tot = w_tot + w_sum;
    x  = xtrans(w_tot / 2.5);
    jz = jtrans(jz_norm);
    qx = qtrans(ABS(qx_pro));
    qy = qtrans(ABS(qy_pro));
    qz = qtrans(ABS(qz_pro));
    ba = atrans(sqrt(q_int / q_maj));
    ca = atrans(sqrt(q_min / q_maj));
    plcircle(x, jz, 0.08);
    plcircle(x, qx, 0.08);
    plbox(   x, qy, -0.08);
    plbox(   x, qz, 0.08);
    plbox(   x, ba, -0.08);
    plbox(   x, ca, 0.08);
    if (! first_time) {
	plvect(x0, jz0, x, jz);
	plvect(x0, qx0, x, qx);
	plvect(x0, qy0, x, qy);
	plvect(x0, qz0, x, qz);
	plvect(x0, ba0, x, ba);
	plvect(x0, ca0, x, ca);
    }
    x0 = x;
    jz0 = jz;
    qx0 = qx;
    qy0 = qy;
    qz0 = qz;
    ba0 = ba;
    ca0 = ca;
    first_time = FALSE;
}

real xtrans(real x)
{
    return ( 2.0 + 16.0 * (x - xrange[0]) / (xrange[1] - xrange[0]));
}

real jtrans(real j)
{
    return (14.0 +  4.0 * (j - jrange[0]) / (jrange[1] - jrange[0]));
}

real qtrans(real q)
{
    return ( 8.0 +  4.0 * (q - qrange[0]) / (qrange[1] - qrange[0]));
}

real atrans(real a)
{
    return ( 2.0 +  4.0 * (a - arange[0]) / (arange[1] - arange[0]));
}

plvect(real x1, real y1, real x2, real y2)
{
    plmove(x1, y1);
    plline(x2, y2);
}
