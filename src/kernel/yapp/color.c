/*
 * COLOR: stubs folr the (4) color functions:
 *		plcolor
 *		plncolor
 *		plpalette
 *		pllut
 *
 */

void plcolor(color)
int color;
{
}

/*
 * PLNCOLORS: return current value of local variable ncolors.
 */

int plncolors()
{
    return 0;
}

/*
 * PLPALETTE: re-initialize color table from user-supplied values.
 */

void plpalette(r, g, b, nc)
real r[], g[], b[];
int nc;
{
    warning("plpalette: Your yapp was not compiled/compilable with -DCOLOR");
}

void pllut(fname, compress)
string fname;
bool compress;
{
    warning("pllut: LUT %s not loaded; yapp not compiled with -DCOLOR",fname);
}
