#include <stdio.h>
#include "vogl.h"

/*
 * font_
 */
void
font_(id)
	int	*id;
{
	font(*id);
}

/*
 * charst_	(same as charstr_)
 */
charst_(s, len, len2) 
	char	*s;
	int	*len, len2;
{
	charstr_(s, len, len2);
}

/*
 * charstr_
 */
charstr_(s, len, len2) 
	char	*s;
	int	*len, len2;
{

	/*
	 * *len is the user passed length, len2 is the compiler's idea
	 * of the length...
	 */
        char            buf[BUFSIZ];
	register char   *p;

	strncpy(buf, s, *len);
	buf[*len] = 0;

	for (p = &buf[*len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	charstr(buf);
}

/*
 * cmov_
 */
void
cmov_(x, y, z)
	float	*x, *y, *z;
{
	cmov(*x, *y, *z);
}

/*
 * cmovs_
 */
void
cmovs_(x, y, z)
	short	*x, *y, *z;
{
	cmov((float)*x, (float)*y, (float)*z);
}

/*
 * cmovi_
 */
void
cmovi_(x, y, z)
	int	*x, *y, *z;
{
	cmov((float)*x, (float)*y, (float)*z);
}

/*
 * cmov2_
 */
void
cmov2_(x, y)
	float	*x, *y;
{
	cmov(*x, *y, 0.0);
}

/*
 * cmov2s_
 */
void
cmov2s_(x, y)
	short	*x, *y;
{
	cmov((float)*x, (float)*y, 0.0);
}

/*
 * cmov2i_
 */
void
cmov2i_(x, y)
	int	*x, *y;
{
	cmov((float)*x, (float)*y, 0.0);
}

#ifdef OLD_GL
/*
 * getwidth_
 */
int
getwidth_()
{
	return((int)getwidth());
}

/*
 * getwid_
 */
int
getwid_()
{
	return((int)getwidth());
}

#endif

/*
 * getheight_
 */
int
getheight_()
{
	return((int)getheight());
}

/*
 * gethei_
 */
int
gethei_()
{
	return((int)getheight());
}

/*
 * strwid_
 *
 * Return the length of a string in pixels
 * Same as strwidth_
 */
strwid_(s, len, len2) 
	char	*s;
	int	*len, len2;
{
	strwidth_(s, len, len2);
}

/*
 * strwidth_
 *
 * Return the length of a string in pixels
 */
int
strwidth_(s, len, len2) 
	char	*s;
	int	*len, len2;
{
        char            buf[BUFSIZ];
	register char   *p;

	strncpy(buf, s, *len);
	buf[*len] = 0;

	for (p = &buf[*len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	return((int)strwidth(buf));
}


/*
 * getcpo_ and getcpos_
 *
 * Return the current chracater position in screen coords.
 */
void
getcpo_(ix, iy)
	short	*ix, *iy;
{
	getcpos(ix, iy);
}


void
getcpos_(ix, iy)
	short	*ix, *iy;
{
	getcpos(ix, iy);
}
