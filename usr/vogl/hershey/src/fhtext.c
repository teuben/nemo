
#include <stdio.h>

/*
 * hsetpath_
 */
void
hsetpath_(fpath, len1, len2)
	char	*fpath;
	int	*len1;
	int	len2;
{
	char	buf[BUFSIZ];

	strncpy(buf, fpath, len2);
	buf[*len1] = 0;

	hsetpath(buf);
}

/*
 * hsetpa_
 */
void
hsetpa_(fpath, len1, len2)
	char	*fpath;
	int	*len1;
	int	len2;
{
	char	buf[BUFSIZ];

	strncpy(buf, fpath, len2);
	buf[*len1] = 0;

	hsetpath(buf);
}

/*
 * hfont_
 */
void
hfont_(fontfile, len1, len2)
	char	*fontfile;
	int	*len1;
	int	len2;
{
	char	buf[BUFSIZ];

	strncpy(buf, fontfile, len2);
	buf[*len1] = 0;

	hfont(buf);
}

/*
 * htextsize_
 */
void
htextsize_(width, height)
	float	*width, *height;
{
	htextsize(*width, *height);
}

/*
 * htexts_
 */
void
htexts_(width, height)
	float	*width, *height;
{
	htextsize(*width, *height);
}

/*
 * hboxtext_
 */
hboxtext_(x, y, l, h, s, length, len)
	float	*x, *y, *l, *h;
	char	*s;
	int	*length;
	int	len;
{
	char		buf[BUFSIZ];
	register char   *p;

	strncpy(buf, s, len);
	buf[*length] = 0;

	hboxtext(*x, *y, *l, *h, buf);
}

/*
 * hboxte_	(same as hboxtext_)
 */
hboxte_(x, y, l, h, s, length, len)
	float	*x, *y, *l, *h;
	char	*s;
	int	*length;
	int	len;
{
	char		buf[BUFSIZ];
	register char   *p;

	strncpy(buf, s, len);
	buf[*length] = 0;

	hboxtext(*x, *y, *l, *h, buf);
}

/*
 * hboxfit_
 */
hboxfit_(l, h, nchars)
	float	*l, *h;
	int	*nchars;
{
	hboxfit(*l, *h, *nchars);
}

/*
 * hboxfi_
 */
hboxfi_(l, h, nchars)
	float	*l, *h;
	int	*nchars;
{
	hboxfit(*l, *h, *nchars);
}

/*
 * htextang_
 */
void
htextang_(ang)
	float	*ang;
{
	htextang(*ang);
}

/*
 * htexta_
 */
void
htexta_(ang)
	float	*ang;
{
	htextang(*ang);
}

/*
 * hdrawchar_
 */
hdrawchar_(s)
	char	*s;
{
	hdrawchar(*s);
}

/*
 * hdrawc_
 */
hdrawc_(s)
	char	*s;
{
	hdrawchar(*s);
}

/*
 * hcharstr_
 */
hcharstr_(s, length, len)
	char	*s;
	int	*length;
	int	len;
{
        char            buf[BUFSIZ];
	register char   *p;

	strncpy(buf, s, len);
	buf[*length] = 0;
	hcharstr(buf);
}

/*
 * hchars_
 */
hchars_(s, length, len)
	char	*s;
	int	*length;
	int	len;
{
        char            buf[BUFSIZ];
	register char   *p;

	strncpy(buf, s, len);
	buf[*length] = 0;
	hcharstr(buf);
}

/*
 * hgetfontheight_
 */
float
hgetfontheight_()
{
	return(hgetfontheight());
}

/*
 * hgetfh_
 */
float
hgetfh_()
{
	return(hgetfontheight());
}

/*
 * hgetfontwidth_
 */
float
hgetfontwidth_()
{
	return(hgetfontwidth());
}

/*
 * hgetfw_
 */
float
hgetfw_()
{
	return(hgetfontwidth());
}

/*
 * hgetdecender_
 */
float
hgetdecender_()
{
	return(hgetdecender());
}

/*
 * hgetde_
 */
float
hgetde_()
{
	return(hgetdecender());
}

/*
 * hgetascender_
 */
float
hgetascender_()
{
	return(hgetascender());
}

/*
 * hgetas_
 */
float
hgetas_()
{
	return(hgetascender());
}

/*
 * hgetfontsize_
 */
void
hgetfontsize_(cw, ch)
	float 	*cw, *ch;
{
	hgetfontsize(cw, ch);
}

/*
 * hgetfs_
 */
void
hgetfs_(cw, ch)
	float 	*cw, *ch;
{
	hgetfontsize(cw, ch);
}

/*
 * hgetcharsize_
 */
void
hgetcharsize_(c, cw, ch)
	char	*c;
	float 	*cw, *ch;
{
	hgetcharsize(*c, cw, ch);
}

/*
 * hgetch_
 */
void
hgetch_(c, cw, ch)
	char	*c;
	float 	*cw, *ch;
{
	hgetcharsize(*c, cw, ch);
}

/*
 * hfixedwidth
 */
void
hfixedwidth_(i)
	int	*i;
{
	hfixedwidth(*i);
}

/*
 * hfixed_
 */
void
hfixed_(i)
	int	*i;
{
	hfixedwidth(*i);
}

/*
 * hcentertext
 */
void
hcentertext_(i)
	int	*i;
{
	hcentertext(*i);
}

/*
 * hcente_
 */
void
hcente_(i)
	int	*i;
{
	hcentertext(*i);
}

/*
 * hrightjustify_
 */
void
hrightjustify_(i)
	int	*i;
{
	hrightjustify(*i);
}

/*
 * hright_
 */
void
hright_(i)
	int	*i;
{
	hrightjustify(*i);
}

/*
 * hleftjustify_
 */
void
hleftjustify_(i)
	int	*i;
{
	hleftjustify(*i);
}

/*
 * hleftj_
 */
void
hleftj_(i)
	int	*i;
{
	hleftjustify(*i);
}

/*
 * hnumchars_
 */
int
hnumchars_()
{
	return(hnumchars());
}

/*
 * hnumch_
 */
int
hnumch_()
{
	return(hnumchars());
}

/*
 * hstrlength_
 */
float
hstrlength_(str, len0, len1)
	char	*str;
	int	*len0;
	int	len1;
{
        char            buf[BUFSIZ];
	register char   *p;

	strncpy(buf, str, len1);
	buf[*len0] = 0;

	return(hstrlength(buf));
}

/*
 * hstrle_
 */
float
hstrle_(str, len)
	char	*str;
	int	len;
{
        char            buf[BUFSIZ];
	register char   *p;

	strncpy(buf, str, len);
	buf[len] = 0;

	for (p = &buf[len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	return(hstrlength(buf));
}
