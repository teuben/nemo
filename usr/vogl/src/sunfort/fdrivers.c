#include <stdio.h>
#include "vogl.h"

/*
 * voutput_
 */
void
voutput_(path, len, len0)
	char	*path;
	int	*len;
	int	len0;
{
	char		buf[BUFSIZ];
	register char	*p;

	strncpy(buf, path, *len);
	buf[*len] = 0;

	for (p = &buf[*len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	voutput(buf);
}

/*
 * voutpu_
 */
void
voutpu_(path, len, len0)
	char	*path;
	int	*len;
	int	len0;
{
	voutput_(path, len, len0);
}

/*
 * vinit_
 */
void
vinit_(dev, len, len0)
	char	*dev;
	int	*len;
	int	len0;
{
	char		buf[BUFSIZ];
	register char	*p;

	strncpy(buf, dev, *len);
	buf[*len] = 0;

	for (p = &buf[*len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	vinit(buf);
}

/*
 * vnewdev_
 */
void
vnewdev_(dev, len, len0)
	char	*dev;
	int	*len;
	int	len0;
{
	char		buf[BUFSIZ];
	register char	*p;

	strncpy(buf, dev, *len);
	buf[*len] = 0;

	for (p = &buf[*len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	vnewdev(buf);
}

/*
 * vnewde_
 */
void
vnewde_(dev, len, len0)
	char	*dev;
	int	*len;
	int	len0;
{
	vnewdev_(dev, len, len0);
}

/*
 * pushdev_
 */
void
pushdev_(dev, len)
	char	*dev;
	int	len;
{
	char		buf[BUFSIZ];
	register char	*p;

	strncpy(buf, dev, len);
	buf[len] = 0;

	for (p = &buf[len - 1]; *p == ' '; p--)
		;

	*++p = 0;

	pushdev(buf);
}

/*
 * 6 character version...
 */
void
pushde_(dev, len)
	char	*dev;
	int	len;
{
	pushdev_(dev, len);
}
	
/*
 * popdev_
 */
void
popdev_()
{
	popdev();
}

/*
 * winopen_
 */
int
winopen_(title, len, len0)
	char	*title;
	int	*len;
	int	len0;
{
	char		buf[BUFSIZ];

	strncpy(buf, title, len0);
	buf[*len] = 0;

	return(winopen(buf));
}
	
/*
 * winope_	(same as winopen_)
 */
int
winope_(title, len1, len2)
	char	*title;
	int	*len1, len2;
{
	char		buf[BUFSIZ];

	strncpy(buf, title, len2);
	buf[*len1] = 0;

	return(winopen(buf));
}

/*
 * gexit_
 */
void
gexit_()
{
	gexit();
}

/*
 * clear_
 */
void
clear_()
{
	clear();
}

/*
 * color_
 */
void
color_(col)
	int	*col;
{
	color(*col);
}

/*
 * colorf_
 */
void
colorf_(col)
	float	*col;
{
	colorf(*col);
}

/*
 * _mapcolor
 */
void
mapcolor_(indx, red, green, blue)
	int	*indx, *red, *green, *blue;
{
	mapcolor((Colorindex)*indx, (short)*red, (short)*green, (short)*blue);
}

/*
 * _mapcol	(same as mapcolor)
 */
void
mapcol_(indx, red, green, blue)
	int	*indx, *red, *green, *blue;
{
	mapcolor((Colorindex)*indx, (short)*red, (short)*green, (short)*blue);
}

/*
 * getplanes_
 */
int
getplanes_()
{
	return((int)getplanes());
}

/*
 * getpla_	(same as getplanes_)
 */
int
getpla_()
{
	return((int)getplanes());
}

/*
 * getvaluator_
 */
int
getvaluator_(dev)
	int	*dev;
{
	return((int)getvaluator((Device)*dev));
}

/*
 * getval_	(same as getvaluator_)
 */
int
getval_(dev)
	int	*dev;
{
	return((int)getvaluator((Device)*dev));
}

/*
 * getbutton_
 */
int
getbutton_(dev)
	int	*dev;
{
	return((int)getvaluator((Device)*dev));
}

/*
 * getbut_	(same as getbutton_)
 */
int
getbut_(dev)
	int	*dev;
{
	return((int)getbutton((Device)*dev));
}

/*
 * gconfig_	
 */
int
gconfig_()
{ }

/*
 * gconfi_	(same as gconfig_)
 */
int
gconfi_()
{ }

/*
 * reshapeviewport_	
 */
int
reshapeviewport_()
{ 
	reshapeviewport();
}

/*
 * reshap_	(same as reshapeviewport_)
 */
int
reshap_()
{
	reshapeviewport();
}

/*
 * winconstraints_	
 */
void
winconstraints_()
{ }

/*
 * wincon_	(same as winconstraints_)
 */
void
wincon_()
{ }

/*
 * shademodel_
 */
int
shademodel_(m)
	int	*m;
{ }

/*
 * shadem_	(same as shademodel_)
 */
int
shadem_(m)
	int	*m;
{ }

/*
 * getgdesc_
 */
int
getgdesc_(m)
	int	*m;
{
	return((int)getgdesc((long)*m));
}



void
getori_(ox, oy)
	int     *ox, *oy;
{
	getorigin(ox, oy);
}

void
getorigin_(ox, oy)
	int     *ox, *oy;
{
	getorigin(ox, oy);
}

/*
 * getgde_	(same as getgde_)
 */
int
getgde_(m)
	int	*m;
{
	return((int)getgdesc((long)*m));
}

/*
 * foreground_
 */
void
foreground_()
{
}

/*
 * foregr_
 */
void
foregr_()
{
}

/*
 * vsetflush_
 */
void
vsetflush_(yn)
	int	*yn;
{
	vdevice.sync = *yn;
}

/*
 * vsetfl_
 */
void
vsetfl_(yn)
	int	*yn;
{
	vdevice.sync = *yn;
}

/*
 * vflush_
 */
void
vflush_()
{
	(*vdevice.dev.Vsync)();
}
