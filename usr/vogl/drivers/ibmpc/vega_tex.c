#include	"vogl.h"


#define		SMALL_CHARS		0
#define		LARGE_CHARS		1

static		int		vega_size = LARGE_CHARS;
extern		unsigned	_cur_color;

void
vega_string(str)
	char	*str;
{
	char	*p;
	register	int		x, y, w;

	x = vdevice.cpVx;
	y = vdevice.sizeSy - vdevice.cpVy - (int) vdevice.hheight;
	w = (int)vdevice.hwidth;
	for (p = str; *p; p++) {
		vegachar(x, y, *p, _cur_color, vega_size);
		x += w;
		}
	vdevice.cpVx = x;
}

void
vega_char(c)
	char	c;
{
	vegachar(vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy - (int)vdevice.hheight, c, _cur_color, vega_size);
	vdevice.cpVx += (int)vdevice.hwidth;
}

vegachar(x, y, c, n, size)
int			x, y, c, n, size;
{
	if (size == LARGE_CHARS) 
		egalchar(x, y, c, n);
	else 
		egaschar(x, y, c, n);

	return(0);
}

vega_font(name)
	char *name;
{
	if (strcmp(name, "small") == 0) {
		vega_size = SMALL_CHARS;
		vdevice.hwidth = 8.0;
		vdevice.hheight = 8.0;
	} else if (strcmp(name, "large") == 0) {
		vega_size = LARGE_CHARS;
		vdevice.hwidth = 8.0;
		vdevice.hheight = 14.0;
	} else {
		vega_size = LARGE_CHARS;
		vdevice.hwidth = 8.0;
		vdevice.hheight = 14.0;
		return (0);
	}
	return (1);
}
