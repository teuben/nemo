#include "vogl.h"

#define S_PIX_ASPECT	1.2 
static	int			old_mode = 3;
extern	unsigned int	_cur_color;
extern	unsigned int	_buffer_segment;

static	unsigned char 	pal[17] = {0, 4, 2, 14, 1, 5, 3, 15,
                        12, 10, 6, 9, 11, 13, 14, 15, 0};

extern	void	
		sig_line(),
		sig_set_colors();

extern	int	sigmaclear(),
		pc_fill(),
		pc_font(),
		pc_getkey(),
		pc_checkkey(),
		pc_locator(),
		pc_string(),
		setmode();

static	int
noop()
{
	return (-1);
}


static	unsigned	int	loc_val;

int
sigma_init()
{

	vdevice.sizeX = 399 * S_PIX_ASPECT;
	vdevice.sizeY = 399;
	vdevice.sizeSx = 639;
	vdevice.sizeSy = 399;
	vdevice.depth = 4;
	_buffer_segment = (unsigned)0xB800;
	old_mode = setmode(0x42);
	sigma_set_colors(pal);
	set_loc(64);
	pc_locinit(vdevice.sizeSx, vdevice.sizeSy);
	return (1);
}


/* 
 * sigma_vclear
 *
 *	Just clears the current viewport.
 */
sigma_vclear()
{
	int     x[4], y[4];

	if (vdevice.maxVx != vdevice.sizeSx
		|| vdevice.maxVy != vdevice.sizeSy
		|| vdevice.minVx != vdevice.sizeSx
		|| vdevice.minVy != vdevice.sizeSy) {
		x[0] = x[3] = vdevice.minVx;
		y[0] = y[1] = vdevice.maxVy;
		y[2] = y[3] = vdevice.minVy;
		x[1] = x[2] = vdevice.maxVx;

		pc_fill(5, x, y);
	} else {
		sigmaclear();
	}

	return(0);
}

/*
 * sigma_exit
 *
 *	Sets the display back to text mode.
 */
sigma_exit()
{
	unshowmouse();
	(void)setmode(3);
	return (1);
}

sigma_draw(x, y)
	int	x, y;
{
	sig_line(vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy, x, vdevice.sizeSy - y, _cur_color);
	vdevice.cpVx = x;
	vdevice.cpVy = y;

	return(0);
}

sigma_char(c)
	int	c;
{
	sigmachar(c, vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy, _cur_color);

	return(0);
}

sigma_color(i)
	int	i;
{
	_cur_color = (unsigned)i;

	return(0);
}

static DevEntry sigmadev = {
	"sigma",
	"large",
	"small",
	noop,
	sigma_char,
	pc_checkkey,
	sigma_vclear,
	sigma_color,
	sigma_draw,
	sigma_exit,
	pc_fill,
	pc_font,
	noop,
	pc_getkey,
	sigma_init,
	pc_locator,
	noop,
	noop,
	noop,
	pc_string,
	noop,
	noop
};

/*
 * _sigma_devcpy
 *
 *	copy the pc device into vdevice.dev.
 */
_sigma_devcpy()
{
	vdevice.dev = sigmadev;

	return(0);
}

