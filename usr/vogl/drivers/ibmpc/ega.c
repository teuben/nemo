#include	"vogl.h"

#define	E_PIX_ASPECT	1.3

static	int	old_mode;
extern	unsigned	int	_buffer_segment;
extern	unsigned	int	_buffer_offset;

extern	int
		vega_backbuf(),
		vega_char(),
		vega_clear(),
		vega_color(),
		vega_draw(),
		vega_setpal(),
		pc_fill(),
		vega_font(),
		vega_frontbuf(),
		pc_getkey(),
		pc_checkkey(),
		pc_locator(),
		vega_string(),
		vega_swapbuf(),
		setmode();

static int
ega_init()
{
	old_mode = setmode(16);
	vdevice.sizeX = 349 * E_PIX_ASPECT;
	vdevice.sizeY = 349;
	vdevice.sizeSx = 639;
	vdevice.sizeSy = 349;
	vdevice.depth = 4;
	_buffer_segment = (unsigned) 0xA000;
	_buffer_offset = 0;
	pc_locinit(vdevice.sizeSx, vdevice.sizeSy);
	zsetup();
	vega_setpal();
	return (1);
}

/* 
 * ega_vclear
 *
 *	Just clears the current viewport.
 */
static
ega_vclear()
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
		vega_clear();
	}

	return(0);
}


/*
 * ega_exit
 *
 *	Sets the display back to text mode.
 */
static
ega_exit()
{
	(void)setmode(old_mode);
	return (1);
}

static	int
noop()
{
	return (-1);
}

static DevEntry egadev = {
	"ega",
	"large",
	"small",
	vega_backbuf,
	vega_char,
	pc_checkkey,
	ega_vclear,
	vega_color,
	vega_draw,
	ega_exit,
	pc_fill,
	vega_font,
	vega_frontbuf,
	pc_getkey,
	ega_init,
	pc_locator,
	noop,
	noop,
	noop,
	vega_string,
	vega_swapbuf,
	noop
};

/*
 * _ega_devcpy
 *
 *	copy the pc device into vdevice.dev.
 */
_ega_devcpy()
{
	vdevice.dev = egadev;

	return(0);
}

