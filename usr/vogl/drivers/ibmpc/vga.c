#include	"vogl.h"

/* 
 * Note, vga can't do double buffering by page swapping...`
 */

#define	V_PIX_ASPECT	1.0

static	int	old_mode;
extern	unsigned	int	_buffer_segment;
extern	unsigned	int	_buffer_offset;

extern	int
		vega_char(),
		vega_clear(),
		vega_color(),
		vega_draw(),
		vega_setpal(),
		pc_fill(),
		vega_font(),
		pc_getkey(),
		pc_checkkey(),
		pc_locator(),
		vega_string(),
		setmode();


static int
vga_init()
{
	old_mode = setmode(18);
	vdevice.sizeX = 479 * V_PIX_ASPECT;
	vdevice.sizeY = 479;
	vdevice.sizeSx = 639;
	vdevice.sizeSy = 479;
	vdevice.depth = 4;
	_buffer_segment = 0xA000;
	_buffer_offset = 0;
	pc_locinit(vdevice.sizeSx, vdevice.sizeSy);
	zsetup();
	vega_setpal();
	return (1);
}


/* 
 * vga_vclear
 *
 *	Just clears the current viewport.
 */
static
vga_vclear()
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
 * vga_exit
 *
 *	Sets the display back to text mode.
 */
static
vga_exit()
{
	(void)setmode(old_mode);
	return (1);
}

static	int
noop()
{
	return (-1);
}

static DevEntry vgadev = {
	"vga",
	"large",
	"small",
	noop,
	vega_char,
	pc_checkkey,
	vga_vclear,
	vega_color,
	vega_draw,
	vga_exit,
	pc_fill,
	vega_font,
	noop,
	pc_getkey,
	vga_init,
	pc_locator,
	noop,
	noop,
	noop,
	vega_string,
	noop,
	noop
};

/*
 * _vga_devcpy
 *
 *	copy the pc device into vdevice.dev.
 */
_vga_devcpy()
{
	vdevice.dev = vgadev;

	return(0);
}

