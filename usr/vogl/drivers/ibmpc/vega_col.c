#include <dos.h>

union REGS inregs, outregs;

static	unsigned		pal[17] = {0, 4, 2, 14, 1, 5, 3, 7,
                               12, 10, 6, 9, 11, 13, 14, 15, 0};

extern	unsigned	int	_cur_color;

vega_color(i)
	int	i;
{
	_cur_color = (unsigned)i;

	return(i);
}

vega_setpal()
{
	unsigned	i;

	for (i = 0; i < 16; i++) {
		inregs.h.ah = 0x10;
		inregs.h.al = 0;
		inregs.h.bl = i;
		inregs.h.bh = pal[i];
		int86(0x10, &inregs, &outregs);
	}

	return(1);
}
