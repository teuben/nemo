/*
 * Sets BIOS Video modes.
 *
 * This file is also a convienient place to declare some Global
 * variables that get used in the *.asm Routines.
 */

unsigned	int	_buffer_segment = 0xA000;	/* Current buffer being drawn to */
unsigned	int	_buffer_offset = 0;	/* Current buffer being drawn to */
unsigned	int	_cur_mode = 3;	/* The current video mode */
unsigned	int	_cur_color = 0;	/* The current color */

#include <dos.h>

union	REGS	regs;


int
setmode(mode)
	int	mode;
{
	int	old_mode;

	_cur_mode = mode;

	regs.h.ah = 15;
	int86(0x10, &regs, &regs);
	old_mode = regs.h.al;
	regs.h.ah = 0;
	regs.h.al = mode;
	int86(0x10, &regs, &regs);
	return (old_mode);
}
