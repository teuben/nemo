#include "vogl.h"

static	int	have_mouse = 0;

extern	int	ismouse(), readmouse();

pc_locinit(x, y)
	int	x, y;
{
	if ((have_mouse = ismouse(x, y))) 
		showmouse();

	return(0);
}

int
pc_locator(x, y)
	int *x, *y;
{
	int ix, iy, b;

	if (!have_mouse) 
		return (-1);


	b = readmouse(&ix, &iy);

	*x = ix;
	*y = vdevice.sizeSy - iy;
	return (b);
}
