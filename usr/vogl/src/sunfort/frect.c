#include "vogl.h"

/*
 * rect_
 */
rect_(x1, y1, x2, y2)
	float 	*x1, *y1, *x2, *y2;
{
	rect(*x1, *y1, *x2, *y2);
}

/*
 * rects_
 */
rects_(x1, y1, x2, y2)
	short 	*x1, *y1, *x2, *y2;
{
	rect((float)*x1, (float)*y1, (float)*x2, (float)*y2);
}

/*
 * recti_
 */
recti_(x1, y1, x2, y2)
	int 	*x1, *y1, *x2, *y2;
{
	rect((float)*x1, (float)*y1, (float)*x2, (float)*y2);
}

/*
 * rectf_
 */
rectf_(x1, y1, x2, y2)
	float 	*x1, *y1, *x2, *y2;
{
	rectf(*x1, *y1, *x2, *y2);
}

/*
 * rectfs_
 */
rectfs_(x1, y1, x2, y2)
	short 	*x1, *y1, *x2, *y2;
{
	rectf((float)*x1, (float)*y1, (float)*x2, (float)*y2);
}

/*
 * rectfi_
 */
rectfi_(x1, y1, x2, y2)
	int 	*x1, *y1, *x2, *y2;
{
	rectf((float)*x1, (float)*y1, (float)*x2, (float)*y2);
}
