#include "vogl.h"

/*
 * prefposition_
 */
void
prefposition_(a, b, c, d)
	int	*a, *b, *c, *d;
{
	prefposition((long)*a, (long)*b, (long)*c, (long)*d);
}

/*
 * prefpo_	(same as prefposition_)
 */
void
prefpo_(a, b, c, d)
	int	*a, *b, *c, *d;
{
	prefposition((long)*a, (long)*b, (long)*c, (long)*d);
}

/*
 * prefsize_
 */
void
prefsize_(x, y)
	int	*x, *y;
{
	prefsize((long)*x, (long)*y);
}

/*
 * prefsi_	(same as prefsize_)
 */
void
prefsi_(x, y)
	int	*x, *y;
{
	prefsize((long)*x, (long)*y);
}
