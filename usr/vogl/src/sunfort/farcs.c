#include "vogl.h"

/*
 * circleprecision_
 */
void
circleprecision_(prec)
	int	*prec;
{
	circleprecision(*prec);
}

/*
 * circpr_
 */
void
circpr_(prec)
	int	*prec;
{
	circleprecision(*prec);
}

/*
 * arcprecision_
 */
void
arcprecision_(prec)
	int	*prec;
{
	circleprecision(*prec);
}

/*
 * arcpre_
 */
void
arcpre_(prec)
	int	*prec;
{
	circleprecision(*prec);
}
/*
 * arc_
 */
void
arc_(x, y, radius, startang, endang)
	Coord	*x, *y, *radius;
	int	*startang, *endang;
{
	arc(*x, *y, *radius, (Angle)*startang, (Angle)*endang);
}

/*
 * arci_
 */
void
arci_(x, y, radius, startang, endang)
	int	*x, *y, *radius;
	int	*startang, *endang;
{
	arci((Icoord)*x, (Icoord)*y, *radius, (Angle)*startang, (Angle)*endang);
}

/*
 * arcs_
 */
void
arcs_(x, y, radius, startang, endang)
	short	*x, *y, *radius;
	short	*startang, *endang;
{
	arcs(*x, *y, *radius, (Angle)*startang, (Angle)*endang);
}

/*
 * arcf_
 */
void
arcf_(x, y, radius, startang, endang)
	Coord	*x, *y, *radius;
	int	*startang, *endang;
{
	arcf(*x, *y, *radius, (Angle)*startang, (Angle)*endang);
}

/*
 * arcfi_
 */
void
arcfi_(x, y, radius, startang, endang)
	int	*x, *y, *radius;
	int	*startang, *endang;
{
	arcfi((Icoord)*x, (Icoord)*y, (Icoord)*radius, (Angle)*startang, (Angle)*endang);
}

/*
 * arcfs_
 */
void
arcfs_(x, y, radius, startang, endang)
	short	*x, *y, *radius;
	short	*startang, *endang;
{
	arcfs(*x, *y, *radius, (Angle)*startang, (Angle)*endang);
}

/*
 * circ_
 */
void
circ_(x, y, radius)
	float	*x, *y, *radius;
{
	circ(*x, *y, *radius);
}

/*
 * circi_
 */
void
circi_(x, y, radius)
	int	*x, *y, *radius;
{
	circi((Icoord)*x, (Icoord)*y, (Icoord)*radius);
}

/*
 * circs_
 */
void
circs_(x, y, radius)
	short	*x, *y, *radius;
{
	circs(*x, *y, *radius);
}

/*
 * circf_
 */
void
circf_(x, y, radius)
	float	*x, *y, *radius;
{
	circf(*x, *y, *radius);
}

/*
 * circfi_
 */
void
circfi_(x, y, radius)
	int	*x, *y, *radius;
{
	circfi((Icoord)*x, (Icoord)*y, (Icoord)*radius);
}

/*
 * circfs_
 */
void
circfs_(x, y, radius)
	short	*x, *y, *radius;
{
	circfs(*x, *y, *radius);
}
