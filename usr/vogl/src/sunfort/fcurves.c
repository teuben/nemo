#include "vogl.h"

/*
 * curvebasis_
 */
void
curvebasis_(basisid)
	int	*basisid;
{
	curvebasis((short)*basisid);
}

/*
 * curveb_	(Same as curvebasis_)
 */
void
curveb_(basisid)
	int	*basisid;
{
	curvebasis((short)*basisid);
}

/*
 * curveprecision_
 */
void
curveprecision_(nsegments)
	int	*nsegments;
{
	curveprecision((short)*nsegments);
}

/*
 * curvep_	(Same as curveprecision_)
 */
void
curvep_(nsegments)
	int	*nsegments;
{
	curveprecision((short)*nsegments);
}

/*
 * rcrv_
 */
void
rcrv_(geom)
	float	geom[4][4];
{
	rcrv(geom);
}

/*
 * rcrvn_
 */
void
rcrvn_(n, geom)
	int	*n;
	float	geom[][3];
{
	rcrvn((long)*n, geom);
}
/*
 * crv_
 */
void
crv_(geom)
	float	geom[4][3];
{
	crv(geom);
}

/*
 * crvn_
 */
void
crvn_(n, geom)
	int	*n;
	float	geom[][3];
{
	crvn((long)*n, geom);
}

/*
 * curveit_
 */
void
curveit_(n)
	int	*n;
{
	curveit((short)*n);
}

/*
 * curvei_
 */
void
curvei_(n)
	int	*n;
{
	curveit((short)*n);
}
