#include "vogl.h"

/*
 * polarview_
 */
polarview_(dist, azim, inc, twist)
	float	*dist;
	int	*azim, *inc, *twist;
{
	polarview(*dist, (Angle)*azim, (Angle)*inc, (Angle)*twist);
}

/*
 * polarv_	(same as polarview_)
 */
polarv_(dist, azim, inc, twist)
	float	*dist;
	int	*azim, *inc, *twist;
{
	polarview(*dist, (Angle)*azim, (Angle)*inc, (Angle)*twist);
}

/*
 * lookat_
 */
lookat_(vx, vy, vz, px, py, pz, twist)
	float  *vx, *vy, *vz, *px, *py, *pz;
	int	*twist;
{
	lookat(*vx, *vy, *vz, *px, *py, *pz, (Angle)*twist);
}

/*
 * perspective_
 */
perspective_(fov, aspect, hither, yon)
	float 	*aspect, *hither, *yon;
	int	*fov;
{
	perspective((Angle)*fov, *aspect, *hither, *yon);
}

/*
 * perspe_	(same as perspective_)
 */
perspe_(fov, aspect, hither, yon)
	float 	*aspect, *hither, *yon;
	int	*fov;
{
	perspective((Angle)*fov, *aspect, *hither, *yon);
}

/*
 * window_
 */
window_(left, right, bottom, top, hither, yon)
	float 	*left, *right, *bottom, *top, *hither, *yon;
{
	window(*left, *right, *bottom, *top, *hither, *yon);
}

/*
 * ortho_
 */
ortho_(left, right, bottom, top, hither, yon)
	float 	*left, *right, *bottom, *top, *hither, *yon;
{
	ortho(*left, *right, *bottom, *top, *hither, *yon);
}

/*
 * ortho2_
 */
ortho2_(left, right, bottom, top)
	float	*left, *right, *bottom, *top;
{
	ortho2(*left, *right, *bottom, *top);
}
