#include "vogl.h"

/*
 * v4f_
 */
void
v4f_(vec)
	float	vec[4];
{
	vcall(vec, 4);
}

/*
 * v3f_
 */
void
v3f_(vec)
	float	vec[3];
{
	vcall(vec, 3);
}

/*
 * v2f_
 */
void
v2f_(vec)
	float	vec[2];
{
	vcall(vec, 2);
}

/*
 * v4d_
 */
void
v4d_(vec)
	double	vec[4];
{
	float	v[4];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	v[3] = vec[3];

	vcall(v, 4);
}

/*
 * v3d_
 */
void
v3d_(vec)
	double	vec[3];
{
	float	v[3];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	
	vcall(v, 3);
}

/*
 * v2d_
 */
void
v2d_(vec)
	double	vec[2];
{
	float	v[2];

	v[0] = vec[0];
	v[1] = vec[1];

	vcall(v, 2);
}

/*
 * v4s_
 */
void
v4s_(vec)
	short	vec[4];
{
	float	v[4];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	v[3] = vec[3];

	vcall(v, 4);
}

/*
 * v3s_
 */
void
v3s_(vec)
	short	vec[3];
{
	float	v[3];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	
	vcall(v, 3);
}

/*
 * v2s_
 */
void
v2s_(vec)
	short	vec[2];
{
	float	v[2];

	v[0] = vec[0];
	v[1] = vec[1];

	vcall(v, 2);
}


/*
 * v4i_
 */
void
v4i_(vec)
	int	vec[4];
{
	float	v[4];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	v[3] = vec[3];

	vcall(v, 4);
}

/*
 * v3i_
 */
void
v3i_(vec)
	int	vec[3];
{
	float	v[3];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	
	vcall(v, 3);
}

/*
 * v2i_
 */
void
v2i_(vec)
	int	vec[2];
{
	float	v[2];

	v[0] = vec[0];
	v[1] = vec[1];

	vcall(v, 2);
}

