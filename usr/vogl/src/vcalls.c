#include "vogl.h"

/*
 * Handles all the v4f v3f v2f calls
 */


/*
 *  vcall
 *
 *	Specify a generic point.
 */
void
vcall(vector, len)
	float	vector[];
	int	len;
{
	Vector	vec;
	static	Coord	tmesh[2][3];	/* 2 vectr register for triangular mesh */
        static	Coord	qstrip[3][3];	/* 3 vectr register for quad mesh */
        int	sjm;			/* temporary variable */

	vec[0] = vector[0];
	vec[1] = vector[1];
	vec[2] = 0.0;
	vec[3] = 1.0;

	if (len == 3) {
		vec[2] = vector[2];
		vec[3] = 1.0;
	} else if (len == 4) {
		vec[2] = vector[2];
		vec[3] = vector[3];
	}

	if (vdevice.save) {
		vdevice.savex = vec[V_X];
		vdevice.savey = vec[V_Y];
		vdevice.savez = vec[V_Z];
	} 
	
	switch (vdevice.bgnmode) {
	case VLINE:
		if (vdevice.save) {
			vdevice.save = 0;
			move(vec[V_X], vec[V_Y], vec[V_Z]);
			break;
		}

		draw(vec[V_X], vec[V_Y], vec[V_Z]);
		break;
	case VPNT:
		pnt(vec[V_X], vec[V_Y], vec[V_Z]);
		break;
	case VCLINE:
		if (vdevice.save) {
			vdevice.save = 0;
			move(vec[V_X], vec[V_Y], vec[V_Z]);
			break;
		}
			
		draw(vec[V_X], vec[V_Y], vec[V_Z]);
		break;
	case VPOLY:
		if (vdevice.save) {
			vdevice.save = 0;
			pmv(vec[V_X], vec[V_Y], vec[V_Z]);
			break;
		} 

		pdr(vec[V_X], vec[V_Y], vec[V_Z]);
		break;

	case VTMESH:
		if( vdevice.save < 2 ) {  /* save the first two calls */
			tmesh[vdevice.save][V_X] = vec[V_X];
			tmesh[vdevice.save][V_Y] = vec[V_Y];
			tmesh[vdevice.save][V_Z] = vec[V_Z];
		} else {  /* plot triangle */
			pmv(tmesh[0][V_X], tmesh[0][V_Y], tmesh[0][V_Z]);
			pdr(tmesh[1][V_X], tmesh[1][V_Y], tmesh[1][V_Z]);
			pdr(     vec[V_X],      vec[V_Y],      vec[V_Z]);
			pclos();
			sjm = vdevice.save % 2;     /* save new point */
			tmesh[sjm][V_X] = vec[V_X];
			tmesh[sjm][V_Y] = vec[V_Y];
			tmesh[sjm][V_Z] = vec[V_Z];
		}
		vdevice.save++;
		break;
	case VQSTRIP:
		if( vdevice.save < 3 ) {  /* save the first three calls */
			qstrip[vdevice.save][V_X] = vec[V_X];
			qstrip[vdevice.save][V_Y] = vec[V_Y];
			qstrip[vdevice.save][V_Z] = vec[V_Z];
		} else { 
			sjm = vdevice.save % 2;
			if (sjm) { /* sjm !=0; plot quadralateral */
				pmv(qstrip[0][V_X], qstrip[0][V_Y], qstrip[0][V_Z]);
				pdr(qstrip[1][V_X], qstrip[1][V_Y], qstrip[1][V_Z]);
				pdr(      vec[V_X],       vec[V_Y],       vec[V_Z]);
				pdr(qstrip[2][V_X], qstrip[2][V_Y], qstrip[2][V_Z]);
				pclos();
				/* save 2 newest points over 2 oldest */
				qstrip[0][V_X] = qstrip[2][V_X]; 
				qstrip[0][V_Y] = qstrip[2][V_Y];
				qstrip[0][V_Z] = qstrip[2][V_Z];
				qstrip[1][V_X] = vec[V_X];
				qstrip[1][V_Y] = vec[V_Y];
				qstrip[1][V_Z] = vec[V_Z];
			} else {  /* sjm == 0; save vec in register 2 */
				qstrip[2][V_X] = vec[V_X];
				qstrip[2][V_Y] = vec[V_Y];
				qstrip[2][V_Z] = vec[V_Z];
			}
		}
		vdevice.save++;
		break;
	default:
		move(vec[V_X], vec[V_Y], vec[V_Z]);
	}
}
		
/*
 * v4f	
 * 	Adds a 4D point to our fake point buffer
 */
void
v4f(vec)
	float	vec[4];
{
	vcall(vec, 4);
}

/*
 * v3f	
 * 	Adds a 3D point to our fake point buffer
 */
void
v3f(vec)
	float	vec[3];
{
	vcall(vec, 3);
}

/*
 * v2f	
 * 	Adds a 2D point to our fake point buffer
 */
void
v2f(vec)
	float	vec[2];
{
	vcall(vec, 2);
}


/*
 * v4d	
 * 	Adds a 4D point to our fake point buffer
 */
void
v4d(vec)
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
 * v3d	
 * 	Adds a 3D point to our fake point buffer
 */
void
v3d(vec)
	double	vec[3];
{
	float	v[3];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];

	vcall(v, 3);
}

/*
 * v2d	
 * 	Adds a 2D point to our fake point buffer
 */
void
v2d(vec)
	double	vec[2];
{
	float	v[2];

	v[0] = vec[0];
	v[1] = vec[1];

	vcall(v, 2);
}



/*
 * v4i	
 * 	Adds a 4D point to our fake point buffer
 */
void
v4i(vec)
	long	vec[4];
{
	float	v[4];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];
	v[3] = vec[3];

	vcall(v, 4);
}

/*
 * v3i	
 * 	Adds a 3D point to our fake point buffer
 */
void
v3i(vec)
	long	vec[3];
{
	float	v[3];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];

	vcall(v, 3);
}

/*
 * v2i	
 * 	Adds a 2D point to our fake point buffer
 */
void
v2i(vec)
	long	vec[2];
{
	float	v[2];

	v[0] = vec[0];
	v[1] = vec[1];

	vcall(v, 2);
}
/*
 * v4s	
 * 	Adds a 4D point to our fake point buffer
 */
void
v4s(vec)
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
 * v3s	
 * 	Adds a 3D point to our fake point buffer
 */
void
v3s(vec)
	short	vec[3];
{
	float	v[3];

	v[0] = vec[0];
	v[1] = vec[1];
	v[2] = vec[2];

	vcall(v, 3);
}

/*
 * v2s	
 * 	Adds a 2D point to our fake point buffer
 */
void
v2s(vec)
	short	vec[2];
{
	float	v[2];

	v[0] = vec[0];
	v[1] = vec[1];

	vcall(v, 2);
}
