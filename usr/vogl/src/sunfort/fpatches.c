#include "vogl.h"

/*
 * patchbasis_
 */
void
patchbasis_(ubasisid, vbasisid)
        int	*ubasisid, *vbasisid;
{
	patchbasis((long)ubasisid, (long)vbasisid);
}

/*
 * patchb_	(same as patchbasis_)
 */
void
patchb_(ubasisid, vbasisid)
        int	*ubasisid, *vbasisid;
{
	patchbasis((long)*ubasisid, (long)*vbasisid);
}

/*
 * patchprecision_
 */
void
patchprecision_(useg, vseg)
        int     *useg, *vseg;
{
	patchprecision((long)*useg, (long)*vseg);
}

/*
 * patchp_	(same as patchprecision_)
 */
void
patchp_(useg, vseg)
        int     *useg, *vseg;
{
	patchprecision((long)*useg, (long)*vseg);
}

/* 
 * patchcurves_
 */
void
patchcurves_(nu, nv)
	int	*nu, *nv;
{
	patchcurves((long)*nu, (long)*nv);
}

/* 
 * patchc_	(same as patchcurves_)
 */
void
patchc_(nu, nv)
	int	*nu, *nv;
{
	patchcurves((long)*nu, (long)*nv);
}

/*
 * patch_
 */
void
patch_(geomx, geomy, geomz)
	Matrix	geomx, geomy, geomz;
{
	patch(geomx, geomy, geomz);
}

/*
 * rpatch_
 */
void
rpatch_(geomx, geomy, geomz, geomw)
	Matrix	geomx, geomy, geomz, geomw;
{
	rpatch(geomx, geomy, geomz, geomw);
}

/*
 * defbasis_
 */
void
defbasis_(id, array)
	int	*id;
	float	array[4][4];
{
	defbasis((short)*id, array);
}

/*
 * defbas_	(same as defbasis)
 */
void
defbas_(id, array)
	int	*id;
	float	array[4][4];
{
	defbasis((short)*id, array);
}
