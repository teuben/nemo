#include <stdio.h>
#include "vogl.h"

/*
 * makeobj_
 */
void
makeobj_(n)
	int	*n;
{
	makeobj((Object)*n);
}

/*
 * makeob_
 */
void
makeob_(n)
	int	*n;
{
	makeobj((Object)*n);
}

/*
 * closeobj_
 */
void
closeobj_()
{
	closeobj();
}

/*
 * closeo_
 */
void
closeo_()
{
	closeobj();
}

/*
 * delobj_
 */
void
delobj_(n)
	int	*n;
{
	delobj((Object)*n);
}

/*
 * genobj_
 */
int
genobj_()
{
	return(genobj());
}

/*
 * getopenobj_
 */
int
getopenobj_()
{
	return(getopenobj());
}

/*
 * getope_
 */
int
getope_()
{
	return(getopenobj());
}

/*
 * callobj_
 */
void
callobj_(n)
	int	*n;
{
	callobj((Object)*n);
}

/*
 * callob_
 */
void
callob_(n)
	int	*n;
{
	callobj((Object)*n);
}
/*
 * isobj_
 */
int
isobj_(n)
	int	*n;
{
	return(isobj((Object)*n));
}
