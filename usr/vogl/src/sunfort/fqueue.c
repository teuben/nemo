#include "vogl.h"
#include "vodevice.h"
#include <stdio.h>

/*
 * qdevice_
 */
void
qdevice_(dev)
	int	*dev;
{
	qdevice((Device)*dev);
}

/*
 * qdevic_	(same as qdevice_)
 */
void
qdevic_(dev)
	int	*dev;
{
	qdevice((Device)*dev);
}

/*
 * unqdevice_
 */
void
unqdevice_(dev)
	int	*dev;
{
	unqdevice((Device)*dev);
}

/*
 * unqdev_	(same as unqdevice_)
 */
void
unqdev_(dev)
	int	*dev;
{
	unqdevice((Device)*dev);
}

/*
 * qread_
 */
int
qread_(val)
	short	*val;
{
	return(qread(val));
}

/*
 * qtest_
 */
int
qtest_()
{
	int	ret = qtest();

	return(ret);
}

/*
 * qreset_
 */
void
qreset_()
{
	qreset();
}

/*
 * isqueued_
 */
int
isqueued_(dev)
	int	*dev;
{
	return((int)isqueued((Device)*dev));
}

/*
 * isqueu_	(same as isqueued_)
 */
int
isqueu_(dev)
	int	*dev;
{
	return((int)isqueued((Device)*dev));
}

/*
 * qenter_
 */
void
qenter_(dev, val)
	int	*dev, *val;
{
	qenter((Device)*dev, (short)*val);
}
