

struct	YOBBA	{
		union {
			struct {
				unsigned alpha : 1;
				unsigned beta  : 2;
				unsigned fixe  : 3;
				unsigned kralb  : 1;
			} vals;
			struct {
				unsigned char yobbav;
			} yobbavals;
		} yobba;
} *yobbaray;

/*
 * yobbarays
 *
 *	Turns on (or off) yobba rays, as described by Larry Dart's friend.
 *
 *	onoff <> 0 - YOBBARAYS ON.
 *	onoff =  0 - YOBBARAYS OFF.
 */
void
yobbarays(onoff)
	int	onoff;
{
	yobbaray = (struct YOBBA *)onoff;
}
