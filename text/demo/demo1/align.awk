#! /usr/bin/awk -f
# 
#   align phases in columns 2 and 3 (time is col 1)
#   which are assumed in degrees. Output from snapinert


#   assumed phase declines !!!

BEGIN {
    xoff = 0;
    yoff = 0;
}

{
    if (NR == 1) {
	# print $1,$2,$3,0,0;
	oldt = $1;
	oldx = $2;
	oldy = $3;
    } else {
	newdx = $2-oldx;
	newdy = $3-oldy;
	if (newdx > 0) {
	    xoff -= 180;
	    yoff -= 180;
	    newdx -=180;
	    newdy -=180;
	} 
	dt = $1-oldt;
	print $1,xoff+$2,yoff+$3,newdx/dt,newdy/dt
	oldt = $1;
	oldx = $2;
	oldy = $3;
    }
}
