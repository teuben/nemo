/* $Header */
#ifndef STANDARD_MACRO
#define STANDARD_MACRO
#define	MAX(x, y) ((x) > (y) ? (x) : (y))
#define	MIN(x, y) ((x) < (y) ? (x) : (y))
#define	ABS(x) ((x) < 0 ? -(x) : (x))
#define ROUND(x)   ((int)(((x)<0)?((x)-0.5):((x)+0.5)))
#define CEILING(x) (((x)>=0)?(int)(x)==(x)?(int)(x):(int)((x)+1):(int)(x))
#define FLOOR(x)   (((x)>=0)?(int)(x):(int)(x)==(x)?(int)(x):(int)((x)-1))
#define EQ(s, t)	(!strcmp(s, t))
#define EQN(s, t, n)	(!strncmp(s, t, n))
#endif STANDARD_MACRO
