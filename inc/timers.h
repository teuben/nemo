
/* timers.c */

/* @todo:  should use the new more portable cycle.h  */

long long readTSC(void);	/* not sure if we keep this public */
void init_timers(int n);
void init_timers2(int n, int mode);
void stamp_timers(int i);
long long diff_timers(int i, int j);
double diff_timers2(int i, int j);
