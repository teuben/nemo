/* -------------------------------------------------------------- *\
|* snapshot_data.h :	JC LAMBERT	 22-Oct-97	 V1.0
|*
|* V1.0  : created
|*
|* Header file 
|*
\* -------------------------------------------------------------- */
extern char  * pos, ** pos_p,    /* position */
             * vel, ** vel_p,    /* velocity */
             * pot, ** pot_p,    /* potential */
             * acc, ** acc_p,    /* acceleration */
             * mass,** mass_p,   /* mass */
             * timu,** time_p,   /* time steps */
             * selt,** selt_p,   /* selected time */
			 * SelectionString123; /* particles selection */

extern int   * nbody,** nbody_p;   /* nbody */
/* -------------------------------------------------------------- *\
|* End of snapshot_data.h
\* -------------------------------------------------------------- */
