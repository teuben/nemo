/* -------------------------------------------------------------- *\
|* file_data.h :	JC LAMBERT	 22-Oct-97	 V1.0
|*
|* V1.0  : created
|*
|* Header file 
|*
\* -------------------------------------------------------------- */

/* Variables used for reading operations */
extern char * io_in[]; 
extern FILE * instr[];
extern bool   read_one[];

/* Variables used for writing operations */
extern char * io_out[]; 
extern FILE * outstr[];
extern bool   save_one[];

/* variables used to keep track of history */
extern bool set_history[];
extern char * history_prog;
extern char * hist_file;        /* history file name */
/* -------------------------------------------------------------- *\
|* End of file_dat.h
\* -------------------------------------------------------------- */
