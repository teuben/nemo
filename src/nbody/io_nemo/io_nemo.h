/* -------------------------------------------------------------- *\
|* io_nemo.h :	JC LAMBERT	 22-Oct-97	 V1.0
|*
|* V1.0  : created
|*
|* Header file for IO_NEMO library
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Extern Shared variables declared in "io_nemo_dmp.c" file
\* -------------------------------------------------------------- */

/* variables for reading */
 extern char * io_in[MAXIO]; 
 extern FILE * instr[MAXIO];
 extern bool   read_one[MAXIO];

/* variables for writing */
 extern char * io_out[MAXIO]; 
 extern FILE * outstr[MAXIO];
 extern bool   save_one[MAXIO];

/* snapshot variables */
 extern char * pos,
             * vel,
             * pot,
             * acc,
             * mass,
             * timu,
             * selt,
             * selp;

 extern int  * nbody;


/* flags */
 extern int N_io, T_io, M_io, X_io, V_io, P_io, A_io,
           F_dim, S_io, R_io, C_io, I_io, ST_io, SP_io ; /*MODIFIED*/

/* -------------------------------------------------------------- *\
|* End of io_nemo.h
\* -------------------------------------------------------------- */
