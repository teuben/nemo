/* --------------------------------------------------------- *\
|* $Id$
|*
|* Functions prototypes
|*
\* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

int get_history_input_file(char *);

int get_old_file(char * curr_file ,
		 char * io_file[] ,
		 bool   io_one[]  ,
		 FILE * io_str[]  ,
		 int);

int get_new_file(char * curr_file ,
		 char * io_file[] ,
		 bool   io_one[]  ,
		 FILE * io_str[]  ,
		 char * mode,
		 int);

#ifdef __cplusplus
}
#endif
