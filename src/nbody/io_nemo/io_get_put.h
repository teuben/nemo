/* --------------------------------------------------------- *\
|* io_get_put.h :
|*
\* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

int put_data_select(char * outfile,
		    int   rtype);

int get_data_select(char * infile,
		    int    rtype);


#ifdef __cplusplus
}
#endif
