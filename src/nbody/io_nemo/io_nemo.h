/* --------------------------------------------------------- *\
|* io_nemo.h :
|*
\* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

char * allocate_pointer(char *, int );

void reajust_ptr();

int io_nemo(char * , char * , ...);

int close_io_nemo(char * );

#ifdef __cplusplus
}
#endif
