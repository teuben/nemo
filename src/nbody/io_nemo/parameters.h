/* --------------------------------------------------------- *\
|* io_init.h :
|*
\* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

int get_case(char * );

char * get_champs(char **);

int chk_parameters(bool , int );

float char2float(char * ,int);

double char2double(char * ,int );

#ifdef __cplusplus
}
#endif
