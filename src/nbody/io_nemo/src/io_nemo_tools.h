/* --------------------------------------------------------- *\
|* $Id$
|*
|* Functions prototypes
|*
\* --------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

float   char2float(char * ,int);

double  char2double(char * ,int );

char *  allocate_pointer(char *, int );

char *  f_ch_to_c(char * ,int );

bool ** chk_select(int * ,int ,int ,string select_pts[]);

char *  get_selected(char *);

#ifdef __cplusplus
}
#endif


