/* --------------------------------------------------------- *\
|* toolsnemo.h :
|*
\* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

int get_data_gen(stream, char *, char *, int, int, int, int, char **);

int get_data_time(stream, char *, int, char **);

int get_data_nbody(stream, char *, int, int **);

int get_data_mass(stream, char *, int, int, char **);

int get_data_phase(stream, char *, int, int, char **, int);

int get_data_pot(stream, char *, int, int, char **);

int get_data_acc(stream, char *, int, int, char **, int);

bool ** chk_select(int * , int, int, string select_pts[]);

#ifdef __cplusplus
}
#endif
