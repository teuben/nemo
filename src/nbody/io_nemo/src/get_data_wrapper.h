/* -------------------------------------------------------------- *\
|* $Id$
|*
|* Wrapper of basic NEMO procedure, it makes me life easiest :)
|*
\* -------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

int get_data_gen(stream, char *, char *, int, int, int, int, char **);

int get_data_time(stream, char *, int, char **);

int get_data_nbody(stream, char *, int, int **);

int get_data_mass(stream, char *, int, int, char **);

int get_data_keys(stream, char *, int, int, char **);

int get_data_phase(stream, char *, int, int, char **, int);

int get_data_pos(stream, char *, int, int, char **, int);

int get_data_vel(stream, char *, int, int, char **, int);

int get_data_pot(stream, char *, int, int, char **);

int get_data_acc(stream, char *, int, int, char **, int);

#ifdef __cplusplus
}
#endif
