//
// C++ declarations of the C NEMO-STARLAB interface functions 
//
// note this currently only works when real==double

#ifdef __cplusplus
extern "C" {
#endif

void check_real(int size);

void put_snap_c(string fname, string hline,
		 int nbody,
		 double time,
                 double *mass, 
                 double *pos, double *vel, double *acc, double *aux, double *phi, int *key);

int get_snap_c(string fname, 
		 double *time,
                 double **mass, 
                 double **pos, double **vel);

#ifdef __cplusplus
}
#endif


