//
// C++ declarations of the C NEMO-STARLAB interface functions 
//
// note this currently only works when real==double

#ifdef __cplusplus
extern "C" {
#endif

void check_real(int size);

void put_snap_c(string fname, 
		 int nbody,
                 double *mass, 
                 double *pos, double *vel);

int get_snap_c(string fname, 
                 double **mass, 
                 double **pos, double **vel);

#ifdef __cplusplus
}
#endif


