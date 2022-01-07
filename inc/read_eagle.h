#include "hdf5.h"

#ifndef _READ_EAGLE_H
#define _READ_EAGLE_H

#define MAX_NAMELEN 500

typedef unsigned long long peanokey;

peanokey peano_hilbert_key(int x, int y, int z, int bits);

/*
  Type to store information about an open snapshot
*/
typedef struct
{
  double         boxsize;
  int            numfiles;
  int            hashbits;
  char           basename[MAX_NAMELEN];
  unsigned char *hashmap;
  int            ncell;
  int            nhash;
  long long     *first_key_in_file[6];
  long long     *last_key_in_file[6];
  long long     *num_keys_in_file[6];
  unsigned int **part_per_cell[6];
  unsigned int **first_in_cell[6];
  long long      numpart_total[6];
  int            num_datasets[6];
  char          *dataset_name[6];
  double         sampling_rate;
  int            split_rank, split_size;
} EagleSnapshot;

/* 
   Type codes for returning the type of a dataset
   (mostly for the use of Fortran/Python/IDL wrappers)
*/
typedef enum e_TypeCode
  {
    t_int       = 0,
    t_long_long = 1,
    t_float     = 2,
    t_double    = 3,
    t_uint      = 4,
    t_ulong_long = 5
  } TypeCode;

/* Return a pointer to the last error message */
char *get_error(void);

/* Set whether we should abort on errors */
void abort_on_error(int flag);

/*
  Function to open a snapshot

  Parameters:
    fname - name of one file in the snapshot

  Return value:
    Success - returns a pointer to a new EagleSnapshot
    Failure - returns a null pointer

  The EagleSnapshot should be deallocated with a call to
  close_snapshot() to avoid memory leaks.
*/
EagleSnapshot *open_snapshot(char *fname);


/*
  Function to close a snapshot

  Parameters:
    snap - pointer to the EagleSnapshot
*/
void close_snapshot(EagleSnapshot *snap);


/*
  Function to select a region in a snapshot

  Parameters:
    snap  - pointer to the EagleSnapshot
    xmin  - minimum x coordinate
    xmax  - maximum x coordinate
    ymin  - minimum y coordinate
    ymax  - maximum y coordinate
    zmin  - minimum z coordinate
    zmax  - maximum z coordinate

 Repeated calls to select_region() can be used to select
 oddly shaped or disjoint regions.
*/
void select_region(EagleSnapshot *snap, 
		  double xmin, double xmax,
		  double ymin, double ymax,
		  double zmin, double zmax);

/*
  Alternate version of select region where integer cell
  coordinates are specified
*/
void select_grid_cells(EagleSnapshot *snap, 
		       int ixmin, int ixmax,
		       int iymin, int iymax,
		       int izmin, int izmax);

/*
  This version allows regions which are not axis aligned
*/
void select_rotated_region(EagleSnapshot *snap,
			   double *centre,
			   double *xvec, double *yvec, double *zvec,
			   double *length);

/*
  Function to set sampling rate
*/
void set_sampling_rate(EagleSnapshot *snap, double rate);

/*
  Clear any selection associated with the specified snapshot

  Parameters:
    snap - pointer to the EagleSnapshot
*/
void clear_selection(EagleSnapshot *snap);


/*
  Count the particles in the selected region

  Parameters:
    snap      - pointer to the eagle snapshot
    itype     - particle type to read

  Return value:
    Success - Number of particles in the selected region
    Failure - a negative number
*/
long long count_particles_with_index(EagleSnapshot *snap, int itype, 
				     int *file_index, int *file_offset,
				     size_t nmax);

/*
  Macros to call count_particles_with_index()
*/
#define count_particles(snap, itype) count_particles_with_index((snap), (itype), NULL, NULL, 0)
#define get_particle_locations count_particles_with_index

/*
  Function to read a dataset for all particles in the
  selected region. Works with 1D datasets and 2D Nx3 datasets.

  Parameters:
    snap      - pointer to the eagle snapshot
    itype     - particle type to read
    name      - HDF5 dataset name relative to the PartTypeX group
    hdf5_type - HDF5 type of the output buffer 'buf'
    buf       - buffer in which to store the result
    n         - size of the buffer 'buf'
    extra_basename - if this is set, read from alternate files <extra_basename>.*.hdf5

  Return value:
    Success - the number of particles read
    Failure - a negative number
*/
long long read_extra_dataset(EagleSnapshot *snap, int itype, char *name, hid_t hdf5_type, void *buf, size_t n,
			     char *extra_basename);
#define read_dataset(snap, itype, name, hdf5_type, buf, n) read_extra_dataset((snap), (itype), (name), (hdf5_type), (buf), (n), NULL)

/*
  Macros for reading data into buffers of specific types.

  These avoid the need to use quantities from hdf5.h in the
  calling program. Parameters have the same meaning as in
  read_dataset().
*/
#define read_dataset_int(snap,itype,name,buf,n) read_extra_dataset(snap,itype,name,H5T_NATIVE_INT,buf,n,NULL)
#define read_dataset_float(snap,itype,name,buf,n) read_extra_dataset(snap,itype,name,H5T_NATIVE_FLOAT,buf,n,NULL)
#define read_dataset_double(snap,itype,name,buf,n) read_extra_dataset(snap,itype,name,H5T_NATIVE_DOUBLE,buf,n,NULL)
#define read_dataset_long_long(snap,itype,name,buf,n) read_extra_dataset(snap,itype,name,H5T_NATIVE_LLONG,buf,n,NULL)
#define read_dataset_unsigned_int(snap,itype,name,buf,n) read_extra_dataset(snap,itype,name,H5T_NATIVE_UINT,buf,n,NULL)
#define read_dataset_unsigned_long_long(snap,itype,name,buf,n) read_extra_dataset(snap,itype,name,H5T_NATIVE_ULLONG,buf,n,NULL)

/*
  Return information about a dataset given its name

  Parameters:
    snap      - pointer to the eagle snapshot
    itype     - particle type to read
    name      - HDF5 dataset name relative to the PartTypeX group
    typecode  - 
    rank      - the rank of the dataset (1 for scalar particle
                properties, 2 for vectors)

 Return value:
    Success - zero
    Failure - non-zero

*/
int get_extra_dataset_info(EagleSnapshot *snap, int itype, char *dset_name, TypeCode *typecode, int *rank, char *extra_basename);
#define get_dataset_info(snap,itype,dset_name,typecode,rank) get_extra_dataset_info(snap, itype, dset_name, typecode, rank, NULL)

/*
  Return the number of datasets available for the specified particle type

  Parameters:
    snap      - pointer to the eagle snapshot
    itype     - particle type to read

  Return value:
    Success - the number of datasets
    Failure - a negative number

*/
int get_dataset_count(EagleSnapshot *snap, int itype);

/*
  Get the name of the specified dataset.

  Parameters:
    snap      - pointer to the eagle snapshot
    itype     - particle type to read
    iset      - index of the dataset
    buf       - pointer to the output buffer
    len       - length of the output buffer

  Return value:
    Success - the length of the string copied to the output buffer
    Failure - a negative number

*/
int get_dataset_name(EagleSnapshot *snap, int itype, int iset, char *buf, size_t len);

/*
  Split the region to be read between processors.
  For use in MPI programs. When using MPI all read_eagle 
  calls are collective.
  
  Parameters:
    ThisTask - MPI rank of this process
    NTask    - Number of MPI ranks
  
  Return value:
    Success - returns zero
    Failure - a negative number

  This must be called after select_region and before read_dataset.
  Each selection can only be split once. I.e. must call clear_selection
  before split_selection can be called again.

*/
int split_selection(EagleSnapshot *snap, int ThisTask, int NTask);


double random_double(void);
void set_random_seed(void);

#endif
