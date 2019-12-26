/*
 *   Transform an orbit and ODM into "BOOM" HDF
 */

#include <nemo.h>
#include <hdf5.h>

string defv[] = {
  "in=\n               Input orbit",
  "odm=\n              Input ODM (an image)",
  "out=junk.hdf\n      Output boom HDF",
  "ndim=3\n            2 or 3 dimensional",            
  "VERSION=0.2\n       25-Dec-2019 XYZ",
  NULL,
};

string usage="Transform an orbit and ODM into 'BOOM' HDF [testing]";

string cvsid="$Id:$";

void nemo_main()
{
  string fin = getparam("in");
  string fodm = getparam("odm");
  string fout = getparam("out");
  float version = 0.3;

  warning("This program only understands the out= keyword");

  hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
  hsize_t     dims[2];
  herr_t      status;

   /* Create a new file using default properties. */
  dprintf(0,"Data %s\n", fout);
  file_id = H5Fcreate(fout, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


   /* Create the data space for the dataset. */
  dims[0] = 4; 
  dims[1] = 6; 
  dataspace_id = H5Screate_simple(2, dims, NULL);

   /* Create the dataset. */
  dataset_id = H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);

  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);

  /* Close the file. */
  status = H5Fclose(file_id);
}

