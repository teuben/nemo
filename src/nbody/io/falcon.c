/*
 *   Process falcon's snapshot files
 *
 *   Weird:   nbody limited to 32bit (4e9)
 *            mass is double, rest float
 *

GROUP "/" {
   ATTRIBUTE "falcON" {
   ATTRIBUTE "num_snapshots" {
   DATASET "history" {
   GROUP "snapshot0" {           group1
      ATTRIBUTE "Nstd" {
      ATTRIBUTE "time" {
      GROUP "std" {              group2
         ATTRIBUTE "N" {
         DATASET "acc" {
         DATASET "mass" {
         DATASET "pos" {
         DATASET "pot" {
         DATASET "vel" {
   GROUP "snapshot1" {

 */


#include <nemo.h>
#include <history.h>

#include <snapshot/snapshot.h>

#include <hdf5.h>

string defv[] = {
  "in=???\n            Input snapshot",
  "out=???\n           Output snapshot",
  "VERSION=0.2\n       21-apr-2024 PJT",
  NULL,
};

string usage="Process falcon HDF5 snapshots";

#define MAXBODY  100000

void rdwt(string fin, string fout);

void nemo_main()
{
  string fin = getparam("in");
  string fout = getparam("out");

  rdwt(fin, fout);

}



void rdwt(string fin, string fout)
{
  
  hid_t       file_id, attr_id, group_id, group2_id, dataset_id;
  herr_t      status;
  hsize_t     dims[2];
  stream      ostr;
  int         i, j, cs = CSCode(Cartesian, NDIM, 2);
  int         ntime, nbody, nbody2;
  double      stime;
  double      mass[MAXBODY];
  float       acc[MAXBODY][3], pot[MAXBODY], pos[MAXBODY][3], vel[MAXBODY][3];
  

  file_id = H5Fopen(fin, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id == H5I_INVALID_HID) error("not valid file");

  status = H5Aexists(file_id, "falcON");
  if (status == 0) error("not a falcON file");
  // @todo:   check if value is 67329

  ostr = stropen(fout, "w");
  put_history(ostr);
    
  attr_id = H5Aopen(file_id, "num_snapshots", H5P_DEFAULT);
  if (attr_id == H5I_INVALID_HID) error("not valid attribute1");  
  dprintf(2,"attr_id: %d\n", attr_id);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &ntime);
  dprintf(1,"ntime = %d\n",ntime);
  H5Aclose(attr_id);

  // @todo  grab the history (DATA)

  for (int itime=0; itime<ntime; itime++) {
    /* open group */
    char snapshot_name[20];
    sprintf(snapshot_name,"snapshot%d", itime);
    group_id = H5Gopen1(file_id, snapshot_name);
    if (group_id == H5I_INVALID_HID) error("not valid group %s",snapshot_name);
    dprintf(2,"group_id: %d\n", group_id);

    // time
    attr_id = H5Aopen(group_id, "time", H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) error("not valid attribute2");      
    dprintf(2,"attr_id: %d\n", attr_id);
    status = H5Aread(attr_id, H5T_IEEE_F64LE, &stime);
    H5Aclose(attr_id);
    
    // Nstd
    attr_id = H5Aopen(group_id, "Nstd", H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) error("not valid attribute2");      
    dprintf(2,"attr_id: %d\n", attr_id);
    status = H5Aread(attr_id, H5T_NATIVE_INT, &nbody);
    H5Aclose(attr_id);
    
    dprintf(0,"time = %g   nbody = %d\n",stime, nbody);

    group2_id = H5Gopen1(group_id, "std"); 
    if (group2_id == H5I_INVALID_HID) error("not valid group2");

    attr_id = H5Aopen(group2_id, "N", H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) error("not valid attribute3");      
    dprintf(2,"attr_id: %d\n", attr_id);
    status = H5Aread(attr_id, H5T_NATIVE_INT, &nbody2);
    dprintf(1,"nbody_std = %d\n",nbody2);
    H5Aclose(attr_id);
    if (nbody2 != nbody)
      warning("Cannot probably process non-std snapshots; nbody_std=%d",nbody2);

    dataset_id = H5Dopen2(group2_id, "mass", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("not valid mass dataset");
    status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
    H5Dclose(dataset_id);
    dprintf(1,"mass[0] = %g\n", mass[0]);

    dataset_id = H5Dopen2(group2_id, "pos", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("not valid pos dataset");
    status = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
    H5Dclose(dataset_id);
    dprintf(1,"pos[0] = %g %g %g\n", pos[0][0], pos[0][1], pos[0][2]);

    dataset_id = H5Dopen2(group2_id, "vel", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("not valid vel dataset");
    status = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel);
    H5Dclose(dataset_id);
    dprintf(1,"vel[0] = %g %g %g\n", vel[0][0], vel[0][1], vel[0][2]);

    dataset_id = H5Dopen2(group2_id, "pot", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("not valid pot dataset");
    status = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot);
    H5Dclose(dataset_id);
    dprintf(1,"pot[0] = %g\n", pot[0]);

    dataset_id = H5Dopen2(group2_id, "acc", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("not valid acc dataset");
    status = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, acc);
    H5Dclose(dataset_id);
    dprintf(1,"acc[0] = %g %g %g\n", acc[0][0], acc[0][1], acc[0][2]);
    
    status = H5Gclose(group2_id);
    status = H5Gclose(group_id);

    put_set(ostr, SnapShotTag);
    put_set(ostr, ParametersTag);
    put_data(ostr, NobjTag, IntType, &nbody, 0);
    put_data(ostr, TimeTag, DoubleType, &stime, 0);
    put_tes(ostr, ParametersTag);
    put_set(ostr, ParticlesTag);
    put_data(ostr, CoordSystemTag, IntType, &cs, 0);
    put_data(ostr, MassTag, DoubleType, mass, nbody, 0);
    put_data(ostr, PosTag, FloatType, pos[0], nbody, NDIM, 0);
    put_data(ostr, VelTag, FloatType, vel[0], nbody, NDIM, 0); 
    put_data(ostr, PotentialTag, FloatType, pot, nbody, 0);   
    put_data(ostr, AccelerationTag, FloatType, acc[0], nbody, NDIM, 0); 
    put_tes(ostr, ParticlesTag);
    put_tes(ostr, SnapShotTag);
  } // for(itime)
  status = H5Fclose(file_id);
  strclose(ostr);
}

