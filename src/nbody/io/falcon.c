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
#include <strlib.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include <hdf5.h>

string defv[] = {
  "in=???\n            Input snapshot",
  "out=???\n           Output snapshot",
  "type=std\n          Particle Group",
  // times=
  "VERSION=0.5\n       21-apr-2024 PJT",
  NULL,
};

string usage="Convert falcon HDF5 snapshot to NEMO snapshot";

// @todo  undef this for trying dynamic arrays
#define MAXBODY  10000000

#ifdef MAXBODY
  float       mass[MAXBODY], pot[MAXBODY];
  float       acc[MAXBODY][3], pos[MAXBODY][3], vel[MAXBODY][3];
#else
  // @todo use mdarray, but figure out if real is float or double
  mdarray1    *mass, *pot;
  mdarray2    *acc, *pos, *vel;
#endif

string file_pipe();


void nemo_main()
{
  string      fin = getparam("in");
  string      fout = getparam("out");
  string      type = getparam("type");
  hid_t       file_id, attr_id, group_id, group2_id, dataset_id;
  herr_t      status;
  stream      ostr;
  int         cs = CSCode(Cartesian, NDIM, 2);
  int         ntime, nbody, nbody1, nbody2;
  double      stime;
  bool        Qpipe = streq(fin,"-");
  char        ntype[32];

  if (Qpipe)
    fin = file_pipe();

  file_id = H5Fopen(fin, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id == H5I_INVALID_HID) error("H5Fopen: not valid file");

  status = H5Aexists(file_id, "falcON");
  if (status == 0) error("H5Aexists: not a falcON file");
  // @todo:   check if value is 67329

  ostr = stropen(fout, "w");
  put_history(ostr);
    
  attr_id = H5Aopen(file_id, "num_snapshots", H5P_DEFAULT);
  if (attr_id == H5I_INVALID_HID) error("num_snapshots: not valid attribute1");  
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
    if (group_id == H5I_INVALID_HID) error("%s: not valid group",snapshot_name);
    dprintf(2,"group_id: %d\n", group_id);

    // time
    attr_id = H5Aopen(group_id, "time", H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) error("time: not valid attribute2");      
    dprintf(2,"attr_id: %d\n", attr_id);
    status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &stime);
    H5Aclose(attr_id);
    
    // Nstd (or type)
    sprintf(ntype,"N%s", type);
    attr_id = H5Aopen(group_id, ntype, H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) error("%s: not valid attribute2",ntype);      
    dprintf(2,"attr_id: %d\n", attr_id);
    status = H5Aread(attr_id, H5T_NATIVE_INT, &nbody);
    H5Aclose(attr_id);
    
    dprintf(0,"time = %g   nbody = %d  type = %s\n",stime, nbody, type);

    group2_id = H5Gopen1(group_id, type); 
    if (group2_id == H5I_INVALID_HID) error("%s: not valid group2", type);

    attr_id = H5Aopen(group2_id, "N", H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) error("N: not valid attribute3");      
    dprintf(2,"attr_id: %d\n", attr_id);
    status = H5Aread(attr_id, H5T_NATIVE_INT, &nbody2);
    dprintf(1,"nbody_std = %d\n",nbody2);
    H5Aclose(attr_id);
    if (nbody2 != nbody)
      warning("Cannot probably process snapshots; nbody_%s=%d",type,nbody2);
    
#ifdef MAXBODY
    if (nbody2 > MAXBODY)
      error("Not enough space for %d particles; MAXBODY=%d", nbody2, MAXBODY);
#else
    if (itime == 0) {
      nbody1 = nbody2;
      mass = (double *) allocate(nbody2*sizeof(double));
      pos  = (float *)  allocate(nbody2*sizeof(double)*3);
      vel  = (float *)  allocate(nbody2*sizeof(double)*3);
      pot  = (float *)  allocate(nbody2*sizeof(double));
      acc  = (float *)  allocate(nbody2*sizeof(double)*3);
    } else {
      if (nbody2 > nbody1) error("Cannot handle increasing particle number yet");
    }
#endif    

    dataset_id = H5Dopen2(group2_id, "mass", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("mass: not valid dataset");
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
    H5Dclose(dataset_id);
    dprintf(1,"mass[0] = %g\n", mass[0]);

    dataset_id = H5Dopen2(group2_id, "pos", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("pos: not valid dataset");
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
    H5Dclose(dataset_id);
    dprintf(1,"pos[0] = %g %g %g\n", pos[0][0], pos[0][1], pos[0][2]);

    dataset_id = H5Dopen2(group2_id, "vel", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("vel: not valid dataset");
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel);
    H5Dclose(dataset_id);
    dprintf(1,"vel[0] = %g %g %g\n", vel[0][0], vel[0][1], vel[0][2]);

    dataset_id = H5Dopen2(group2_id, "pot", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("pot: not valid dataset");
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pot);
    H5Dclose(dataset_id);
    dprintf(1,"pot[0] = %g\n", pot[0]);

    dataset_id = H5Dopen2(group2_id, "acc", H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) error("acc: not valid dataset");
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, acc);
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
    put_data(ostr, MassTag, FloatType, mass, nbody, 0);
    put_data(ostr, PosTag, FloatType, pos[0], nbody, NDIM, 0);
    put_data(ostr, VelTag, FloatType, vel[0], nbody, NDIM, 0); 
    put_data(ostr, PotentialTag, FloatType, pot, nbody, 0);   
    put_data(ostr, AccelerationTag, FloatType, acc[0], nbody, NDIM, 0); 
    put_tes(ostr, ParticlesTag);
    put_tes(ostr, SnapShotTag);
  } // for(itime)
  status = H5Fclose(file_id);
  strclose(ostr);

   if (Qpipe) remove(fin);
}

string file_pipe()
{
  stream istr, ostr;
  char tempname[MAXPATHLEN];
  char buffer[1024*1024];          /* 1 MB buffer */
  int k, n, ndat = 1024*1024;      /* that buffer size */
  int fds;

  strcpy(tempname, "/tmp/scrNemo.XXXXXX");
#if 1
  fds = mkstemp(tempname);
#endif
  ostr = stropen(tempname,"w!");

  fds = dup(fileno(stdin));
  istr = (stream ) fdopen(fds, "r");

  k = 0;
  while (1) {
    n = fread(buffer, 1, ndat, istr);
    if (n<0) error("file_pipe: error reading data from pipe");
    fwrite(buffer, 1, n, ostr);
    k+=n;
    if (n<ndat) break;
  }
  dprintf(0,"file_pipe: read %d bytes from HDF pipe, written to %s\n",k,tempname);
  strclose(ostr);
  strclose(istr);
    
  return scopy(tempname);
}
