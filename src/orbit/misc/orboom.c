/*
 *   Transform an orbit and ODM into "BOOM" HDF
 */

#include <nemo.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <orbit.h>
#include <image.h>
#include <hdf5.h>

string defv[] = {
  "in=???\n            Input snapshot series",
  "out=???\n           Output BOOM",
  "ndim=3\n            2 or 3 dimensional",
  "mode=0\n            0, 1, 2, 3",
  "size=2\n            Box will be from -size : size",
  "npix=32\n           Number of pixels along box vertex",
  "times=all\n         Which times to use from the snapshot [not yet]",
  "body=\n             If set, only select this body for the orbit (0=first)",
  "VERSION=0.4\n       27-dec-2019 XYZ",
  NULL,
};

string usage="Transform an orbit and ODM into 'BOOM' HDF [testing]";

string cvsid="$Id:$";

//https://support.hdfgroup.org/ftp/HDF5/examples/src-html/c.html
// this combines h5_crtdat.c h5_rdwt.c h5_crtatt.c

void scandata(string fin, string fout); // 0
void crtdat(string fout);  // 1
void rdwt(string fin);     // 2
void crtatt(string fin);   // 3


//  make sure we got the 0.5 pixel thing right here , see snapgrid
int xbox(real x, int size, int npix)
{
  if (-size < x && x < size)
    return (int) floor((x-size)/npix);
  return -1;
}

void nemo_main()
{
  string fin = getparam("in");
  string fout = getparam("out");
  int mode = getiparam("mode");

  if (mode==0) {
    scandata(fin,fout);
  } else if (mode==1)
    crtdat(fin);
  else if (mode==2)
    rdwt(fin);
  else if (mode==3)
    crtatt(fin);
}


#define MAXSNAP 10000

void scandata(string fin, string fout)
{
  stream instr = stropen(fin,"r");
  stream outstr = stropen(fout,"w");
  real tsnap, tsnaps[MAXSNAP];
  Body *bp, *btab[MAXSNAP];
  int nbody0 = 0, nbody, bits;
  real sum;
  int i, j, nsnap = 0;

  for (i=0; i<MAXSNAP; i++) btab[i] = NULL;
  
  dprintf(0,"scanning %s\n",fin);
  get_history(instr);
  for (;;) {
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag))
      break;
    get_snap(instr, &btab[nsnap], &nbody, &tsnap, &bits);
    if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0)
      continue;       /* just skip it's probably a diagnostics */
    if ((bits & TimeBit) == 0)
      tsnap = 0.0;
    if (nbody0==0) nbody0 = nbody;
    if (nbody0 != nbody)
      error("Cannot process snapshots with changing nbody (%d,%d)",nbody0,nbody);
    dprintf(1,"Processing %g\n",tsnap);
    tsnaps[nsnap] = tsnap;
    nsnap++;
    if (nsnap == MAXSNAP) {
      warning("Can only handle %d snapshots",MAXSNAP);
      break;
    }
  }

  dprintf(0,"Found %d snapshots with %d bodies each\n",nsnap,nbody0);

  put_history(outstr);
  orbitptr optr;
  allocate_orbit( &optr, 3, nsnap);

  int npix = getiparam("npix");
  real size = getrparam("size");

  imageptr iptr;
  int ix,iy,iz, iorb = -1;
  real dmin, dmax;
  create_cube (&iptr, npix, npix, npix);
  Dx(iptr) = Dy(iptr) = Dz(iptr) = size / npix;
  Xmin(iptr) = Ymin(iptr) = Zmin(iptr) = -size;
  if (hasvalue("body")) {
    iorb = getiparam("body");
    dprintf(0,"Only processing body %d\n",iorb);
  }
  
  for (i=0; i<nbody0; i++) {    // loop over all bodies

    if (iorb >= 0 && iorb != i) continue;
    
    for(ix=0; ix<npix; ix++)           // clear the image
    for(iy=0; iy<npix; iy++)
    for(iz=0; iz<npix; iz++)
      CubeValue(iptr,ix,iy,iz) = 0.0;


    for (j=0; j<nsnap; j++) {          // assemble an orbit
      Torb(optr, j) = tsnaps[j];
      Xorb(optr, j) = Pos(btab[j])[0];
      Yorb(optr, j) = Pos(btab[j])[1];
      Zorb(optr, j) = Pos(btab[j])[2];
      Uorb(optr, j) = Vel(btab[j])[0];
      Vorb(optr, j) = Vel(btab[j])[1];
      Worb(optr, j) = Vel(btab[j])[2];

      ix = xbox(Pos(btab[j])[0],  size, npix);     // simple ODM
      if (ix<0) continue;
      iy = xbox(Pos(btab[j])[1],  size, npix);
      if (iy<0) continue;
      iz = xbox(Pos(btab[j])[2],  size, npix);
      if (iz<0) continue;
      CubeValue(iptr,ix,iy,iz) += 1.0;   // could imagine a better normalization
    }
    dmin = dmax = CubeValue(iptr,0,0,0);
    for(ix=0; ix<npix; ix++)           // clear the image
    for(iy=0; iy<npix; iy++)
      for(iz=0; iz<npix; iz++) {
	if (CubeValue(iptr,ix,iy,iz) < dmin) dmin = CubeValue(iptr,ix,iy,iz);
	if (CubeValue(iptr,ix,iy,iz) > dmax) dmax = CubeValue(iptr,ix,iy,iz);
      }
    MapMin(iptr) = dmin;
    MapMax(iptr) = dmax;

    write_orbit(outstr,optr);
    write_image(outstr,iptr);
  }
  strclose(outstr);
}




#if 0
GROUP "/" {
   DATASET "BOOM" {
      DATATYPE  H5T_STRING {
         STRSIZE H5T_VARIABLE;
         STRPAD H5T_STR_NULLTERM;
         CSET H5T_CSET_UTF8;
         CTYPE H5T_C_S1;
      }
      DATASPACE  SCALAR
      DATA {
      (0): "0.3"
      }
   }
   DATASET "Np" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SCALAR
      DATA {
      (0): 2
      }
   }
   GROUP "particle0" {
   ....
#endif     


void crtdat(string fout)
{  

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



void crtdat1(string fout)
{  

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


void rdwt(string fin)
{
  
   hid_t       file_id, dataset_id;  /* identifiers */
   herr_t      status;
   int         i, j, dset_data[4][6];

   /* Initialize the dataset. */
   for (i = 0; i < 4; i++)
      for (j = 0; j < 6; j++)
         dset_data[i][j] = i * 6 + j + 1;

   /* Open an existing file. */
   file_id = H5Fopen(fin, H5F_ACC_RDWR, H5P_DEFAULT);

   /* Open an existing dataset. */
   dataset_id = H5Dopen2(file_id, "/dset", H5P_DEFAULT);

   /* Write the dataset. */
   status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                     dset_data);

   status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    dset_data);

   /* Close the dataset. */
   status = H5Dclose(dataset_id);

   /* Close the file. */
   status = H5Fclose(file_id);
}


void crtatt(string fin)
{
   hid_t       file_id, dataset_id, attribute_id, dataspace_id;  /* identifiers */
   hsize_t     dims;
   int         attr_data[2];
   herr_t      status;

   /* Initialize the attribute data. */
   attr_data[0] = 100;
   attr_data[1] = 200;

   /* Open an existing file. */
   file_id = H5Fopen(fin, H5F_ACC_RDWR, H5P_DEFAULT);

   /* Open an existing dataset. */
   dataset_id = H5Dopen2(file_id, "/dset", H5P_DEFAULT);

   /* Create the data space for the attribute. */
   dims = 2;
   dataspace_id = H5Screate_simple(1, &dims, NULL);

   /* Create a dataset attribute. */
   attribute_id = H5Acreate2 (dataset_id, "Units", H5T_STD_I32BE, dataspace_id, 
                             H5P_DEFAULT, H5P_DEFAULT);

   /* Write the attribute data. */
   status = H5Awrite(attribute_id, H5T_NATIVE_INT, attr_data);

   /* Close the attribute. */
   status = H5Aclose(attribute_id);

   /* Close the dataspace. */
   status = H5Sclose(dataspace_id);

   /* Close to the dataset. */
   status = H5Dclose(dataset_id);

   /* Close the file. */
   status = H5Fclose(file_id);
}
