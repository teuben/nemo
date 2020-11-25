/*
 *   LMTINFO:     some info from the LMT netCDF4 data files
 *
 */

#include <nemo.h>
#include <netcdf.h>


string defv[] = {
    "in=???\n            Input netCDF4 file",
    "VERSION=0.1\n       24-nov-2020 PJT",
    NULL,
};

string usage = "netCDF4 info and bench reduction procedures";



#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


void nemo_main(void)
{
  string infile;
  int ncid, retval;

  infile = getparam("in");

  if ((retval = nc_open(infile, NC_NOWRITE, &ncid)))
    ERR(retval);



  
#if 0
  // read the first SpecFile 
  read_spec_file(&S, OTF.i_filename[0]);

  // copy over obs header variables
  C.obsnum = S.obsnum;
  printf("%d\n",C.obsnum);
  printf("%s\n",S.source);
  strncpy(C.source,S.source,18);
  printf("%s\n",C.source);
  C.x_position = S.x_position;
  C.y_position = S.y_position;
  // set up convolution array for the gridding process.
  initialize_convolve_function(&CF, OTF.resolution_size, OTF.cell_size, OTF.rmax, OTF.nsamples);
  printf("4\n");
  if(OTF.otf_select == 1)
    initialize_jinc_filter(&CF, OTF.otf_jinc_a, OTF.otf_jinc_b, OTF.otf_jinc_c);
  else if(OTF.otf_select == 2)
    initialize_gauss_filter(&CF, OTF.otf_jinc_b);
  else
    initialize_box_filter(&CF, OTF.cell_size/2.);

  /* prints the convolution function
  printf("n_cell= %d\n",CF.n_cells);
  printf("r (as)   c\n");
  for(i=0;i<256;i++)
    printf("%5.2f %8.4f\n",i*CF.delta, CF.array[i]);
  */

  // initialize cube and axes
  n[0] = 2 * (int)(floor((OTF.x_extent+OTF.cell_size/2.)/OTF.cell_size)) + 1;
  n[1] = 2 * (int)(floor((OTF.y_extent+OTF.cell_size/2.)/OTF.cell_size)) + 1;
  n[2] = S.nchan;

  printf("5\n");
  initialize_cube(&C, n);
  // note that we add one to crpix's per fits convention
  initialize_cube_axis(&C, Z_AXIS, S.CRVAL, S.CRPIX+1., S.CDELT, S.CTYPE, "km/s");
  initialize_cube_axis(&C, X_AXIS, 0.0, (n[0]-1.)/2.+1., OTF.cell_size, "X", "arcsec");
  initialize_cube_axis(&C, Y_AXIS, 0.0, (n[1]-1.)/2.+1., OTF.cell_size, "Y", "arcsec");

  initialize_plane(&Weight, n);
  initialize_plane_axis(&Weight, PLANE_X_AXIS, 0.0, (n[0]-1.)/2.+1., OTF.cell_size, "X", "arcsec");
  initialize_plane_axis(&Weight, PLANE_Y_AXIS, 0.0, (n[1]-1.)/2.+1., OTF.cell_size, "Y", "arcsec");
  printf("6\n");

  free_spec_file(&S);
  printf("axes initialized\n");

  for(i=0;i<16;i++)
    printf("%d %d\n",i,OTF.use_pixels[i]);

  for(ifile=0;ifile<OTF.nfiles;ifile++)
    {
      // read the new specfile for gridding
      printf("file %d %s\n",ifile,OTF.i_filename[ifile]);
      read_spec_file(&S, OTF.i_filename[ifile]);

      // now we do the gridding
      for(i=0;i<S.nspec;i++)
	{
	  if(OTF.use_pixels[S.Pixel[i]] == 1)
	    {
	      if(S.RMS[i] < OTF.rms_cutoff)
		{
		  spectrum = get_spectrum(&S,i);
		  ix = cube_axis_index(&C, X_AXIS, S.XPos[i]);
		  iy = cube_axis_index(&C, Y_AXIS, S.YPos[i]);
		  if( (ix>=0) && (iy>=0) )
		    {
		      for(ii=-CF.n_cells; ii<=CF.n_cells; ii++)
			for(jj=-CF.n_cells; jj<=CF.n_cells; jj++)
			  {
			    x = S.XPos[i]-C.caxis[X_AXIS][ix+ii];
			    y = S.YPos[i]-C.caxis[Y_AXIS][iy+jj];
			    distance = sqrt(x*x+y*y);
			    if ((S.RMS[i] != 0.0) && (OTF.noise_sigma > 0.0)) {
			      rmsweight = 1.0 /(S.RMS[i] * S.RMS[i]);
			    } else {
			      rmsweight = 1.0;
			    }
			    weight = get_weight(&CF, distance) * rmsweight;
			    iz = cube_z_index(&C, C.caxis[X_AXIS][ix+ii], C.caxis[Y_AXIS][iy+jj]);
			    for(k=0;k<C.n[Z_AXIS];k++)
			      C.cube[iz+k] = C.cube[iz+k] + weight * spectrum[k];
			    izp = plane_index(&Weight, C.caxis[X_AXIS][ix+ii], C.caxis[Y_AXIS][iy+jj]);
			    Weight.plane[izp] = Weight.plane[izp] + weight;
			  }		    
		    }
		}
	    }
	}
      free_spec_file(&S);
    }

  printf("Cube Completed\n");

  // compute averages for each map point; if no data assign NAN
  for(i=0;i<C.n[X_AXIS];i++)
    {
      x = C.caxis[X_AXIS][i];
      for(j=0;j<C.n[Y_AXIS];j++)
	{
	  y = C.caxis[Y_AXIS][j];
	  izp = plane_index(&Weight, x, y);
	  iz = cube_z_index(&C, x, y);
	  if(Weight.plane[izp]>0.)
	    {
	      for(k=0;k<C.n[Z_AXIS];k++)
		C.cube[iz+k] = C.cube[iz+k] / Weight.plane[izp];
	    }
	  else
	    {
	      for(k=0;k<C.n[Z_AXIS];k++)
		C.cube[iz+k] = NAN;
	    }
	}
    }

  printf("Weighting Completed\n");

  // dumping the spectrum at 0,0 for fun... 
  izp = plane_index(&Weight, 0.0, 0.0);
  printf("Weight of %f %f is %f\n",0.0,0.0,Weight.plane[izp]);
  iz = cube_z_index(&C, 0.0, 0.0);
  for(i=0;i<S.nchan;i++)
    printf("%d %8.3f %6.2f\n ",i, C.caxis[Z_AXIS][i],C.cube[iz+i]);

  printf("write to %s\n",OTF.o_filename);
  // finally write the data cube as FITS file
  write_fits_cube(&C, OTF.o_filename);

#endif  
}
