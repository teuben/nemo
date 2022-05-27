/*
 *   LMTINFO:     some info from the LMT netCDF4 data files
 *
 */

#include <nemo.h>
#include <netcdf.h>


string defv[] = {
    "in=???\n            Input netCDF4 file",
    "mode=1\n            1=ifproc  2=roach   3=SpecFile",
    "VERSION=0.2a\n      27-may-2022 PJT",
    NULL,
};

string usage = "netCDF4 info and bench reduction procedures";



#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


void nemo_main(void)
{
  string infile;
  int mode, ncid, retval;

  warning("Experimental program to interrogate LMT netcdf files");

  infile = getparam("in");
  mode = getiparam("mode");

  if ((retval = nc_open(infile, NC_NOWRITE, &ncid)))
    ERR(retval);

  dprintf(0,"File: %s  ncid=%d\n",infile,ncid);

  int obsnum_id, source_id,source_x,source_y,crval_id,crpix_id,cdelt_id,ctype_id,caxis_id;
  int obsnum;

  if (mode == 1) { // IFPROC files
    
    // Header.Source.SourceName = "IRC+10216 
    if ((retval = nc_inq_varid(ncid, "Header.Source.SourceName", &source_id)))
      ERR(retval);
    dprintf(0,"Header.Obs.ObsNum: %d\n",source_id);
    
    if ((retval = nc_inq_varid(ncid, "Header.IfProc.CalObsNum", &obsnum_id)))
      ERR(retval);
    if((retval = nc_get_var_int(ncid,obsnum_id, &obsnum)) != NC_NOERR)
      ERR(retval);

    dprintf(0,"Header.IfProc.CalObsNum: %d %d\n",obsnum_id, obsnum);
    
  } else if (mode == 2) { // ROACH files

    // Header.Telescope.ObsNum = 79447 ;
    if ((retval = nc_inq_varid(ncid, "Header.Telecope.ObsNum", &obsnum_id)))
      ERR(retval);
    if((retval = nc_get_var_int(ncid,obsnum_id, &obsnum)) != NC_NOERR)
      ERR(retval);
    
    dprintf(0,"Header.Telescope.ObsNum: %d %d\n",obsnum_id,obsnum);

    // Header.Mode.numchannels = 2048 ;
    // Header.Telescope.ObsNum = 79447 ;
    // Header.Telescope.source_name = "IRC+10216" ;
    // Header.Telescope.obspgm = "Cal" ;

  } else if (mode == 3) { // SpecFile

    // Header.Telescope.ObsNum = 79447 ;
    if ((retval = nc_inq_varid(ncid, "Header.Obs.ObsNum", &obsnum_id)))
      ERR(retval);
    if((retval = nc_get_var_int(ncid,obsnum_id, &obsnum)) != NC_NOERR)
      ERR(retval);

    // Data.Spectra    float Data.Spectra(nspec, nchan) ;
    int nspec_id, nchan_id;
    size_t nspec, nchan;
    if ((retval = nc_inq_dimid(ncid, "nspec", &nspec_id)))
      ERR(retval);
    if ((retval = nc_inq_dimid(ncid, "nchan", &nchan_id)))
      ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nspec_id, &nspec)))
      ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nchan_id, &nchan)))
      ERR(retval);

    int data_id;
    float *data = (float *)malloc(nspec*nchan*sizeof(float));

    if ((retval = nc_inq_varid(ncid, "Data.Spectra", &data_id)))
      ERR(retval);
    if ((retval = nc_get_var_float(ncid, data_id, data)))
      ERR(retval);
    double psum = 0.0;
    double nsum = 0.0;
    for (int i=0; i<nspec*nchan; i++) {
      if (data[i] < 0.0)
	nsum += data[i];
      else
	psum += data[i];
    }
    double sratio = (psum+nsum)/(psum-nsum);
    
    dprintf(0,"Header.Obs.ObsNum: %d %d\n",obsnum_id,obsnum);
    dprintf(0,"nspec=%d nchan=%d\n",nspec,nchan);
    dprintf(0,"data sumn=%f sump=%f sratio=%f\n",nsum,psum,sratio);
    

  }

}
