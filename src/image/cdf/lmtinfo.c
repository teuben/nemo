/*
 *   LMTINFO:     some info from the LMT netCDF4 data files
 *
 */

#include <nemo.h>
#include <netcdf.h>


string defv[] = {
    "in=???\n            Input netCDF4 file",
    "mode=1\n            1=ifproc  2=roach",
    "VERSION=0.1\n       24-nov-2020 PJT",
    NULL,
};

string usage = "netCDF4 info and bench reduction procedures";



#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


void nemo_main(void)
{
  string infile;
  int mode, ncid, retval;

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

  }

}
