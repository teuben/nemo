# setup manybody environment for csh

#   you really should hardcode your absolute address in here
setenv MANYBODY `pwd`

setenv ACSROOT                $MANYBODY/acs

setenv NEMO                   $MANYBODY/nemo_cvs
if (-e $NEMO/nemo_start.csh) then
   source $NEMO/nemo_start.csh
endif

setenv STARLAB                $MANYBODY/starlab_cvs
if (-e $STARLAB/starlab_start.csh) then
   source $STARLAB/starlab_start.csh
endif

if (-e $MANYBODY/intel/bin/iccvars.csh) then
   source $MANYBODY/intel/bin/iccvars.csh
endif


set path=(. $MANYBODY/opt/bin $path)
rehash
