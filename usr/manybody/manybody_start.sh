# setup manybody environment for sh

# you really should hardcode your address in here
export MANYBODY=`pwd`

export NEMO=$MANYBODY/nemo_cvs
export STARLAB=$MANYBODY/starlab_cvs
export ACSROOT=$MANYBODY/acs


source $NEMO/nemo_start.sh
source $STARLAB/starlab_start.sh

PATH=$MANYBODY/bin:$PATH
