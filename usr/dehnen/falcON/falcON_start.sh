#
# FALCON bash startup, needs to be sourced in the shell
# Written for NEMO V4.1
#
# $FALCON must be set.
#
if [ -z $FALCON ]; then
  echo "Environment variable FALCON is not set, falcON cannot be set up"
  echo "add something like: export FALCON=/xxx/yyy"
  echo "to denote the root directory of falcON to your .bashrc file"
else
  # set enviroment variables FALCONLIB and ACCPATH
  export  FALCONLIB=$FALCON/lib
  export WDUTILSLIB=$FALCON/utils/lib

  if [ -z $ACCPATH ]; then
    export ACCPATH=$FALCON/acc/
  else    
    export ACCPATH=$FALCON/acc/:$ACCPATH
  fi

  # add FALCON/bin to PATH
  export PATH=$FALCON/bin:$PATH

  # add FALCONLIB to LD_LIBRARY_PATH (linux) and DYLD_LIBRARY_PATH (darwin)
  if [ -z "$LD_LIBRARY_PATH" ]; then
    export LD_LIBRARY_PATH=$FALCONLIB:$WDUTILSLIB
  else
    export LD_LIBRARY_PATH=$FALCONLIB:$WDUTILSLIB:$LD_LIBRARY_PATH
  fi
  if [ -z "$DYLD_LIBRARY_PATH" ]; then
    export DYLD_LIBRARY_PATH=$FALCONLIB:$WDUTILSLIB
  else
    export DYLD_LIBRARY_PATH=$FALCONLIB:$WDUTILSLIB:$DYLD_LIBRARY_PATH
  fi
fi
