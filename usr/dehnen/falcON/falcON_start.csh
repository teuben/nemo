#
# this file is normally sourced from falcON_start
#
# $FALCON must be set.
#
if ($?FALCON == 0) then
  echo "Environment variable FALCON is not set, falcON cannot be set up"
  echo "add something like: setenv FALCON /xxx/yyy"
  echo "to denote the root directory of falcON to your .cshrc file"
  goto done
endif

# set enviroment variables FALCONLIB and ACCPATH
setenv FALCONLIB $FALCON/lib
setenv WDUTILSLIB $FALCON/utils/lib

if($?ACCPATH) then
    setenv ACCPATH $FALCON"/acc/:"$ACCPATH
else
    setenv ACCPATH $FALCON"/acc/"
endif

# add FALCON/bin to PATH
setenv PATH $FALCON"/bin:"$PATH
rehash

if ($?prompt) then
  echo "falcON loaded from $FALCON"
endif

# add FALCONLIB to LD_LIBRARY_PATH and DYLD_LIBRARY_PATH
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${FALCONLIB}":"${WDUTILSLIB}":"${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${FALCONLIB}":"${WDUTILSLIB}
endif

set uname = (`uname -s`)       # for now, we only support Linux and Darwin
if($uname == Darwin) then
    if($?DYLD_LIBRARY_PATH) then
      setenv DYLD_LIBRARY_PATH ${FALCONLIB}":"${WDUTILSLIB}":"${DYLD_LIBRARY_PATH}
    else
      setenv DYLD_LIBRARY_PATH $LD_LIBRARY_PATH
    endif
endif

done:

