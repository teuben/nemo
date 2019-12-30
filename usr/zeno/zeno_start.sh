#   this is the ZENO rc file assuming you have installed 
#   this via NEMO in $NEMO/usr/zeno/zeno

export ZENOPATH="$NEMO/usr/zeno/zeno"
export ZCC="gcc"
if test $NEMOSYS = "Linux"; then
    export ZCCFLAGS="-std=gnu99 -DLINUX  -I$ZENOPATH/inc"
else
    export ZCCFLAGS="-std=gnu99 -DMACOSX -I$ZENOPATH/inc"
fi

export ZLDFLAGS="-L$ZENOPATH/lib"
export ZENO_SAFE_SELECT="true"
export ZENO_MSG_OPTION="all"


export PATH=$ZENOPATH/bin:$PATH
