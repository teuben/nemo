setenv ZENOPATH         "$NEMO/usr/zeno/zeno"
setenv ZCC              "gcc"
#setenv ZCCFLAGS        "-DMACOSX -I$ZENOPATH/inc"
setenv ZCCFLAGS         "-DLINUX  -I$ZENOPATH/inc"
setenv ZLDFLAGS         "-L$ZENOPATH/lib"
setenv ZENO_SAFE_SELECT "true"
setenv ZENO_MSG_OPTION  "all"
#
set path=($ZENOPATH/bin $path)
rehash
