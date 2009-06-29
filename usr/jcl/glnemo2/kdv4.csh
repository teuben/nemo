#!/bin/csh -f

setenv KDEVELOP /r5data/jcl/local/latest-kdevelop
setenv KDEDIRS /r5data/jcl/local/latest-kdevelop:/usr/bin
kbuildsycoca
setenv QTDIR /r5data/jcl/local/qt-latest
setenv PATH ${QTDIR}/bin:${KDEVELOP}/bin:${PATH}
setenv LD_LIBRARY_PATH ${QTDIR}/lib:${LD_LIBRARY_PATH}
rehash
exec kdevelop $*

#
