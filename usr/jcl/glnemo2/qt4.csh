#!/bin/csh -f

setenv QTDIR /r5data/jcl/local/qt-latest
setenv LD_LIBRARY_PATH ${QTDIR}/lib:${LD_LIBRARY_PATH}
setenv PATH ${QTDIR}/bin:${PATH}
rehash
#
