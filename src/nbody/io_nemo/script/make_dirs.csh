#!/bin/csh -f
#
# $Id$
#

set ARCH = ${OSTYPE}

set dirs = ( ${ARCH}/obj ${ARCH}/bin )

if ( ! -d $ARCH ) then
    mkdir $ARCH
endif

foreach i ( $dirs )
    if ( ! -d $i ) then
	mkdir $i
    endif
end
#
