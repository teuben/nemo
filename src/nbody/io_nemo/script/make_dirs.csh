#!/bin/csh -f
#
#
#

set ARCH = ${OSTYPE}

set dirs = ( ${ARCH}/obj ${ARCH}/obj-debug ${ARCH}/bin )

if ( ! -d $ARCH ) then
    mkdir $ARCH
endif

foreach i ( $dirs )
    if ( ! -d $i ) then
	mkdir $i
    endif
end
#
