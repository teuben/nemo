#!/bin/csh -f
#
#
#

set DESTDIR = compile

set dirs = ( ${DESTDIR}/obj ${DESTDIR}/obj-debug ${DESTDIR}/bin )

if ( ! -d $DESTDIR ) then
    mkdir $DESTDIR
endif

foreach i ( $dirs )
    if ( ! -d $i ) then
	mkdir $i
    endif
end
#
