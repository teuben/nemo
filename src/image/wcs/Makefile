#  Makefile for testpos.c and worldpos.c 
#  Copyright (C) 1994
#  Associated Universities, Inc. Washington DC, USA.
# 
#  This library is free software; you can redistribute it and/or modify it
#  under the terms of the GNU Library General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or (at your
#  option) any later version.
# 
#  This library is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
#  License for more details.
# 
#  You should have received a copy of the GNU Library General Public License
#  along with this library; if not, write to the Free Software Foundation,
#  Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
# 
#  Correspondence concerning AIPS should be addressed as follows:
#         Internet email: aipsmail@nrao.edu
#         Postal address: AIPS Group
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA

# Substitute your favorite ANSI-C compiler (e.g., 'acc' for Suns):
#CC =	gcc -pedantic -lm
CC  =	acc

all:	runtest

worldpos.o:	worldpos.c
	${CC} -c worldpos.c

testpos:	worldpos.o testpos.c
	${CC}  testpos.c worldpos.o -o ${@}

testpos.out:	testpos
	testpos > ${@}

runtest:	testpos.out
	more testpos.out
	@echo " "
	@echo "Compare test run to reference output (if present):"
	-if [ -f testpos.ref ] ;then diff testpos.out testpos.ref ; fi

ITEMS   = worldpos.c testpos.c testpos.ref COPYING.LIB-2.0 Makefile
PACKAGE = worldpos.tar
ANONFTP = /fits_1/ftp/fits/src/wcs/

export:	testpos.out
	rm -Rf worldpos
	mkdir worldpos
	cp -p testpos.out testpos.ref
	cp -p ${ITEMS} worldpos
	rm -f ${PACKAGE}* testpos.ref
	tar cf ${PACKAGE} worldpos
	rm -Rf worldpos
	tar tvf ${PACKAGE}
	gzip ${PACKAGE}
	ls -lg ${PACKAGE}.gz
	cp -p ${PACKAGE}.gz ${ANONFTP}
	ls -lg ${ANONFTP}

clean:	
	rm -f *~ core *.o testpos testpos.out \
		${PACKAGE}* ${ANONFTP}${PACKAGE}.gz
	rm -Rf worldpos

LIST	= testpos.out Makefile testpos.c worldpos.c

print:	${LIST}
	# substitute your favorite print utility for 'lwf.sh':
	lwf.sh ${LIST}
