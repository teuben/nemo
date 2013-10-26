SH=/bin/sh

all:
	@echo Nothing to build

.PHONY : cvsdist
cvsdist:
	@( CVSROOT=`cat CVS/Root` && \
	export CVSROOT && \
	VERSION="`grep '^# Version' pipestatus | \
	   sed 's/^[^0-9]*\\([0-9].*\\)[^0-9]*$$/\\1/'`" && \
	VERSION_CVS=`echo $${VERSION} | tr . -` && \
	export VERSION VERSION_CVS && \
	$(MAKE) ChangeLog && \
	cp ChangeLog /tmp/pipestatus.ChangeLog.$${VERSION} && \
	echo "***** Exporting files for pipestatus-$${VERSION}..." && \
	cd /tmp && \
	rm -rf /tmp/pipestatus-$${VERSION} && \
	cvs export -fd pipestatus-$${VERSION} -r pipestatus-$${VERSION_CVS} pipestatus && \
	cd pipestatus-$${VERSION} && \
	mv /tmp/pipestatus.ChangeLog.$${VERSION} ChangeLog && \
	chmod -R a+rX,g-w . && \
	cd .. && \
	echo "***** Taring and zipping pipestatus-$${VERSION}.tar.gz..." && \
	tar cvvf pipestatus-$${VERSION}.tar pipestatus-$${VERSION} && \
	gzip -9f pipestatus-$${VERSION}.tar && \
	echo "***** Done making /tmp/pipestatus-$${VERSION}.tar.gz")

.PHONY: ChangeLog
ChangeLog:
	@(echo "***** Making new ChangeLog..."; \
	rm -f ChangeLog; \
	AWK=gawk rcs2log -i 2 -r -d'2007-09-27<' . \
	| sed -e 's,/cvsroot/runawk/,,g' \
	      -e 's,cheusov@[^>]*,vle@gmx.net,g' \
	      -e 's,author  <[^>]*,Aleksey Cheusov <vle@gmx.net,g' \
              -e 's,\(.*\)<\([^@<>]\+\)@\([^@<>]\+\)>\(.*\),\1<\2 at \3}\4,g' \
	> ChangeLog;)

.PHONY: test
test:
	@printf " ---------------\nSH=%s\n" ${SH}; \
	set -e; ${SH} ./selftest; echo succeeded

SHELLS=/bin/sh /bin/ksh /bin/bash /bin/dash /bin/zsh /bin/mksh /bin/pdksh \
	/usr/bin/sh /usr/bin/ksh /usr/bin/bash \
	/usr/bin/dash /usr/bin/zsh /usr/bin/mksh /usr/bin/pdksh \
	/usr/pkg/bin/ksh /usr/pkg/bin/bash /usr/pkg/bin/dash \
	/usr/pkg/bin/zsh /usr/bin/mksh /usr/pkg/bin/pdksh \
	/usr/pkg/heirloom/bin/sh

.PHONY: test_all
test_all:
	@set -e; \
	for sh in ${SHELLS}; do \
		if test -x $$sh; then \
			$(MAKE) test SH=$$sh; \
		fi; \
	done
