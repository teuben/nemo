#

PDIRS = `find . -type d -name .git`

# the Makefile tree in the $NEMO/usr tree is not as well structured as
#     the one in $NEMO/src

help:

install:

libs:

bins:

clean:

pull:
	@for i in $(PDIRS);do\
	(cd $$i/..; echo "===== `pwd`"; git pull); done

