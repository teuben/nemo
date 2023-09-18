include $(NEMOLIB)/makedefs


URL1 = https://github.com/jfwallin/JSPAM


help:
	@echo no help yet


JSPAM:
	git clone $(URL1)

clean:
	rm -rf JSPAM

install: JSPAM
	(cd JSPAM/fortran; make ; cp basic_run $(NEMOBIN)/jspam)


# example runs

run1:
	jspam
	(rm a.snap; tabtos a.101 a.snap block1=pos,vel)
