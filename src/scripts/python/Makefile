# 	$NEMO/src/scripts/Makefile

include $(NEMOLIB)/makedefs


TESTS = test_unsio.py test_gala.py test_galpy.py test_pynbody.py test_yt.py test_amuse.py


help:
	@echo special NEMO stuff
	@echo No more specific help here

clean:
	rm -f p100

tests:
	-@for p in $(TESTS); do\
	(echo $$p; ./$$p > /dev/null); done
