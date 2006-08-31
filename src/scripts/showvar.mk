# show a NEMO makedef variable, return 0 if not present
# normally this makefile is used by the shell script 'showvar'

include  $(NEMOLIB)/makedefs

# default
VAR = NEMO

var:
ifneq ($(strip $($(VAR))),)
	@echo $($(VAR))
else
	@echo 0
endif
