
#

zeno:  zeno_jeb zeno_pjt
	@echo 'You will need to create a symlink zeno -> zeno_XYZ'
	@echo Pick from:
	@echo ln -s zeno_jeb zeno
	@echo ln -s zeno_pjt zeno

#  original one
zeno_jeb:
	git clone https://github.com/joshuabarnes/zeno/ zeno_jeb

#  my current repo has a branch for linux
zeno_pjt:
	git clone https://github.com/teuben/zeno zeno_pjt
	(cd zeno_pjt; git co teuben1)
