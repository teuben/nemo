rem
rem This is another variant on building with djgpp....
rem
make -f makefile.dj2 DEVICES="-DPOSTSCRIPT -DHPGL -DGRX" MCFLAGS=-g2 DOBJS="../drivers/ps.o ../drivers/hpdxy.o ../drivers/grx.o" RANLIB=ranlib
