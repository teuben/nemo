cc -c -o pg.cleanup.o -DPGDISP -O -I/usr/openwin/include cleanup.c
cc -c -O -I/usr/openwin/include pgdisp.c
cc -c -o pg.figcurs.o -DPGDISP -O -I/usr/openwin/include figcurs.c
cc -c -o pg.getdata.o -DPGDISP -O -I/usr/openwin/include getdata.c
cc -c -o pg.getvisuals.o -DPGDISP -O -I/usr/openwin/include getvisuals.c
cc -c -o pg.handlexevent.o -DPGDISP -O -I/usr/openwin/include handlexevent.c
cc -c -o pg.proccom.o -DPGDISP -O -I/usr/openwin/include proccom.c
cc -c -o pg.resdb.o -DPGDISP -O -I/usr/openwin/include resdb.c
cc -c -O -I/usr/openwin/include exposelgwin.c
cc -c -O -I/usr/openwin/include getcolors.c
cc -c -O -I/usr/openwin/include initlgluts.c
cc -c -O -I/usr/openwin/include initlgwin.c
cc -c -O -I/usr/openwin/include initlock.c
cc -c -O -I/usr/openwin/include initwmattr.c
cc -c -O -I/usr/openwin/include mainloop.c
cc -c -O -I/usr/openwin/include resizelgwin.c
cc -c -O -I/usr/openwin/include returnbuf.c
cc -c -O -I/usr/openwin/include waitevent.c
cc -c -O -I/usr/openwin/include updatelgtitle.c
cc -o pgdisp -O -I/usr/openwin/include -Bstatic pg.cleanup.o pgdisp.o pg.figcurs.o pg.getdata.o pg.getvisuals.o pg.handlexevent.o pg.proccom.o pg.resdb.o exposelgwin.o getcolors.o initlgluts.o initlgwin.o initlock.o initwmattr.o mainloop.o resizelgwin.o returnbuf.o waitevent.o updatelgtitle.o -lX11 -lm
