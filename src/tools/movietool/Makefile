# Makefile for movietool - slightly adapted for NEMO's hierarchical Make calls
# Location:  $NEMO/src/tools/movietool
# cc -> $(CC)
AUDIODIR = `pwd`
AUDIOPLAYER = play
BINDIR = $(NEMOBIN)
MANEXT = 1
MANDIR = $(NEMO)/man/man${MANEXT}
CC = /bin/cc -Bdynamic
# END OF LOCAL NEMO STUFF

LIBS = -lsuntool -lsunwindow -lpixrect
# For profiling:
# LIBS = -lsuntool_p -lsunwindow_p -lpixrect_p

# Change the AUDIOPLAYER and AUDIODIR definitions,
# or comment out the next 2 lines if you don't have audio 
#AUDIOPLAYER = /usr/local/bin/play
#AUDIODIR = /usr/local/common/audio

# Change the CDPLAYER definition,
# or comment out the next line if you don't have a CD-ROM drive attached.
# If you do have it, do "cd cdrom;make install" to get a CD-player.
CDPLAYER = /usr/local/bin/movietool_cdplayer

# PROFILE_FLAG=-pg
CFLAGS = -O ${PROFILE_FLAG} -DAUDIOPLAYER=\"${AUDIOPLAYER}\" -DAUDIODIR=\"${AUDIODIR}\" -DCDPLAYER=\"${CDPLAYER}\"

# Installation directories etc.
#BINDIR = /usr/local/bin
#MANEXT = l
#MANDIR = /usr/man/man${MANEXT}

OBJS = main.o frames.o zoom.o pr_load.o pr_decode.o pw_rop.o paneledit.o \
	audio.o cdrom.o

movietool: movietool.o movietool.h ${OBJS}
	$(CC) ${CFLAGS} -o movietool movietool.o ${OBJS} ${LIBS}

play: play.c
	$(CC) ${CFLAGS} -o play play.c

pw_decode.o: pw_decode.c
	${CC} -O4 ${PROFILE_FLAG} -c pw_decode.c

pr_decode.o: pr_decode.c
	${CC} -O4 ${PROFILE_FLAG} -c pr_decode.c

pr_load.o: pr_load.c
	${CC} -O4 ${PROFILE_FLAG} -c pr_load.c

zoom.o: zoom.c
	${CC} -O4 ${PROFILE_FLAG} -c zoom.c

${OBJS}: movietool.h

movietool.o: movietool.icon

install: movietool
	install -s movietool ${BINDIR}
	cp movietool.man ${MANDIR}/movietool.${MANEXT}

clean:
	rm -f core *.o movietool
