/*				Play.c					*/
/*	This program will play sound files in u-law format on the 	*/
/*	Sparcstation I.  It is a modification of Sun's sound.c program	*/
/*	distributed with SunOS.  These modifications were done by	*/
/*	Richard Gopstein (..!rutgers!soleil!gopstein), September 1989	*/

#include <sys/types.h>
#include <sys/time.h>
#include <sys/dir.h>
#include <sys/file.h>

#include <fcntl.h>
#include <errno.h>
#include <stdio.h>

#include <sys/ioctl.h>
#include <sbusdev/audioreg.h>
#include <sun/audioio.h>

#define WRITE_SIZE			4096

#define BUFFER_SIZE	819200

#define NOTIFIER_LOOP_DELAY	50000

int	Audio_fd;
int	fd;
int	Compress = 1;
int     play_flag = 1;

#define MODULATION_U_LAW	1
#define MODULATION_A_LAW	2
#define MODULATION_BINARY	3

struct sound_buffer {
	int	bits_per_sample;
	int	samples_per_second;
	int	channels;
	int	modulation;
	int	allocated;
	int	size;
	int	cursor_1;
	int	cursor_2;
	int	cursor_start;
	int	cursor_end;
	int	display_position;
	int	last_displayed;
	int	io_position;
	char	directory[256];
	char	filename[256];
	char	data[BUFFER_SIZE];
} Buffer;


main(argc,argv)
	int	argc;
	char	**argv;
{
        int volume;

	if ((argc != 2) && (argc != 4)) {
	  fprintf(stderr, "Usage: play [-v volume] soundfile\n");
		exit(1);
	}

	argv++;

	if (argc == 4) {

	  if (strcmp(*argv++, "-v") != 0) {
	    fprintf(stderr, "Usage: play [-v volume] soundfile\n");
	    exit(1);
	  } else {
	    volume = atoi(*argv++) * 3;
	  }

	} else {

          volume = 21;
	}

	buffer_initialize(&Buffer);

        file_open(*argv);

	audio_open();

        main_play_volume_proc(volume);

        while(play_flag) {
          main_play();
		nap(NOTIFIER_LOOP_DELAY);
	}
}

/*  This is the routine that is called from the main loop that writes
 *  sound to the device.
 *
 *  It needs the following state data:
	- cursor start and end (calculate in proc)
 *	- current position in the buffer
 */

main_play_volume_proc(value)
	int	value;
{
	int	gr;
	int	ger;

	ger = value/2;
	gr = value - ger;
	if (ger < -10) {
		ger = -10;
		gr = value - ger;
	}
	if (gr > 12) {
		gr = 12;
		ger = value - gr;
	}

	if ((gr + ger) != value) {
		fprintf(stderr,
		"Error: gr (%d) + ger (%d) != value (%d)\n",
		    gr,ger,value);
	}

	audio_set_gr(gr);
	audio_set_ger(ger);
}

main_play()

{
	int	bytes;
	int	start;
	int	end;
	int	queue_size;
	int	display_position;
	int	rtn;

	/*  Calculate the cursor positions every time we enter this loop.
	 *  This will let us change the cursors on the fly and have them
	 *  take effect immediately.
	 */
	if (Buffer.cursor_1 < Buffer.cursor_2) {
		Buffer.cursor_start = Buffer.cursor_1;
		Buffer.cursor_end = Buffer.cursor_2;
	} else {
		Buffer.cursor_start = Buffer.cursor_2;
		Buffer.cursor_end = Buffer.cursor_1;
	}

	/*  See if there is enough room in the device for another block
	 *  of sound.  We want to write in large chunks and keep as much
	 *  data enqueued as possible.
	 */
		if (ioctl(Audio_fd,AUDIOWRITEQ,&bytes) < 0) {
			perror("AUDIOWRITEQ ioctl failed");
		}

		if (ioctl(Audio_fd,AUDIOGETQSIZE,&queue_size) < 0) {
			perror("AUDIOGETQSIZE ioctl failed");
		}

	start = Buffer.io_position;
	end = start + WRITE_SIZE;
	if (((queue_size - bytes) > WRITE_SIZE) &&
	    (start < Buffer.cursor_end)) {
		if (end > Buffer.cursor_end) {
			end = Buffer.cursor_end;
		}
		Buffer.io_position = end;

			rtn = write(Audio_fd,&Buffer.data[start],
			    end - start);
			if (rtn < 0) {
		  perror("audio write");
			}
		}

		if (ioctl(Audio_fd,AUDIOWRITEQ,&bytes) < 0) {
			perror("AUDIOWRITEQ ioctl failed");
		}


	/*  XXX We don't want to return until we are done playing the sound.
	 */

	if ((end >= Buffer.cursor_end) && (bytes == 0)) {


			main_play_off();
		}
	}


/*  This routine is used when we stop playing sound, for whatever reason.
 */
main_play_off()

{
	int	dummy;

		if (ioctl(Audio_fd,AUDIOSTOP,&dummy) < 0) {
			perror("AUDIOSTOP ioctl failed");
		}

        play_flag = 0;

}

file_open(path)
char *path;
{


	if ((fd = open(path,0)) < 0) {
		perror("Error opening input file");
		return;
	}

	Buffer.size = read(fd,Buffer.data,BUFFER_SIZE);
	Buffer.display_position = 0;
	Buffer.samples_per_second = 8192;
	Buffer.bits_per_sample = 8;
	Buffer.channels = 1;
	Buffer.modulation = MODULATION_U_LAW;
	Buffer.cursor_1 = 0;
	Buffer.cursor_2 = Buffer.size - 1;
	Compress = 1;

}

audio_open()
{
	struct	audio_ioctl	tmp;

	if ((Audio_fd = open("/dev/audio",2)) < 0) {
		perror("/dev/audio");
		exit(1);
	}

	/* Initialize the three gain registers to 0db. */
	audio_set_ger(0);
	audio_set_gr(0);
	audio_set_gx(0);

	/*
	 * Tell the chip to set the gains according to the
	 * register values we just set.
	 */

	tmp.control = AUDIO_MAP_MMR1;
	tmp.data[0] = AUDIO_MMR1_BITS_LOAD_GX |
	    AUDIO_MMR1_BITS_LOAD_GR |
	    AUDIO_MMR1_BITS_LOAD_GER;

		if (ioctl(Audio_fd,AUDIOSETREG,&tmp) < 0) {
			perror("Set REG");
		}

	/*  Initialize the MMR2 register to send the output to the builtin
	 *  speaker.  This is the default, but we will be defensive
	 *  anyway.  Note that we read the register and turn on the
	 *  appropriate bit, rather thanjust setting it.  This keeps us
	 *  from stompin any useful bits that are already set.  It turns
	 *  out that there is one -- the bit that selects AINB for
	 *  recording.
	 */
		tmp.control = AUDIO_MAP_MMR2;
		if (ioctl(Audio_fd,AUDIOGETREG,&tmp) < 0) {
			perror("Set REG");
		}

		tmp.data[0] |= AUDIO_MMR2_BITS_LS;
		if (ioctl(Audio_fd,AUDIOSETREG,&tmp) < 0) {
			perror("Set REG");
		}
	}


	/*	These are tables of values to be loaded into various
		gain registers.

		Note that for the ger entry for -8db, we use the data
		sheet value for -7.5db.  The data sheet gives values for
		-8db which are wrong and produce too much gain.
	*/

static	unsigned char ger_table[][2] = {
		0xaa,	0xaa,	/* -10db */
		0x79,	0xac,
		/*0x41,	0x91,*/
		0x31,	0x99,	/* -7.5db */
		0x9c,	0xde,
		0x74,	0x9c,	/* -6db */
		0x6a,	0xae,
		0xab,	0xdf,
		0x64,	0xab,
		0x2a,	0xbd,
		0x5c,	0xce,
		0x00,	0x99,	/* 0db */
		0x43,	0xdd,
		0x52,	0xef,
		0x55,	0x42,
		0x31,	0xdd,
		0x43,	0x1f,
		0x40,	0xdd,	/* 6db */
		0x44,	0x0f,
		0x31,	0x1f,
		0x10,	0xdd,
		0x41,	0x0f,
		0x60,	0x0b,
		0x42,	0x10,	/* 12db */
		0x11,	0x0f,
		0x72,	0x00,
		0x21,	0x10,
		0x22,	0x00,
		0x00,	0x0b,
		0x00,	0x0f,	/* 18db */
};


static	unsigned char gr_gx_table[][2] = {
		0x8b,	0x7c,	/* -18db */
		0x8b,	0x35,
		0x8b,	0x24,
		0x91,	0x23,
		0x91,	0x2a,
		0x91,	0x3b,
		0x91,	0xf9,	/* -12db */
		0x91,	0xb6,
		0x91,	0xa4,
		0x92,	0x32,
		0x92,	0xaa,
		0x93,	0xb3,
		0x9f,	0x91,	/* -6db */
		0x9b,	0xf9,
		0x9a,	0x4a,
		0xa2,	0xa2,
		0xaa,	0xa3,
		0xbb,	0x52,
		0x08,	0x08,	/* 0db */
		0x3d,	0xac,
		0x25,	0x33,
		0x21,	0x22,
		0x12,	0xa2,
		0x11,	0x3b,
		0x10,	0xf2,	/* 6db */
		0x02,	0xca,
		0x01,	0x5a,
		0x01,	0x12,
		0x00,	0x32,
		0x00,	0x13,
		0x00,	0x0e,	/* 12db */
};



audio_set_ger(value)
	int	value;
{
	struct	audio_ioctl	tmp;


	if ((value < -10) || (value > 18)) {
		fprintf(stderr,
		    "GER value %d out of range; %d <= GER <=  %d\n",
		    value,0,18);
		return;
	}

	/*  Add 10 to the value to get the index into the table.
	 */
	tmp.control = AUDIO_MAP_GER;
	tmp.data[0] = ger_table[value + 10][1];
	tmp.data[1] = ger_table[value + 10][0];

		if (ioctl(Audio_fd,AUDIOSETREG,&tmp) < 0) {
			perror("Set REG");
		}
	}


audio_set_gr(value)
	int	value;
{
	struct	audio_ioctl	tmp;

	if ((value < -18) || (value > 12)) {
		fprintf(stderr,
		    "GR value %d out of range; %d <= GR <=  %d\n",
		    value,0,12);
		return;
	}

	tmp.control = AUDIO_MAP_GR;
	tmp.data[0] = gr_gx_table[value + 18][1];
	tmp.data[1] = gr_gx_table[value + 18][0];


		if (ioctl(Audio_fd,AUDIOSETREG,&tmp) < 0) {
			perror("Set REG");
		}

	}


audio_set_gx(value)
	int	value;
{
	struct	audio_ioctl	tmp;

	if ((value < -18) || (value > 12)) {
		fprintf(stderr, 
		    "GX value %d out of range; %d <= GX <=  %d\n",
		    value,0,12);
		return;
	}

	/*  We add 18 to get the index into the table, since entry 0 represents
	 *  -18db.
	 */
	tmp.control = AUDIO_MAP_GX;
	tmp.data[0] = gr_gx_table[value + 18][1];
	tmp.data[1] = gr_gx_table[value + 18][0];

		if (ioctl(Audio_fd,AUDIOSETREG,&tmp) < 0) {
			perror("Set REG");
		}
	}

buffer_initialize(buffer)
	struct	sound_buffer	*buffer;
{
	extern char *getwd();

	buffer->size = 0;
	buffer->filename[0] = 0;
	buffer->cursor_1 = 0;
	buffer->cursor_2 = 0;

	bzero(buffer->data, buffer->size);

	if (!getwd(Buffer.directory))
		strcpy(Buffer.directory, ".");
}


static
nap(usec)
	int usec;
{
	struct timeval timeout;

	timeout.tv_sec = 0;
	timeout.tv_usec = usec;

	(void) select(0, (fd_set *) 0, (fd_set *) 0, (fd_set *) 0,
		&timeout);
}


