void audio (audiofile)
char *audiofile;	/* Filename relative to AUDIODIR */
{

/* The AUDIODIR is the directory of *.au audio files */
/* The AUDIOPLAYER is the program that will play audio files */

/* It would be good to test for the existence of /dev/audio.
 * I don't know if there is a faster way than opening the device,
 * which I anyway believe gives an audible click when being opened.
 * For lack of anything better, I just check for AUDIODIR and AUDIOPLAYER */

#if defined AUDIODIR && defined AUDIOPLAYER
	char buf[200];
	sprintf (buf, "%s %s/%s &", AUDIOPLAYER, AUDIODIR, audiofile);
	system (buf);	/* execute the command in the background (& flag) */
#endif
}
