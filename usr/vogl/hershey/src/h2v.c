#include <stdio.h>
#include <fcntl.h>
#include "h2v.h"

/*
 *	Read the raw hersh font data and write the binary font files
 *	as specified in h2v.h or as specified in an index file.
 */

extern char	*malloc();

HTAB	hersh[MAX_CHARS];

void	readdata(), readindex(), writefont();

/*
 * main driver - if argc > 2 we are creating a file from an
 * index, otherwise use the table in h2v.h
 */
main(argc, argv)
	int	argc;
	char	**argv;
{
	FILE	*fp;
	FTAB	table;
	int	i;
	
	if (argc != 2 && argc != 4) {
		fprintf(stderr, "Usage: h2v datafile [indexfile fontfile]\n");
		exit(1);
	}
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "h2v: can't open hersh data file %s\n", argv[1]);
		exit(1);
	}

	readdata(fp);

	if (argc == 4) {
		readindex(argv[2], argv[3], &table);
		writefont(&table);
	} else
		for (i = 0; i < sizeof(fonts) / sizeof(FTAB); i++)
			writefont(&fonts[i]);

	exit(0);
}

/*
 *  readdata
 *
 *  Reads the raw hersh data
 */
void
readdata(fp)
	FILE	*fp;
{
	int	i = 0, j, x, y;
	int	charno, pairs;
	char	buf[MAX_BUF];
	
	while (getcharacter(fp, &charno, &pairs, buf)) {
		hersh[charno - 1].ch = malloc(2 * pairs + 1);
		strcpy(hersh[charno - 1].ch, buf);
		hersh[charno - 1].len = strlen(hersh[charno - 1].ch);
	}

	fclose(fp);
}

/*
 *  readindex
 * 
 *  Read an index file into index tab.
 */
void
readindex(name, fname, tab)
	char	*name, *fname;
	FTAB	*tab;
{
	
	FILE	*fp;
	int	i;

	if ((fp = fopen(name, "r")) == NULL) {
		fprintf(stderr, "h2v: can't open index file\n");
		exit(1);
	}

	tab->name = fname;

	i = 0;
	while (fscanf(fp, "%d %d", &tab->ent[i], &tab->ent[i + 1]) == 2)
		if ((i += 2) >= MAX_ENTS - 2) {
			fprintf(stderr, "h2v: indexfile to big - increase MAX_ENTS\n");
			exit(1);
		}
	
	tab->ent[i] = 0;

	fclose(fp);
}

/*
 * writefont
 *
 *	output a font to file name based on font table tab
 */
void
writefont(tab)
	FTAB	*tab;
{
	int	l ,f, j, j1, x, y, k1, k2, k;
	short	i, nchars, asdecw[3]; 
	short	start, end, nvects, fd;
	char	*p;
	HTAB	*curch;

	fprintf(stderr, "Font name: %s\n", tab->name);

#ifdef PC
	if ((fd = open(tab->name, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0644)) < 0) {
#else
	if ((fd = open(tab->name, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0) {
#endif
		fprintf(stderr, "Can't open output file: %s\n", tab->name);
		exit(1); 
	}

	asdecw[0] = asdecw[2] = -100;
	asdecw[1] = 100;
	nvects = nchars = 0;

	lseek(fd, (long)(5 * sizeof(short)), 0); /*  Leave room for stuff at top */

	for (i = 0; (start = tab->ent[i]) != 0; i += 2) {
		end = tab->ent[i + 1];
#ifdef DEBUG
		fprintf(stderr, "Char: %d to %d\n", start, end);
#endif
		do {
			curch = &hersh[start - 1];
			nchars++;
			if (curch->ch == (char *)NULL) {
				fprintf(stderr, "h2v: character %d not available\n", start);
				exit(1);
			}
			asdecw[2] = MAX(asdecw[2], curch->ch[1] - curch->ch[0]);

#ifdef DEBUG
			fprintf(stderr, "Char: %d length %d\n", start, curch->len);
#endif
			for (p = &curch->ch[2]; *p; p++) {
				x = *p++;
				if (x != ' ') {
					asdecw[0] = MAX(asdecw[0], COORD(*p));
					asdecw[1] = MIN(asdecw[1], COORD(*p));
				}
			}

			nvects += curch->len / 2;
			if (write(fd, &curch->len, sizeof(short)) != sizeof(short)) {
				fprintf(stderr,"h2v: ERROR writing character length to file\n");
				exit(1);
			}
			if (write(fd, curch->ch, (unsigned)curch->len) != curch->len) {
				fprintf(stderr,"h2v: ERROR writing character data to file\n");
				exit(1);
			}
			start++;
		} while (start <= end);
	}

#ifdef DEBUG
	fprintf(stderr,"nchars: %d, nvects: %d\n", nchars, nvects);
	fprintf(stderr,"ascender: %d, decender: %d maxwidth: %d\n",
		asdecw[0], asdecw[1], asdecw[2]);
#endif

	lseek(fd, 0L, 0);
	if (write(fd, &nchars, sizeof(nchars)) != sizeof(nchars)) {
		fprintf(stderr,"Error writing to file\n");
		exit(1);
	}

	if (write(fd, &nvects, sizeof(nvects)) != sizeof(nvects)) {
		fprintf(stderr,"Error writing to file\n");
		exit(1);
	}
	if (write(fd, asdecw, sizeof(asdecw)) != sizeof(asdecw)) {
		fprintf(stderr, "Error writing to file\n");
		exit(1);
	}

	close(fd);
}
