/*
 *
 *  cf: reblock a fits file that was written with fortran unformatted
 *      2880 blocked I/O
 *  TODO: the functionality of this program should  be added to
 *	'unfio' with an out= option
 *
 *	june 1995	first version ?			pjt
 *	june 1996	put in NEMO, but simple main() user interface	pjt
 */


#include <stdio.h>

main(int ac,char *av[])
{
    FILE *fi, *fo;
    char buf[2888];
    int n;

    if (ac!=3) usage(av[0]);

    fi = fopen(av[1],"r");
    if (fi==NULL) {
        printf("Cannot open File %s for input\n",av[1]);
        exit(-1);
    }
    fo = fopen(av[2],"w");
    if (fo==NULL) {
        printf("Cannot open File %s for input\n",av[2]);
        exit(-1);
    }

    n=0;
    while ( fread(buf,1,2888,fi)==2888 ) {
        n++;
        if (fwrite(&buf[4],1,2880,fo) != 2880) {
            perror("Error during write");
            exit(-1);
        }
    }
    printf("%s: converted %d blocks\n",av[1],n);
    fclose(fi);
    fclose(fo);
    exit(0);
}

usage(char *name)
{
    printf("Usage: %s infile outfile\n",name);
    printf("Reblocks a 2888-block based FITS file to 2880-blocked\n");
    printf("under the assumption that the first and last 'word'\n");
    printf("of each 2888-byte block is the block count (2880)\n");
    exit(0);
}
