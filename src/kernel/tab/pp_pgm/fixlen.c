#include <stdio.h>
int main(int argc, char **argv)
{
  FILE *finp, *fout;
  int ich, nch, nline;
  char buf[82];
  if (argc < 3) {
    printf("usage fixlen inpfile outfile");
    exit(1); 
  }
  finp = fopen(argv[1],"r");
  if (finp == NULL) {
    printf("bad input file: %s\n", argv[1]);
    exit(1);
  }
  fout = fopen(argv[2],"w");
  if (fout == NULL) {
    printf("bad output file: %s\n", argv[2]);
    exit(1);
  }
  nch = nline = 0;
  buf[80] = '\n'; buf[81] = 0;
  while ((ich = getc(finp)) > 0) {
    if (ich == '\n') {
      for (; nch < 80; ++nch) buf[nch] = ' ';
      fputs (buf, fout);
      nch = 0; ++nline;
    } else if (nch < 80 && ich >= 32) {
      buf[nch] = ich; ++nch;
    }
  }
  fclose(finp); fclose(fout);
  printf("%d lines\n", nline);
  return 0;
}
