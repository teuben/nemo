#include <stdio.h>
#include <stdlib.h>
#include "catread.h"
int main()
{
  static char fmt998[] = "              %6d%10.6f%8.4f%8.3f/%s%9.4f\n";
  static double q[7];
  static float smin[100], tempx[100], qx[100], wid[100];
  static int moltag[100];
  static char line[82], buf[82], cfqhi[16], *molnam[100], cfqlow[16];
  FILE *finp, *flist;
  double freq, sval, clight;
  float fqlow, fqhi, dwid;
  int imol, ierr, kmol, iver, nline;

  /* program to list catalog between frequency limits using strength */

  clight = 29979.2458;
  puts("enter input file");
  gets(buf);
  finp = fopen(buf, "r");
  puts("enter list file");
  gets(buf);
  flist = fopen(buf, "w");
  fqlow = 0.; fqhi = 1.e3;
  fgets(line, 82, finp); sscanf(line, "%f%f%f", &fqlow, &fqhi, &dwid);
  fqlow *= clight; fqhi *= clight;
  sprintf(cfqlow,"%13.4f", (double)fqlow);
  sprintf(cfqhi ,"%13.4f", (double)fqhi);
  moltag[0] = 0; kmol = 0;
  while (kmol < 100) {
    smin[kmol] = -3; wid[kmol] = dwid; tempx[kmol] = -0.75;
    if (fgets(buf, 82, finp) == NULL) break;
    if (sscanf(buf,"%d%f%f%f", &moltag[kmol], &smin[kmol], 
	       &wid[kmol], &tempx[kmol]) < 4) break;
    molnam[kmol]= catdir(moltag[kmol], &nline, q, &iver);
    nline = catfrq(moltag[kmol], cfqlow, line);
    if (nline <= 0 || strcmp(line, cfqhi) > 0) continue;
    qx[kmol] = (q[0] - q[1]) * 8.0039228;
    printf(fmt998, moltag[kmol], (double)qx[kmol], (double)wid[kmol],
	   (double)tempx[kmol], &molnam[kmol * 16], (double)smin[kmol]);
    do {
      line[13]=0;
      freq = atof(line) / clight; sval = atof(&line[21]);
      if (sval >= smin[kmol])
	fprintf(flist,"%13.6f%f8.4f%s\n", freq, sval, &line[21]);
      ++nline;
      ierr = catrd(moltag[kmol], nline, line);
    }while (ierr == 0 && strcmp(line, cfqhi) <= 0) ;
    ++kmol;
  }
  fputc('\n', flist);
  for (imol = 0; imol < kmol; ++imol) {
    printf(fmt998, moltag[imol], (double)qx[imol], (double)wid[imol],
	   (double)tempx[imol], &molnam[imol * 16], (double)smin[imol]);
  }
  puts("last species read");
  return 0;
}
