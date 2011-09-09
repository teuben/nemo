/*
 * pure C program showing how you can use popen() to perform
 * a NEMO task without the need to link NEMO code.
 * Not very efficient, but it works.
 */

#include <stdio.h>

#define N 1000

void main() {
    float x[N], y[N], z[N], pot[N];
    int i, n;
    FILE *fp;

    /* 
       before you run this program, create p100.tab as follows 
           mkplummer - 100 | snapprint - x,y,z > p100.tab  
    */
    n = 100;
    fp = fopen("p100.tab","r");
    for (i=0; i<n; i++)  fscanf(fp,"%g %g %g\n",&x[i],&y[i],&z[i]);
    fclose(fp);


    /* silly, but write them out again */

    fp = fopen("tmpxyz.tab","w");
    for (i=0; i<n; i++)   fprintf(fp,"%g %g %g\n",x[i],y[i],z[i]);
    fclose(fp);

    /* use popen(2) to read the returned potentials */

    fp = popen("xyz2pot tmpxyz.tab","r");
    for (i=0; i<n; i++)  {
      fscanf(fp,"%g\n",&pot[i]);
      printf("%d : %f %f %f   %f\n",i+1,x[i],y[i],z[i],pot[i]); 
    }
    pclose(fp);

}


