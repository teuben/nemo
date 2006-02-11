/* 
 * in which order arguments are evaluated in C
 * i.e. never depend on it!
 *
 * See e.g. http://www.linuxjournal.com/article/8850
 */

#include <stdio.h>

int *inc(int *x)
{
  printf("%d\n", *x);
  (*x)++;
  return x;
}

void printboth(int x, int y)
{
  printf("x: %d, y: %d\n", x, y);
}

int main(int argc, char *argv[])
{
  int x = 0;
  printboth(*inc(&x) + 10, 
	    *inc(&x) + 100);
  printf("Linux i386:   x: 12, y: 101\n");
  printf("Mac OS 10.4:  x: 11, y: 102\n");
  printf("Solaris 8:    x: 11, y: 102\n");
  return 0;
}
