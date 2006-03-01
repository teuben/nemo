/*
 * unformatted fortran I/O using C
 */

extern int unfsize(int);
extern int unfswap(bool);
extern int unfscan(stream fp);
extern int unfread(stream fp, char *buf, int bufsize);
/* TODO::                                ^^^
                                         int should become size_t
 */
