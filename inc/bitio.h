/*
 *   manipulate bits  
 *   in linux see also /usr/include/sys/param.h:   setbit, clrbit, isset, isclr
 *
 */

#define SETBIT(word, n) ((word) |=  (1 << (n)))   // set bit n of word
#define GETBIT(word, n) ((word) &   (1 << (n)))   // get bit n of word
#define CLRBIT(word, n) ((word) &= ~(1 << (n)))   // clear bit n of word

