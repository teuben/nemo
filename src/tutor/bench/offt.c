/*
 *  compile this with      
	gcc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
    and your off_t will be 64 bit

    Reminders:

	off_t lseek(int fildes, off_t offset, int whence);
	ssize_t write(int fd, const void *buf, size_t count);

 *
 */

#include <sys/types.h>
#include <unistd.h>

main() {			// normal       -D....        true
				// ia-32 mode	ia-64 mode	ia-64
    short  i0 = 1;              // 2            2               2
    off_t  i1 = 1;		// 4 		8		8
    size_t i2 = 1;		// 4		4		8
    long   i3 = 1;		// 4		4		8
    long long   i4 = 1;		// 8		8		8
    ssize_t     i5  = 1;	// 4		4		8
    
    int l0 = sizeof(i0);
    int l1 = sizeof(i1);
    int l2 = sizeof(i2);
    int l3 = sizeof(i3);    
    int l4 = sizeof(i4);
    int l5 = sizeof(i5);

    printf("%d %d %d %d %d %d\n",l0,l1,l2,l3,l4,l5);
}
