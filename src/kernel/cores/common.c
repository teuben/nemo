/*
 * common.c:    routines to access and manage common blocks.
 *		these are similar to malloc'd buffers, but do not
 *		get free'd as this can be very wastful.
 *
 *	8-jan-99    Created for random I/O buffers; id not used yet.
 *     13-jan-98    free old buffer when shrinking buffersize
 *     28-nov-00    fixed return type of get_common
 *     20-jun-01    gcc3
 */
#include <stdinc.h>

static byte *random_buffer = NULL;
static int   random_buffer_size = 10240;     /* # bytes for buffer */
static int   random_buffer_lock = 0;

void set_common(int id, int byte_size)
{
    if (byte_size < 1) {
        warning("set_common: bad size=%d",byte_size);
        return;
    }
    if (byte_size < random_buffer_size) {
        free(random_buffer);
        random_buffer = NULL;
    }
    if (random_buffer == NULL) {
        random_buffer = (byte *) allocate(byte_size);
        random_buffer_size = byte_size;
    }
    if (byte_size > random_buffer_size) {
        random_buffer = (byte *) reallocate(random_buffer,byte_size);
        random_buffer_size = byte_size;
    }
}

int get_common(int id, int elt_size, int bucket_size)
{
    int n = random_buffer_size/(elt_size * bucket_size);
    if (n < 1) 
        error("get_common(%d): bad (buffer_size=%d)/(elt=%d * bucket=%d) \n",
                id,random_buffer_size,elt_size,bucket_size);
    return 0;
}

byte *open_common(int id)
{
    if (random_buffer_lock) error("open_common: already in use");
    if (random_buffer == NULL)
        random_buffer = (byte *) allocate(random_buffer_size);

    random_buffer_lock = 1;
    return random_buffer;
}

void close_common(int id)
{
    random_buffer_lock = 0;
}

