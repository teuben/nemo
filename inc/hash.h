/*


 Copyright (C) 1990 Texas Instruments Incorporated.

 Permission is granted to any individual or institution to use, copy, modify,
 and distribute this software, provided that this complete copyright and
 permission notice is maintained, intact, in all copies and supporting
 documentation.

 Texas Instruments Incorporated provides this software "as is" without
 express or implied warranty.


 *  
 * I/O support for internal defmacros
 *   Defmacros which execute in their own process use stdin and stdout.
 *   Defmacros which execute as cpp functions use these routines instead.
 *
 * Edit history
 * Created: LGO 30-Mar-89 -- Initial design and implementation.
 *          PJT 28-Apr-92 -- renamed this defmacio.h to hash.h
 *                           and stripped most of the bottem half
 *                           Boolean is now NEMO's bool
 *	    pjt 21-sep-93 -- ansi
 *
 */

#include <stdinc.h>
#include <stdlib.h>	/* for getmem() */


/************************************************************/

/*  I moved these definitions in here from hash.c because the SAS C
 *  compiler (MVS) complains about not knowing about type Hash_Table
 *  in declaration    extern struct Hash_Table* init_Hash_Table(); 
 *  This looks like a bug their compiler ... DKM
 */

#define BUCKET_SIZE 8

struct hash_pair {
  char* key;
  void* value;
};

struct bucket {
  struct hash_pair data[BUCKET_SIZE];
};

struct Hash_Table {
  struct bucket* table;
  unsigned char* items_in_buckets;
  int entry_count;			
  int bucket_index;			
  int current_bucket;
  int current_index;
};

/**************************************************************/


extern struct Hash_Table* init_Hash_Table(void);
extern void*              get_hash(struct Hash_Table *, string);
extern bool               put_hash(struct Hash_Table *, string, void *);
extern void*              next_hash(struct Hash_Table *, bool);

/* -- still some functions to complete here */
/* get_entry_count  */

