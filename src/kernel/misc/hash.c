/*


 Copyright (C) 1990 Texas Instruments Incorporated.

 Permission is granted to any individual or institution to use, copy, modify,
 and distribute this software, provided that this complete copyright and
 permission notice is maintained, intact, in all copies and supporting
 documentation.

 Texas Instruments Incorporated provides this software "as is" without
 express or implied warranty.


 *
 * Edit history
 * Created: LGO 30-Mar-89 -- Initial design and implementation.
 *          PJT 28-apr-92 -- Testbed for standalone; away from 'cpp'
 *			     EOS->NULL
 *	    pjt 22-feb-94 -- ansi headers - ANSI warning (see *CAST* comment)
 *	    pjt 20-jun-01    gcc 3
 * Hash functions
 */

#include <stdinc.h>	/* just for getmem() */
#include <strlib.h>
#include <hash.h>

static int hash_primes[] = {3, 7, 13, 19, 29, 41, 53, 67, 83, 97, 113, 137,
			    163, 191, 223, 263, 307, 349, 401, 461, 521, 587,
			    653, 719, 773, 839, 911, 983, 1049, 1123, 1201,
			    1279, 1367, 1459, 1549, 1657, 1759, 1861, 1973,
			    2081, 2179, 2281, 2383, 2503, 2617, 2729, 2843,
			    2963, 3089, 3203, 3323, 3449, 3571, 3697, 3833,
			    3967, 4099, 4241, 4391, 4549, 4703, 4861, 5011,
			    5171, 5333, 5483, 5669, 5839, 6029, 6197, 6361,
  			    6547, 6761, 6961, 7177, 7393, 7517, 7727, 7951,
			    8101, 8209, 16411, 32771, 65537, 131301, 262147,
			    524287};

void*
next_hash(
   struct Hash_Table* h,
   bool firstp)
{
  if(firstp) {
    h->current_bucket = 0;
    h->current_index = 0;
  }
  while (h->current_bucket < hash_primes[h->bucket_index]) {
    if (h->current_index < h->items_in_buckets[h->current_bucket])     /*CAST*/
      return h->table[h->current_bucket].data[h->current_index++].value;
    else {
      h->current_bucket += 1;
      h->current_index = 0;
    }
  }
  return NULL;
}

/* sxhash -- Hash function for char*        (seems local to me; pjt)
 * Input:	Character string
 * Output:	unsigned long hash value
 */
static unsigned long sxhash(char *s)
{
  register unsigned long hash = *s++;
  if(*s != 0) {
    hash = (hash << 7) ^ *s++;
    if (*s != 0) {
      hash = (hash << 7) ^ *s++;
      if (*s != 0) {
	hash = (hash << 7) ^ *s++;
	while (*s != 0) { /* rotate hash left 7 bits */
	  int rest = hash >> 25;
	  hash = ((hash << 7) | rest) ^ *s++;
	}
      }
    }
  }
  return hash;
}

struct Hash_Table*
init_Hash_Table() {
  struct Hash_Table* h = (struct Hash_Table*) getmem(sizeof(struct Hash_Table));
  int prime, i;
  h->entry_count = 0;	
  h->bucket_index = 0;	
  prime  = hash_primes[h->bucket_index];	
  h->items_in_buckets = (unsigned char*) getmem (prime);	
  for (i = 0; i < prime; i++)		
    h->items_in_buckets[i] = 0;
  h->table = (struct bucket*) getmem(sizeof(struct bucket)*prime);
  return h;
}

static void
resize (struct Hash_Table* h, int n)
{
  unsigned char* temp1;
  struct bucket* temp2;			
  long old_prime, i, j, new_prime, hash;
  old_prime = hash_primes[h->bucket_index];

  while (hash_primes[h->bucket_index]*BUCKET_SIZE < n)
    h->bucket_index++;

 retry:
  new_prime = hash_primes[h->bucket_index];
  temp1 = (unsigned char*) getmem (new_prime);		
  for (i = 0; i < new_prime; i++)		
    temp1[i] = 0;				
  temp2 = (struct bucket*) getmem(sizeof (struct bucket) * new_prime); 
  for (i = 0; i < old_prime; i++) {	
    for (j = 0; j < h->items_in_buckets[i]; j++) { 	/*CAST*/
      char* key = h->table[i].data[j].key;
      hash = sxhash(key) % new_prime;
      if (temp1[hash] >= BUCKET_SIZE) {	  /* Bucket Overflow */
	free (temp1);
	free (temp2);
	h->bucket_index++;
	goto retry;
      }
      temp2[hash].data[temp1[hash]].key = key;	
      temp2[hash].data[temp1[hash]].value = h->table[i].data[j].value;	
      ++temp1[hash];
    }
  }
  free (h->items_in_buckets);
  free (h->table);			
  h->items_in_buckets = temp1;		
  h->table = temp2;			
}

int
get_entry_count (struct Hash_Table* h)
{
  return (h->entry_count);
}

bool
put_hash (
    struct Hash_Table* h,
    char* key,
    void* value)
{
  long prime,hash;
  int next, i;
 retry:
  prime = hash_primes[h->bucket_index];
  hash = sxhash(key) % prime;
  for (i = 0; i < h->items_in_buckets[hash]; i++)	/*CAST*/
    if (!strcmp(key,h->table[hash].data[i].key))
      return FALSE;
  if (h->items_in_buckets[hash] >= BUCKET_SIZE) {
    resize (h, hash_primes[h->bucket_index+1]*BUCKET_SIZE);
    goto retry;
  }
  next = h->items_in_buckets[hash]++;
  h->table[hash].data[next].key = key;
  h->table[hash].data[next].value = value;
  h->entry_count++;			
  return TRUE;
}

void*
get_hash (
     struct Hash_Table* h,
     char* key)
{
  long prime, hash;
  int i, len;
  prime = hash_primes[h->bucket_index];
  hash = sxhash(key) % prime;
  len = h->items_in_buckets[hash];
  for (i = 0; i < len; i++) {
    if (!strcmp(key,h->table[hash].data[i].key))
      return (h->table[hash].data[i].value);
  }
  return NULL;
}



#if defined(TESTBED)

#include <getparam.h>

string defv[] = {
    "key=???\n      Input table with keys",
    "val=???\n      Input table with values",
    "s=\n           Search keys to test for",
    "show=f\n       Show all values in hashed order?",
    "VERSION=1.0\n  28-apr-92 PJT",
    NULL,
};

string usage="lookup table using hash table";

#define LENKEY 64
#define LENVAL 256

nemo_main()
{
    struct Hash_Table  *h;
    int i, nmax, nhash=0;
    char key[LENKEY], val[LENVAL], *scopy(), *cp;
    string input, *keys, *vals, *burststring(), *sp, *s;
    stream inkey, inval;
    bool first=TRUE;

    input = getparam("key");
    inkey = stropen(input,"r");
    nmax = file_lines(input,10000);
    if (nmax<=0) error("nmax=%d; no good support for pipes",nmax);

    input = getparam("val");
    inval = stropen(input,"r");


    keys = (string *) allocate(nmax*sizeof(string));
    vals = (string *) allocate(nmax*sizeof(string));

    h = init_Hash_Table();

    for (i=0; i<nmax; i++) {
        if(fgets(key,LENKEY,inkey) == NULL)
            error("Too many lines in key file");
        if(fgets(val,LENVAL,inval) == NULL)
            error("Too many lines in key file");
        key[strlen(key)-1]=NULL;        /* strip NULL */
        val[strlen(val)-1]=NULL;
        keys[i] = scopy(key);
        vals[i] = scopy(val);
        if (!put_hash(h,keys[i],vals[i]))
            warning("Record %d: cannot store hash value for %s",i+1,keys[i]);
        else {
            nhash++;
            dprintf(2,"Entering %s=%s\n",keys[i],vals[i]);
        }
    }
    dprintf(1,"Entered %d keys\n",nhash);

    if (hasvalue("s")) {
        sp = burststring(getparam("s")," ,");
        for (s=sp; *s; s++)
            printf("%s -> %s\n",*s,get_hash(h,*s));
    }

    if(getbparam("show")) {
        while( (cp=next_hash(h,first)) ){
            printf("%s\n",cp);
            first = FALSE;
        }
    }        
}
#endif
