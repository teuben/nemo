/*
 *	although never called, it forces the linker to load
 *	most (all) filestruct routines for a dynamic object
 *	loader.
 */

#include <filestruct.h>
 
static void fs_linker(dstr)
stream dstr;
{   
   int    i;

   error("fs_linker cannot be called");
   (void) allocate(1);
   (void) stropen("Dummy","r");
   (void) get_tag_ok(dstr,"Dummy");
   (void) get_string(dstr,"Dummy");
   (void) get_set(dstr,"Dummy");
   (void) get_data(dstr,"Dummy","Dummy",&i,0);
   (void) get_tes(dstr,"Dummy");
   (void) strclose(dstr);
}

