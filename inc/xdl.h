/* for NEMO adapted : see src/kernel/loadobj/xdl */

#ifndef __XDL_H__
#define __XDL_H__

int xdl_init( /* char *__executable_name */ );

int xdl_link( /* char **__null_terminated_object_file_name_array */ );

void *xdl_symbol( /* char *__function_name */ );

void xdl_finish( /* void */);

#endif  __XDL_H__

