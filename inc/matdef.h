/* MATDEF.H	
 *		defines how matrices are stored in image(5NEMO) programs
 *		This include file is automatically included by image.h
 *		and perhaps by a few others, e.g. contour.c
 *
 *	#define CDEF			c-definition
 *	#define FORDEF			fortran-definition
 *
 *	The default is CDEF
 *
 */
#define CDEF 
#if !defined(FORDEF) && !defined(CDEF)
# define FORDEF
#endif

