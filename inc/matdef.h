/* MATDEF.H	
 *		defines how matrices are stored in image(5NEMO) programs
 *		This include file is automatically included by image.h
 *		and perhaps by a few others, e.g. contour.c
 *
 *              The use of Karma-style intelligent array indexing is
 *              also defined here.
 *
 *	#define CDEF			c-definition        = row major
 *	#define FORDEF			fortran-definition  = column major
 *
 *	The default is CDEF
 * 
 *   C:  data[row][col]     row major,     row after row
 *   F:  data(row,col)      column major,  column after column
 *
 *
 *      #define USE_IARRAY
 */


#define CDEF 

#if !defined(FORDEF) && !defined(CDEF)
# define FORDEF
#endif

// #define USE_IARRAY
