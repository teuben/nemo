/* MATDEF.H	
 *		defines how n-dimensional images are stored in image(5NEMO) 
 *              programs. Mostly 2 and 3 dim images. 
 *              See also mdarray.h for arrays of higher dimension.
 *		This include file is automatically included by image.h
 *		and perhaps by a few others, e.g. contour.c
 *
 *              The use of Karma-style intelligent array indexing is
 *              also defined here.
 *
 *	#define CDEF			c-definition        = row major
 *	#define FORDEF			fortran-definition  = column major
 *
 *	As of January 2013 the default is FORDEF
 * 
 *   C:  data[row][col]     row major,     row after row
 *   F:  data(row,col)      column major,  column after column
 *
 *
 *      #define USE_IARRAY
 */


//#define CDEF

#if !defined(FORDEF) && !defined(CDEF)
# define FORDEF
#endif

// #define USE_IARRAY
