/*
 * try and decipher the convention of fortran symbols
 * as they become visible for outsiders like C and C++
 *
 * Note: g77/f2c can be controlled with the
 *	-fno-underscoring
 * or
 *	-fno-second-underscore
 * compile flag. Those situations are hard to control/check
 * here?
 */

#if !defined(VMS) && !defined(aix)

# if defined(unicos)
#define f2c_ucase
# else
#define f2c_trail
# endif

#endif

