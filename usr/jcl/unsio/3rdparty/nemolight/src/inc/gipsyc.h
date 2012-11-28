/*
 *                            COPYRIGHT (c) 1988
 *                      Kapteyn Astronomical Institute
 *       University of Groningen  -  9700 AV Groningen, The Netherlands
 *
 *      This software is furnished under a license for use  only  on  a
 *      single  computer  system  and  may  be  copied  only  with  the
 *      inclusion of the above copyright notice.
 *
 *      This software, or any other copies thereof, may not be provided
 *      or  otherwise made available to any other person except for use
 *      on such system and to the  one  who  agrees  to  these  license
 *      terms.
 *
 *      Title to and ownership of  the  software  shall  at  all  times
 *      remain the property of the University of Groningen.
 *
 *      The information in this document is subject to  change  without
 *      notice  and  should  not  be  construed  as a commitment by the
 *      Kapteyn Astronomical Institute.
 *
 *      The Kapteyn Astronomical Institute  assumes  no  responsibility
 *      for the use or reliability of this software.
 *
 *      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *
 * Include File:   GIPSYC.H
 *
 * Purpose:        This include file should be included in all C sources
 *                 of GIPSY level one, two and three routines.
 *
 * Author:         K.G. Begeman
 *
 * Contents:       The include file defines the ANSI-F77 types:
 *                 character (for VMS only)
 *                 integer
 *                 logical
 *                 real
 *                 double
 *                 complex
 *                 true and false are defined
 *                 the common blocks for logical function undefined()
 *                 the system value for BLANK
 *
 * Updates:        21/sep/88 KGB creation date
 *
 */

#if defined(VMS)            /* get predefined VMS character descriptor struct */
#include <descrip.h>
#define character struct dsc$descriptor_s
#endif

#define integer   long                     /* equivalent with F77 INTEGER     */
#define logical   long                     /*   .....     ..   .  LOGICAL     */
#define real      float                    /*   .....     ..   .  REAL        */
#define complex   struct { float dr, di; } /*   .....     ..   .  COMPLEX     */
#define true      1                        /* this is true (also in F77)      */
#define false     0                        /* this is false (also in F77)     */

#if defined(VMS)
#define UNDF  undef                        /* link with routine UNDEFINED     */
#define CUNDF cundef
#else
#define UNDF  undef_                       /* link with routine UNDEFINED     */
#define CUNDF cundef_
#endif

extern char UNDF, CUNDF;                   /* COMMON BLOCK for defined()      */

#if defined(VMS)
#define BLANKVAL 0x00008000                /* Universal value for BLANK       */
#else
#define BLANKVAL 0x7fffffff                /* Universal value for BLANK       */
#endif

