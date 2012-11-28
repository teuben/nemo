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
 * Subroutine:     SETFBLANK
 *
 * Purpose:        Subroutine to set a data value to the universal BLANK.
 *
 * Files:          SETFBLANK.C, SETFBLANK.DC2
 *
 * Language:       ANSI C
 *
 * Author:         K.G. Begeman
 *
 * Address:        KGB@HGRRUG5
 *
 * Fortran call:   CALL SETFBLANK( DATA )
 *
 *                 DATA   REAL    Data will contain the value for BLANK
 *                                on return.
 *
 * C call:         void setfblank_( data )
 *
 *                 real  *data    data will contain the value for BLANK
 *                                on return.
 *
 * Note:           The logical function FBLANK / fblank_() returns whether
 *                 a data value is BLANK or not. See FBLANK.DC2.
 *
 * Updates:        24/oct/88 KGB Document created.
 *                 10/sep/08 PJT g++ friendly
 *
 */

#include <stdio.h>
#include "gipsyc.h"

void setfblank_( void *value )
{
   integer *ivalue = (integer *) value;
   *ivalue = BLANKVAL;
}

