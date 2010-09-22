// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/inc/timer.h
///
/// \author  Song Ho Ahn (song.ahn@gmail.com), Walter Dehnen (wd11@le.ac.uk)
///                                                                             
/// \date    2003-2010
/// \version 13-jan-2003 created SHA
/// \version 13-jan-2006 updated SHA
/// \version 05-jun-2010 adapted WD
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2003-2010 Song Ho Ahn, Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_timer_h
#define WDutils_included_timer_h

extern "C" {
#ifdef WIN32   // Windows system specific
#  include <windows.h>
#else          // Unix based system specific
#  include <sys/time.h>
#endif
}

#ifndef WDutils_included_cstdlib
#define WDutils_included_cstdlib
#  include <cstdlib>
#endif

namespace WDutils {
  /// High Resolution timer.
  /// this timer measures elapsed time with 1 micro-second accuracy
  class Timer {
  private:
#ifdef WIN32
    LARGE_INTEGER frequency;               ///< ticks per second
    LARGE_INTEGER startCount;              ///< value of counter at start
    LARGE_INTEGER endCount;                ///< value of counter at end
#else
    timeval startCount;                    ///< value of counter at start
    timeval endCount;                      ///< value of counter at end
#endif
    bool    stopped;                       ///< stop flag 
  public:
    /// ctor
    Timer()
      : stopped(true)  {}
    /// start timer.
    /// startCount will be set here
    void start()
    {
      stopped = false; // reset stop flag
#ifdef WIN32
      QueryPerformanceCounter(&startCount);
#else
      gettimeofday(&startCount, NULL);
#endif
    }
    /// stop the timer.
    /// endCount will be set at this point.
    void stop()
    {
      stopped = true;  // set timer stopped flag
#ifdef WIN32
      QueryPerformanceCounter(&endCount);
#else
      gettimeofday(&endCount, NULL);
#endif
    }
    /// compute elapsed time with micro-second resolution.
    double getElapsedTime()
    {
      if(!stopped) stop();
#ifdef WIN32
      return 
	double(endCount.QuadPart - startCount.QuadPart) /
	double(frequency.QuadPart);
#else
      return
	double(endCount.tv_sec-startCount.tv_sec) +
	double(endCount.tv_usec-startCount.tv_usec) * 1.e-6;
#endif
    }
  };// class Timer
} // namespace WDutils
// /////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_timer_h
