// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/inc/timer.h
///
/// \author  Song Ho Ahn (song.ahn@gmail.com), Walter Dehnen (wd11@le.ac.uk)
///                                                                             
/// \date    2003-2010
/// \version 13-jan-2003  SHA  created
/// \version 13-jan-2006  SHA  updated
/// \version 05-jun-2010  WD   adapted
/// \version 06-oct-2010  WD   improved
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
#  define WDutils_included_cstdlib
#  include <cstdlib>
#endif


#define TIMER_SUPPORT_FOR_DEPRACATED

namespace WDutils {
  /// High Resolution timer.
  /// this timer measures elapsed wall-clock time with micro-second accuracy
  class Timer {
  private:
#ifdef WIN32
    LARGE_INTEGER frequency;               ///< ticks per second
    LARGE_INTEGER oldCount;                ///< counter at last call to take()
#else
    timeval oldCount;                      ///< counter at last call to take()
#endif
#ifdef TIMER_SUPPORT_FOR_DEPRACATED
    double taken;                          ///< time taken in seconds
#endif
    /// read clock, compute elapsed time in @a t, return clock reading
#ifdef WIN32
    LARGE_INTEGER getTime(double&t) const
    {
      LARGE_INTEGER newCount;
      QueryPerformanceCounter(&newCount);
      t = double(newCount.QuadPart - oldCount.QuadPart)
	/ double(frequency.QuadPart);
      return newCount;
    }
#else
    timeval getTime(double&t) const
    {
      timeval newCount;
      gettimeofday(&newCount, NULL);
      t = double(newCount.tv_sec -oldCount.tv_sec )
	+ double(newCount.tv_usec-oldCount.tv_usec) * 1.e-6;
      return newCount;
    }
#endif
  public:
    /// start timer: take clock reading and reset counter
    void start()
    { 
#ifdef WIN32
      QueryPerformanceCounter(&oldCount);
#else
      gettimeofday(&oldCount, NULL);
#endif
    }
    /// ctor: start timer
    Timer()
#ifdef TIMER_SUPPORT_FOR_DEPRACATED
      : taken(0.0)
#endif
    {
#ifdef WIN32
      QueryPerformanceFrequency(&frequency);
#endif
      start();
    }
    /// stop timer and return elapsed time: take wallclock-time, compute
    /// elapsed time, and reset start counter.
    /// \return wallclock-time (in sec) since last call to start() or stop()
    double stop()
    {
#ifndef TIMER_SUPPORT_FOR_DEPRACATED
      double  taken;
#endif
      oldCount = getTime(taken);
      return taken;
    }
    /// take time since last call to start() or stop()
    /// \note Does not reset the start counter.
    /// \return wallclock-time (in sec) since last call to start() or stop()
    double take() const
    {
      double  t;
      getTime(t);
      return t;
    }
#ifdef TIMER_SUPPORT_FOR_DEPRACATED
    /// elapsed time in seconds with micro-second resolution.
    /// \note deprecated and provided for backwards compatibility only
    double getElapsedTime() const
    { return taken; }
#endif
  };// class Timer
} // namespace WDutils
// /////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_timer_h
