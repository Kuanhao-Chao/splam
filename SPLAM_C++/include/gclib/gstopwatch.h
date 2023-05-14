#ifndef __GSTOPWATCH_H
#define __GSTOPWATCH_H
#include "GBase.h"

#ifdef _WIN32
typedef struct {
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} stopWatch;

class GStopWatch {

private:
  stopWatch timer;
  LARGE_INTEGER frequency;
  double LIToSecs( LARGE_INTEGER & L);
public:
  GStopWatch();
  void startTimer( );
  void stopTimer( );
  double getElapsedTime();
};

#else
#include <sys/time.h>

typedef struct {
  timeval start;
  timeval stop;
} stopWatch;

class GStopWatch {

private:
  stopWatch timer;
public:
  GStopWatch() {};
  void startTimer( );
  void stopTimer( );
  double getElapsedTime();
};

#endif

#endif
