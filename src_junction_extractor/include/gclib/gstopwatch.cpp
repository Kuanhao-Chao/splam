#include "gstopwatch.h"

#ifdef _WIN32
double GStopWatch::LIToSecs( LARGE_INTEGER & L) {
  return ((double)L.QuadPart /(double)frequency.QuadPart);
}

GStopWatch::GStopWatch(){
  timer.start.QuadPart=0;
  timer.stop.QuadPart=0;  
  QueryPerformanceFrequency( &frequency );
}

void GStopWatch::startTimer( ) {
    QueryPerformanceCounter(&timer.start);
}

void GStopWatch::stopTimer( ) {
    QueryPerformanceCounter(&timer.stop);
}


double GStopWatch::getElapsedTime() {
  LARGE_INTEGER time;
  time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
    return LIToSecs( time) ;
}
#else
//Linux code:
void GStopWatch::startTimer( ) {
  gettimeofday(&(timer.start),NULL);
}

void GStopWatch::stopTimer( ) {
  gettimeofday(&(timer.stop),NULL);
}

double GStopWatch::getElapsedTime() {  
  timeval res;
  timersub(&(timer.stop),&(timer.start),&res);
  return res.tv_sec + res.tv_usec/1000000.0; // 10^6 uSec per second
}

#endif
