#ifndef _GRESUSAGE_
#define _GRESUSAGE_
#include "GBase.h"
#if defined _WIN32 && ! defined __CYGWIN__
  #define	RUSAGE_SELF	0		/* calling process */
  #define	RUSAGE_CHILDREN	-1		/* terminated child processes */
  #define	RUSAGE_THREAD	1

  struct rusage {
	struct timeval ru_utime;	/* user time used */
	struct timeval ru_stime;	/* system time used */
	long ru_maxrss;
	long ru_majflt;
  };
#else
 #include <sys/resource.h>
#endif
#include <time.h>

// report the memory usage of the current process in bytes
size_t getCurrentMemUse(); //current memory usage of the program (RSS)
size_t getPeakMemUse(); //maximum memory usage (RSS) for the program until now

void printMemUsage(FILE* fout=stderr); //in kilobytes

double get_usecTime();

class GResUsage {
  protected:
	bool started;
	bool stopped;
	size_t start_mem;
	size_t stop_mem;
	struct rusage start_ru;
	struct rusage stop_ru;
	struct timespec start_ts;
	struct timespec stop_ts;
	void stopCheck(const char* s);
  public:
	GResUsage(bool do_start=false):started(false),
	   stopped(false), start_mem(0), stop_mem(0) {
          if (do_start) start();
	}

	double start(); //returns microseconds time using clock_gettime(CLOCK_MONOTONIC)
	double stop(); //stop the stopwatch, returns the current time in microseconds
	double elapsed(); //microseconds elapsed between start and stop (wallclock time)
	double u_elapsed(); //microseconds of user time elapsed
	double s_elapsed(); //microseconds of system time elapsed
	double memoryUsed(); //memory increase between start and stop in KB (can be negative)
};

#endif
