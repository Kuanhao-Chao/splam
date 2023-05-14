#include "GResUsage.h"

#if defined(__APPLE__) && defined(__MACH__)
  #include <AvailabilityMacros.h>
  #include <sys/resource.h>
  #include <mach/mach.h>
  #include <mach/task_info.h>

  #ifndef MAC_OS_X_VERSION_10_12
    #define MAC_OS_X_VERSION_10_12 101200
  #endif
  #if MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_12
    #define G_gettime(s) clock_gettime(CLOCK_MONOTONIC, &s);
  #else
    #include <mach/mach_time.h>
	#define MACHGT_NANO (+1.0E-9)
	#define MACHGT_GIGA UINT64_C(1000000000)
	void mach_gettime( struct timespec* t) {
	  // be more careful in a multithreaded environement
	  static double machgt_timebase = 0.0;
	  static uint64_t machgt_timestart = 0;
	  if (!machgt_timestart) {
		mach_timebase_info_data_t tb;
		tb.numer=0;tb.denom=0;
		mach_timebase_info(&tb);
		machgt_timebase = tb.numer;
		machgt_timebase /= tb.denom;
		machgt_timestart = mach_absolute_time();
	  }
	 ;
	  double diff = (mach_absolute_time() - machgt_timestart) * machgt_timebase;
	  t->tv_sec = diff * MACHGT_NANO;
	  t->tv_nsec = diff - (t->tv_sec * MACHGT_GIGA);
	}
    #define G_gettime(s) mach_gettime(&s)
  #endif
#else
 #ifdef _WIN32
    //Windows implementation:
    #include <psapi.h>
	LARGE_INTEGER
	getFILETIMEoffset() {
	    SYSTEMTIME s;
	    FILETIME f;
	    LARGE_INTEGER t;
	    s.wYear = 1970;  s.wMonth = 1; s.wDay = 1;
	    s.wHour = 0; s.wMinute = 0; s.wSecond = 0;
	    s.wMilliseconds = 0;
	    SystemTimeToFileTime(&s, &f);
	    t.QuadPart = f.dwHighDateTime;
	    t.QuadPart <<= 32;
	    t.QuadPart |= f.dwLowDateTime;
	    return (t);
	}

	static void	usage_to_timeval(FILETIME *ft, struct timeval *tv) {
	    ULARGE_INTEGER time;
	    time.LowPart = ft->dwLowDateTime;
	    time.HighPart = ft->dwHighDateTime;

	    tv->tv_sec = time.QuadPart / 10000000;
	    tv->tv_usec = (time.QuadPart % 10000000) / 10;
	}

    //implementation of getrusage for Windows
	int	getrusage(int who, rusage *usage) {
	    FILETIME creation_time, exit_time, kernel_time, user_time;
	    PROCESS_MEMORY_COUNTERS pmc;

	    memset(usage, 0, sizeof(rusage));

	    if (who == RUSAGE_SELF) {
	        if (!GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time,
	                             &kernel_time, &user_time)) {
	            GMessage("Error: GetProcessTimes() failed!\n");
	            return -1;
	        }

	        if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
	            GMessage("Error: GetProcessMemoryInfo() failed!\n");
	            return -1;
	        }

	        usage_to_timeval(&kernel_time, &usage->ru_stime);
	        usage_to_timeval(&user_time, &usage->ru_utime);
	        usage->ru_majflt = pmc.PageFaultCount;
	        usage->ru_maxrss = pmc.PeakWorkingSetSize / 1024;
	        return 0;
	    } else if (who == RUSAGE_THREAD) {
	        if (!GetThreadTimes(GetCurrentThread(), &creation_time, &exit_time,
	                            &kernel_time, &user_time)) {
	            GMessage("Error: GetThreadTimes() failed!\n");
	            return -1;
	        }
	        usage_to_timeval(&kernel_time, &usage->ru_stime);
	        usage_to_timeval(&user_time, &usage->ru_utime);
	        return 0;
	    } else {
	        return -1;
	    }
	}


	void win_gettime(struct timespec* ts) {
	    LARGE_INTEGER           t;
	    FILETIME            f;
	    double                  microseconds;
	    static LARGE_INTEGER    offset;
	    static double           frequencyToMicroseconds;
	    static int              initialized = 0;
	    static BOOL             usePerformanceCounter = 0;

	    if (!initialized) {
	        LARGE_INTEGER performanceFrequency;
	        initialized = 1;
	        usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
	        if (usePerformanceCounter) {
	            QueryPerformanceCounter(&offset);
	            frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
	        } else {
	            offset = getFILETIMEoffset();
	            frequencyToMicroseconds = 10.;
	        }
	    }
	    if (usePerformanceCounter) QueryPerformanceCounter(&t);
	    else {
	        GetSystemTimeAsFileTime(&f);
	        t.QuadPart = f.dwHighDateTime;
	        t.QuadPart <<= 32;
	        t.QuadPart |= f.dwLowDateTime;
	    }

	    t.QuadPart -= offset.QuadPart;
	    microseconds = (double)t.QuadPart / frequencyToMicroseconds;
	    t.QuadPart = microseconds;
	    ts->tv_sec = t.QuadPart / 1000000;
	    ts->tv_nsec = (t.QuadPart % 1000000)*1000;
	}
   #define G_gettime(s) win_gettime(&s)
 #else  //assume Linux compatible
   #define G_gettime(s) clock_gettime(CLOCK_MONOTONIC, &s)
 #endif
#endif

// Returns the peak (maximum so far) resident set size (physical
// memory use) measured in bytes
size_t getPeakMemUse() {
#if defined(_WIN32)
	// -- Windows
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.PeakWorkingSetSize;
#elif defined(__APPLE__) && defined(__MACH__) && defined(MACH_TASK_BASIC_INFO)
	  struct mach_task_basic_info info;
		mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
		if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
			(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		   return (size_t)info.resident_size_max;
			else
		   return (size_t)0L;		// Can't access?
#else // defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	// asssume BSD, Linux, or OSX
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
	#if defined(__APPLE__)
	  return (size_t)(rusage.ru_maxrss);
	#else //linux returns this in kilobytes
	  return (size_t)(rusage.ru_maxrss*1024L);
	#endif
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes
 */
size_t getCurrentMemUse() {
#if defined(_WIN32)
	// -- Windows
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
#if defined MACH_TASK_BASIC_INFO
	struct mach_task_basic_info info;
#else
  struct task_basic_info info;
  #define MACH_TASK_BASIC_INFO TASK_BASIC_INFO
  #define MACH_TASK_BASIC_INFO_COUNT TASK_BASIC_INFO_COUNT
#endif

	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		return (size_t)0L;		// Can't access?
	return (size_t)info.resident_size;
#else
	//-- assume Linux
	long progsize = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	//if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	if ( fscanf( fp, "%*s%ld", &progsize ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	int page_size=sysconf(_SC_PAGESIZE);
	return ((size_t)progsize * (size_t)page_size);
#endif
}

void printMemUsage(FILE* fout) {
  double rs= getCurrentMemUse();
  rs/=1024;
  fprintf(fout, "Resident Size: %6.1fMB\n", rs);
}

double get_usecTime() {
	struct timespec start_ts;
	G_gettime(start_ts);
	return (((double)start_ts.tv_sec)*1000000.0 + ((double)start_ts.tv_nsec)/1000.0);
}

double GResUsage::start() {
	started=true;
	stopped=false;
	start_mem=getCurrentMemUse();
	double mem=(double)start_mem/1024;
	GMessage("   start_mem=%.2f\n", mem);
	getrusage(RUSAGE_SELF, &start_ru);
	G_gettime(start_ts);
	double tm=start_ts.tv_sec*1000000.0 + start_ts.tv_nsec/1000.0;
	return tm;
}

double GResUsage::stop() {
	if (started!=true)
		GError("Error: calling GResUsage::stop() without starting it first?\n");
	stopped=true;
	G_gettime(stop_ts);
	getrusage(RUSAGE_SELF, &stop_ru);
	double tm=stop_ts.tv_sec*1000000.0 + stop_ts.tv_nsec/1000.0;
	stop_mem=getCurrentMemUse();
	double mem=(double)stop_mem/1024;
	GMessage("   stop_mem=%.2f\n", mem);

	return tm;
}

#define RUSAGE_STOPCHECK

void GResUsage::stopCheck(const char* s) {
	if (!started || !stopped)
      GError("Error: GResUsage::%s() cannot be used before start&stop\n", s);
}

double GResUsage::elapsed() {
	stopCheck("elapsed");
	double st=start_ts.tv_sec*1000000.0 + start_ts.tv_nsec/1000.0;
	double et=stop_ts.tv_sec*1000000.0 + stop_ts.tv_nsec/1000.0;
	return (et-st);
}

double GResUsage::u_elapsed() {
	stopCheck("u_elapsed");
	double st=start_ru.ru_utime.tv_sec*1000000.0 + start_ru.ru_utime.tv_usec;
	double et=stop_ru.ru_utime.tv_sec*1000000.0 + stop_ru.ru_utime.tv_usec;
	return (et-st);
}

double GResUsage::s_elapsed() {
	stopCheck("s_elapsed");
	double st=start_ru.ru_stime.tv_sec*1000000.0 + start_ru.ru_stime.tv_usec;
	double et=stop_ru.ru_stime.tv_sec*1000000.0 + stop_ru.ru_stime.tv_usec;
	return (et-st);
}

double GResUsage::memoryUsed() { //in kilobytes
	stopCheck("memoryUsed");
	return (((double)stop_mem-(double)start_mem)/1024.0);
}

