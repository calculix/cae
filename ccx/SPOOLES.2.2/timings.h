#ifndef _TIMINGS_
#define _TIMINGS_
#include <sys/time.h>
static struct timeval  TV ;
static struct timezone TZ ;
#define MARKTIME(t) \
   gettimeofday(&TV, &TZ) ; \
   t = (TV.tv_sec + 0.000001*TV.tv_usec)
#endif

