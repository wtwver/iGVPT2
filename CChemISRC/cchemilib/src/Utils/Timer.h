/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
 Copyright (C) 2020 Abdulrahman Allouche (University Lyon 1)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef H_TIMER
#define H_TIMER

#include <sys/time.h>

typedef struct {
  struct timeval begin;
  struct timeval end;
  time_t seconds;
  suseconds_t microseconds;
} TimerType;

#define timer_init(t)				\
  do {						\
    (t).seconds=0;				\
    (t).microseconds=0;				\
  } while(0)
    
#define timer_start(t) gettimeofday(&((t).begin), NULL)

#define timer_stop(t)						\
  do {								\
    gettimeofday(&((t).end), NULL);				\
    (t).seconds += (t).end.tv_sec - (t).begin.tv_sec;		\
    (t).microseconds += (t).end.tv_usec - (t).begin.tv_usec;	\
  } while(0)

#define timer_get(t) (((double)((t).seconds))*1000000. + (double)((t).microseconds))

#endif /* H_TIMER */
