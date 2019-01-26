#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

/*
 * header only library to take snapshots of a process status
 * like memory, cpu time, (etc.. in the future)
 */

/*
 /proc/self/statm for current resident set size

 On Linux, the /proc pseudo-file system contains a directory for each running or zombie process. 
 The /proc/[pid]/stat, /proc/[pid]/statm, and /proc/[pid]/status pseudo-files for the process 
 with id [pid] all contain a process's current resident set size IN PAGES, among other things. 
 But the /proc/[pid]/statm pseudo-file is the easiest to read since it contains a single line of text 
 with white-space delimited values:

    total program size
    resident set size
    shared pages
    text (code) size
    library size
    data size (heap + stack)
    dirty pages

 To get the field for the current process, open /proc/self/statm and parse the second integer value. Multiply the field by the page size from sysconf()

 http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
 
 http://elinux.org/Runtime_Memory_Measurement
 */
 
#ifndef __SELFIE_H__
#define __SELFIE_H__

typedef struct statm {
	long total_program_size;
    long resident_set_size;
    long shared_pages;
    long code_size;
    long library_size;
    long data_size;
    long dirty_pages;
} statm;

int getstatm(statm * self) {
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return -2;		// Can't open? */
	if ( fscanf( fp, "%ld%ld%ld%ld%ld%ld%ld", &(self->total_program_size), &(self->resident_set_size), &(self->shared_pages), &(self->code_size), &(self->library_size), &(self->data_size), &(self->dirty_pages) ) != 7 )
	{
		fclose( fp );
		return -4;		// Can't read? */
	}
	fclose( fp );
	return 0;
}

int getstatm_qd(statm * self) { //quick and dirty
	FILE* fp = fopen( "/proc/self/statm", "r" );
	fscanf( fp, "%ld%ld%ld%ld%ld%ld%ld", &(self->total_program_size), &(self->resident_set_size), &(self->shared_pages), &(self->code_size), &(self->library_size), &(self->data_size), &(self->dirty_pages) );
	return 0;
}

typedef struct metrics {
	time_t wall_clock;
	statm memory;
	struct rusage usage; //maxrss[kB], utime[s.us], stime[s.us]
						 // http://linux.die.net/man/2/getrusage
} metrics;

int getmetrics_qd(metrics * self) { //quick and dirty (no actual err returned)
	self->wall_clock=time(NULL);
	getstatm_qd(&(self->memory));
	getrusage(RUSAGE_SELF,&(self->usage));
	return 0;
}

int getmetrics(metrics * self) {
	self->wall_clock=time(NULL);
	return (getstatm(&(self->memory)) +	getrusage(RUSAGE_SELF,&(self->usage))); // -1 err in getrusage, -2 or-4 err getstatm
}

void displaymetrics(metrics set) { //display on stdout (video)
	long pagesize = sysconf( _SC_PAGESIZE);
	printf("wall clock [s]                       \t%ld\n",set.wall_clock);
    printf("CPU user time [s]                    \t%ld.%06ld\n",set.usage.ru_utime.tv_sec,set.usage.ru_utime.tv_usec);
    printf("CPU system time [s]                  \t%ld.%06ld\n",set.usage.ru_stime.tv_sec,set.usage.ru_stime.tv_usec);
    printf("allocated memory [B]                 \t%ld\n",set.usage.ru_maxrss*1024);
	printf("resident data size (heap + stack) [B]\t%ld\n",set.memory.data_size*pagesize);
	printf("resident (code) text size [B]        \t%ld\n",set.memory.code_size*pagesize);
	printf("resident (code) library size [B]     \t%ld\n",set.memory.library_size*pagesize);
    printf("resident set size (data + code) [B]  \t%ld\n",set.memory.resident_set_size*pagesize);
	printf("total program size (res + sw) [B]    \t%ld\n",set.memory.total_program_size*pagesize);
    printf("shared pages                         \t%ld\n",set.memory.shared_pages);
    printf("dirty pages                          \t%ld\n",set.memory.dirty_pages);
}

void pipemetrics(metrics set) { //pipe output to stderr
	long pagesize = sysconf( _SC_PAGESIZE);
	fprintf(stderr,"wall clock [s]                       \t%ld\n",set.wall_clock);
    fprintf(stderr,"CPU user time [s]                    \t%ld.%06ld\n",set.usage.ru_utime.tv_sec,set.usage.ru_utime.tv_usec);
    fprintf(stderr,"CPU system time [s]                  \t%ld.%06ld\n",set.usage.ru_stime.tv_sec,set.usage.ru_stime.tv_usec);
    fprintf(stderr,"allocated memory [B]                 \t%ld\n",set.usage.ru_maxrss*1024);
	fprintf(stderr,"resident data size (heap + stack) [B]\t%ld\n",set.memory.data_size*pagesize);
	fprintf(stderr,"resident (code) text size [B]        \t%ld\n",set.memory.code_size*pagesize);
	fprintf(stderr,"resident (code) library size [B]     \t%ld\n",set.memory.library_size*pagesize);
    fprintf(stderr,"resident set size (data + code) [B]  \t%ld\n",set.memory.resident_set_size*pagesize);
	fprintf(stderr,"total program size (res + sw) [B]    \t%ld\n",set.memory.total_program_size*pagesize);
    fprintf(stderr,"shared pages                         \t%ld\n",set.memory.shared_pages);
    fprintf(stderr,"dirty pages                          \t%ld\n",set.memory.dirty_pages);
}

#endif
