/* 
Program:   accudate.c
Version:   1.0.1
Date:      10-DEC-2002
Author:    David Mathog, Biology Division, Caltech
email:     mathog@caltech.edu
url:       ftp://saf.bio.caltech.edu/pub/software/linux_or_unix_tools/
Copyright: 2002 David Mathog and California Institute of Technology (Caltech)

Description:

     "date" replacment which emits date and time down to milliseconds
     to stdout in format:

        2002-12-01 13:28:36.716
	
Changes:

  1.0.2 12-DEC-2002, Oops, removed a colon by accident going to 1.0.1
  1.0.1 10-DEC-2002, Added delta time capability.  -t0 switch
     emits current time to stdout as seconds microseconds. -df 
     switch followed by seconds microseconds emits delta time
     in ddd:hh:mm:ss.mmm and -dt emits it as seconds.mmm.  Note
     that seconds is converted to long long and then emitted,
     and microseconds to an int.  On 32 bit systems and above
     it should usually be possible to convert back and forth
     using this method.  long long isn't in older ANSI C though.
  1.0.0 06-DEC-2002, first release

License terms:
    You may run this program on any platform. You may
    redistribute the source code of this program subject to
    the condition that you do not first modify it in any way.
    You may  distribute binary versions of this program so long
    as they were compiled from unmodified source code.  There
    is no charge for using this software.  You may not charge
    others for the use of this software.

Miscellaneous:
    This should be relatively portable on Unix systems.  Compiles
    and builds on RedHat 7.3 and Solaris 8 with
    
    gcc -Wall -ansi -pedantic -o accudate accudate.c
    
    You may ignore the warnings about "long long" and "ll" not being
    supported - so long as _your_ compiler supports them.
    
    Since the program itself takes time to run the time point
    indicated by milliseconds is fairly arbitrary.  On our Sparc
    IIi machine a series of accudate commands emit times .01
    seconds apart, and on our Athlon 2200 RedHat machine the
    times are .002 seconds apart.  On future machines the run
    time should drop below 1 millisecond, at which point accudate
    will be accurate to the last decimal shown.

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>
 
/* definitions */
#define VERSTRING "1.0.2 12-DEC-2002"
enum timemodes {NORMAL,TZERO,DELTA,DELTASECONDS};

/* function prototypes */
void emit_help(void);
void insane(char *string);
int  lcl_strcasecmp(char *s1, char *s2);
void process_command_line_args(int argc,char **argv, enum timemodes *tmode, struct timeval *start_time);
void setzerotime(struct timeval *val,int *numarg,int argc,char **argv,char * label); 


/* functions */
void setzerotime(struct timeval *val,int *numarg,int argc,char **argv,char * label){
long long external_seconds;
int       external_microseconds;
      (*numarg)++;
      if( ( *numarg  >= argc ) || (argv[*numarg] == NULL)){
        (void) fprintf( stderr, "%s: missing seconds\n",label);
        exit(EXIT_FAILURE);
      }
      if(sscanf(argv[*numarg],"%lld",&external_seconds) != 1){
        (void) fprintf(stderr,"Bad integer argument/parameter [%s %s] \n",label,argv[*numarg]);
        exit(EXIT_FAILURE);
      }
      (*numarg)++;
      if( ( *numarg  >= argc ) || (argv[*numarg] == NULL)){
        (void) fprintf( stderr, "%s: missing microseconds\n",label);
        exit(EXIT_FAILURE);
      }
      if(sscanf(argv[*numarg],"%d",&external_microseconds) != 1){
        (void) fprintf(stderr,"Bad integer argument/parameter [%s %s] \n",label,argv[*numarg]);
        exit(EXIT_FAILURE);
      }
      if(external_seconds < 0){
        (void) fprintf(stderr,"Illegal value for seconds [%s %s] \n",label,argv[*numarg]);
        exit(EXIT_FAILURE);
      }
      if(external_microseconds < 0){
        (void) fprintf(stderr,"Illegal value for microseconds [%s %s] \n",label,argv[*numarg]);
        exit(EXIT_FAILURE);
      }
      val->tv_sec  = external_seconds;
      val->tv_usec = external_microseconds;
}

void emit_help(void){
(void) fprintf(stderr,"accudate command summary:\n\n");
(void) fprintf(stderr,"   [no args]     write time as YYYY-MM-DD:hh:mm:ss.lll  to stdout\n");
(void) fprintf(stderr,"   -t0           write time as \"seconds microseconds\" to stdout\n");
(void) fprintf(stderr,"   -df seconds microseconds\n");
(void) fprintf(stderr,"                 write elapsed time as \"DDDD-hh:mm.ss.lll\" to stdout\n");
(void) fprintf(stderr,"   -ds seconds microseconds\n");
(void) fprintf(stderr,"                 write elapsed time as \"ssssss.lll\" to stdout\n");
(void) fprintf(stderr,"   -h            print this help message (also -help --h --help -? --?)\n");
(void) fprintf(stderr,"   -i            emit version, copyright, license and contact information\n\n");
exit(EXIT_FAILURE);
}

void insane(char *string){
 (void) fprintf(stderr,"%s\n",string);
 exit(EXIT_FAILURE);
}

int  lcl_strcasecmp(char *s1, char *s2){
int c1;
int c2;
  for(; ;s1++,s2++){
    c1=toupper(*s1);
    c2=toupper(*s2);
    if(c1 < c2)return -1;
    if(c1 > c2)return  1;
    if(c1 == 0)return  0;  /*c2 also is 0 in this case */
  }
}

void process_command_line_args(int argc,char **argv, enum timemodes *tmode, struct timeval *start_time){

 int numarg=0;
 
  *tmode=NORMAL;

 while( ++numarg < argc){
    if( (lcl_strcasecmp(argv[numarg], "-h")==0)     ||
        (lcl_strcasecmp(argv[numarg], "--h")==0)    ||
        (lcl_strcasecmp(argv[numarg], "-?")==0)     ||
        (lcl_strcasecmp(argv[numarg], "--?")==0)    ||
        (lcl_strcasecmp(argv[numarg], "-help")==0)  ||
        (lcl_strcasecmp(argv[numarg], "--help")==0) ){
      emit_help();
    }
    else if(lcl_strcasecmp(argv[numarg], "-i")==0){
      (void)fprintf(stderr,"Version:   %s\n",VERSTRING);
      (void)fprintf(stderr,"bugs to:   mathog@caltech.edu\n");
      (void)fprintf(stderr,"Copyright: 2002 David Mathog and California Institute of Technology\n");
      (void)fprintf(stderr,"License terms:\n");
      (void)fprintf(stderr,"    You may run this program on any platform. You may\n");
      (void)fprintf(stderr,"    redistribute the source code of this program subject to\n");
      (void)fprintf(stderr,"    the condition that you do not first modify it in any way.\n");
      (void)fprintf(stderr,"    You may  distribute binary versions of this program so long\n");
      (void)fprintf(stderr,"    as they were compiled from unmodified source code.  There\n");
      (void)fprintf(stderr,"    is no charge for using this software.  You may not charge\n");
      (void)fprintf(stderr,"    others for the use of this software.\n");
      exit(EXIT_SUCCESS);
    }
    else if(lcl_strcasecmp(argv[numarg], "-ds")==0){
      setzerotime(start_time,&numarg,argc,argv,"-ds");
      *tmode=DELTASECONDS;
      continue;
    }
    else if(lcl_strcasecmp(argv[numarg], "-df")==0){
      setzerotime(start_time,&numarg,argc,argv,"-df");
      *tmode=DELTA;
      continue;
    }
    else if(lcl_strcasecmp(argv[numarg], "-t0")==0){
      *tmode=TZERO;
      continue;
    }
    else
      (void) fprintf(stderr,"Unknown command line argument: %s\n",argv[numarg]);
      exit(EXIT_FAILURE);
      continue;
    }
}

int main(int argc, char *argv[]){

        struct  tm  *time_structure;
        time_t time_val;
	struct timeval start_time;
	struct timeval thetv;
        struct timeval *tp=&thetv;
	int millisec;
	enum timemodes tmode;
	long long external_seconds;
	time_t    delta_s;
	int       external_microseconds,delta_us;
	int       days,hours,minutes,seconds;
        
        process_command_line_args(argc,argv, &tmode, &start_time);
	
        (void) time(&time_val);
        time_structure = localtime(&time_val);
        tp->tv_sec     = 0;
	tp->tv_usec    = 0;
        if (gettimeofday(tp,(struct timezone *)NULL) != 0) {
                (void) fprintf(stderr,"Error calling gettimeofday\n");
                exit(EXIT_FAILURE);
        }
	millisec = tp->tv_usec/1000;
	
	switch (tmode){
	  case NORMAL:
           (void) fprintf(stdout,"%.4d-%.2d-%.2d:%.2d:%.2d:%.2d.%.3d\n",
               1900+time_structure->tm_year,
               1+time_structure->tm_mon,
               time_structure->tm_mday,
               time_structure->tm_hour,
               time_structure->tm_min,
               time_structure->tm_sec,
               millisec);
	    break;
	  case TZERO:
	   external_seconds = tp->tv_sec;      /* implicit conversion, usually safe */
	   external_microseconds= tp->tv_usec; /* implicit conversion, usually safe */
           (void) fprintf(stdout,"%lld %d\n",external_seconds,external_microseconds);
	    break;
	  case DELTA:
            delta_us = tp->tv_usec - start_time.tv_usec;
	    if(delta_us < 0)delta_us += 1000000;
	    millisec=delta_us/1000;
	    delta_s= time_val - start_time.tv_sec;
	    external_seconds=delta_s;
	    days=external_seconds/(24*60*60);
	    if(days>0)external_seconds -= days*24*60*60;
	    hours=external_seconds/(60*60);
	    if(hours>0)external_seconds -= hours*60*60;
	    minutes=external_seconds/60;
	    if(minutes>0)external_seconds -= minutes*60;
	    seconds=external_seconds;
            (void) fprintf(stdout,"%.4d-%.2d:%.2d:%.2d.%.3d\n",days,hours,minutes,seconds,millisec);
	    break;
	  case DELTASECONDS:
            delta_us = tp->tv_usec - start_time.tv_usec;
	    if(delta_us < 0)delta_us += 1000000;
	    millisec=delta_us/1000;
	    delta_s= time_val - start_time.tv_sec;
	    external_seconds=delta_s;
           (void) fprintf(stdout,"%7.7lld.%.3d\n",external_seconds,millisec);
	    break;
	}
        exit(EXIT_SUCCESS);
}                                                                       

