/*
 *  RUN:   generic edit config file and run a program in a (new) run directory
 *
 *    all arguments after "--" on the commandline are passed onto the editor
 *    which then replaces all "key=val" instances as directed.
 *
 */

#include <nemo.h>

string defv[] = {
  "exe=???\n          Executable to run",
  "ini=???\n          Ini file to edit",
  "dir=???\n          Directory in which to run",
  "log=\n             Logfile, if needed (> logfile)",
  "stdin=f\n          If true, read ini from stdin (< inifile)",
  "VERSION=0.1\n      8-may-2011 PJT",
  NULL,
};

string usage="Generic frontend to edit an ini file and run programs";

string cvsid="$Id$";


void nemo_main(void);
void goto_rundir(string name);
void make_rundir(string name);
void run_program(string cmd);
void edit_file(string inifile, int argc, string *argv);

void nemo_main()
{
  string exefile = getparam("exe");
  string inifile = getparam("ini");
  string dirname = getparam("dir");
  string logfile = getparam("log");
  bool usestdin = getbparam("stdin");
  char dname[256], cmd[256];
  int i, argc;
  string *argv;


  make_rundir(dirname);

#if 0
  sprintf(dname,"%s/%s",outdir,parfile);
  parstr = stropen(dname,"w");
  fprintf(parstr,"%d\n",getiparam("numbs"));
  fprintf(parstr,"%s\n",getparam("aexpn"));
  fprintf(parstr,"%s\n",getparam("adiv"));
  fprintf(parstr,"%s\n",getparam("om0"));
  fprintf(parstr,"%s\n",getparam("vsca"));
  fprintf(parstr,"%s\n",getparam("hubble"));
  fprintf(parstr,"%d\n",getiparam("mond"));
  fprintf(parstr,"%d\n",getiparam("freq"));
  fprintf(parstr,"%d\n",getiparam("brand"));
  fprintf(parstr,"%s\n",getparam("au0"));
  strclose(parstr);
#endif

  sprintf(cmd,"cp %s %s/%s", inifile, dirname, inifile);
  run_program(cmd);

  goto_rundir(dirname);

  argv = getargv(&argc);
  edit_file(inifile, argc, argv);

  sprintf(cmd,"%s ",exefile);
  if (usestdin) {
    strcat(cmd,"<"); 
    strcat(cmd,inifile);
    strcat(cmd," "); 
  }
  if (hasvalue("log")) {
    strcat(cmd,">"); 
    strcat(cmd,logfile);
  }
  run_program(cmd);

}

void goto_rundir(string name)
{
  dprintf(0,"CHDIR: %s\n",name);
  if (chdir(name))
    error("Cannot change directory to %s",name);
}

void make_rundir(string name)
{
  dprintf(0,"MKDIR: %s\n",name);
  if (mkdir(name, 0755))
    warning("Run directory %s already exists",name);
}

void run_program(string cmd)
{
  dprintf(0,"EXEC: %s\n",cmd);
  system(cmd);
}

void edit_file(string inifile, int argc, string *argv)
{
  int i;

  for (i=1; i<argc; i++) {
    dprintf(0,"EDIT: %s\n",argv[i]);
  }
}
