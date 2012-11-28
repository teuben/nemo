/*
 * COMMAND:   interactive command parser
 *
 *  MIR_READLINE:  use miriad's save/load feature (not implemented)
 *
 *
 * TODO:
 *    - command_read
 */


#include <nemo.h>
#include <extstring.h>
#include <table.h>
#include <command.h>
#ifdef HAVE_LIBREADLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

extern string *burststring(string,string);
extern void freestrings(string *);
extern int nemo_file_lines(string, int);

#define VALID_TYPES "irs."

static int todo_readline=0;
static char prompt[64];

/*
 *
 */

command *command_init(string name) 
{
  command *c;
  c = (command *) allocate(sizeof(command));
  c->name = strdup(name);
  c->ncmd = 0;
#ifdef HAVE_LIBREADLINE
    dprintf(0,"[readline library installed]\n");
#endif
#if (MIR_READLINE==1)
  if (todo_readline==0)
    dprintf(0,"[readline library installed]\n");
  if (todo_readline) {
    todo_readline=1;
    ini_readline();
  } else
    todo_readline++;
#endif
  return c;
}

/*
 *
 */

void command_register(command *c, string cmd, string type, string help)
{
  int i, clen, n = c->ncmd;
  char *cp;

  for (i=0; i<n; i++) {     /* check if not already known */
    if (streq(c->cmd[i],cmd))
      error("Command %s was already registered",cmd);
  }
  if (n==MAXCMD) error("Too many commands (%d)",MAXCMD);
  c->cmd[n] = strdup(cmd);
  c->help[n] = strdup(help);
  c->type[n] = strdup(type);
  c->nargs[n] = strlen(type);
  c->ncmd++;

  cp = c->type[n];
  clen = strlen(cp);
  if (clen>0) {
    for (i=0; i<clen; i++) {
      if (strchr(VALID_TYPES,cp[i]) == NULL)
	error("Command arguments to %s are %s, should only contain \"%s\"",
	      cmd,type,VALID_TYPES);
      if (cp[i] == '.') {
	*cp = 0;
	break;
      }
    }      
    c->nargs[n] = strlen(c->type[n]);
  }
}

static string cmds[MAXCMD];
static int ncmds = 0;
static int curcmd = 0;

void command_read(string fname)
{
  stream instr;
  int n;
  if (fname==0 || *fname==0) {
    warning("command_read: blank filename");
    return;
  }

  n = nemo_file_lines(fname,0);
  dprintf(0,"Found %d lines in %s\n",n,fname);
  if (n < 1) return;

  instr = stropen(fname,"r");
  ncmds = curcmd = 0;
  while ( (cmds[ncmds] = getaline(instr) )) {
    ncmds++;
    if (ncmds == MAXCMD) error("Too many commands in file for now: %d",MAXCMD);
  }
  strclose(instr);
  dprintf(0,"Read %d command lines from %s\n",ncmds,fname);
}

/*
 *  TODO:
 *   - breaking up line into argv[] should listen to quotes
 */

string *command_get(command *c)
{
  char *s, line[1024];
  string *argv;
  int i, icmd, n, na;

  sprintf(prompt,"%s> ",c->name);
  
 again:

  if (ncmds) {
    strcpy(line,cmds[curcmd]);
    curcmd++;
    dprintf(0,"COMMAND_READ-%d> %s\n",curcmd,line);
    if (curcmd==ncmds) {
      dprintf(1,"This was the last command; freeing buffers\n");
      for (i=0; i<ncmds; i++)
	free(cmds[i]);
      curcmd = ncmds = 0;
    }
  } else {
#ifdef HAVE_LIBREADLINE
    for(;;) {
      if ((s=readline(prompt)) != (char *)NULL) {
	/* stripwhite(s); */
	if (*s) {
	  strcpy(line,s);
	  free(s);
	  s=0;
	  break;
	}
      }
      if (s) free(s);
    }
    add_history(line);
#else
    printf("%s",prompt);
    fflush(stdout);
    clearerr(stdin);
    if (fgets(line,1024,stdin) == NULL)
      return NULL;
#endif
  }
  n = strlen(line);
  if (n>0 && line[n-1]=='\n')
    line[n-1]=0;

  if (line[0] == 0)                             /* empty line */
    goto again;

  if (line[0] == '#')                           /* comment line */
    goto again;

  if (line[0] == '.')                           /* internal quit command */
    return NULL;

  if (line[0] == '?') {                         /* internal help command */
    dprintf(0,"Internal commands: \n");
    dprintf(0,"?          this help\n");
    dprintf(0,"!CMD       shell escape\n");
    dprintf(0,"#          comment line\n");
    dprintf(0,".          quit\n");
    dprintf(0,"<FILE      read commands from FILE\n\n"); 
    dprintf(0,"Registered commands: \n");
    for (i=0; i<c->ncmd; i++) 
      dprintf(0,"%-10s %s\n",c->cmd[i],c->help[i]);
    dprintf(0,"\n");
    goto again;
  }

  if (line[0] == '<') {                         /* internal read command */
    for (i=1; line[i]; i++)
      if (!isspace(line[i])) break;
    icmd = i;
    for (; line[i]; i++)
      if (isspace(line[i])) line[i] = 0;
    command_read(&line[icmd]);
    goto again;
  }

  if (line[0] == '!') {                         /* shell escape */
    system(&line[1]);
    goto again;
  }

  /* having come here, parse input into an argv[] vector of strings */

  argv = burststring(line," \n\t");
  na = xstrlen(argv,sizeof(string))-2;
  if (na<0) goto again;
  for (i=0; i<c->ncmd; i++) {
    if (streq(c->cmd[i],argv[0])) {
      dprintf(1,"Found matching command %s, needs %d, got %d\n",
	      argv[0],c->nargs[i],na);
      if (na < c->nargs[i]) {
	warning("Not enough arguments for %s (need %d)",argv[0],c->nargs[i]);
	freestrings(argv);
	goto again;
      }
      return argv;         /* return this argv[] vector to the user */
    }
  }
  warning("%s: not a valid command, try '?'",argv[0]);
  freestrings(argv);
  goto again;
}

/*
 *
 */

void command_close(command *c)
{
  dprintf(1,"Closing command[%s]\n",c->name);
#if (MIR_READLINE==1)
  todo_readline--;
  if (todo_readline==0)
    end_readline();
#endif
}



/* -------------------------------------------------------------------------------- */

#ifdef TESTBED

string defv[] = {
  "n=10\n         some integer",
  "name=test123\n basename and prompt for command processor",
  "VERSION=1\n    22-dec-03 PJT",
  NULL,
};

string usage = "testing command processing";

void do_a(int a) {
  printf("a::a=%d\n",a);
}

void do_b(int a, double b) {
  printf("b::a=%d\n",a);
  printf("b::b=%g\n",b);
}
void do_c(string a) {
}

void do_d1(int a) {
  printf("d1::a=%d\n",a);
}

void do_d2(int a , double b) {
  printf("d2::a=%d\n",a);
  printf("d2::b=%g\n",b);
}

void do_e(void) {
  printf("e::\n");
}

nemo_main()
{
  int i, n = getiparam("n");
  int na, ivar;
  
  string name = getparam("name");
  command *cmd;
  string *argv;
  
  cmd = command_init(name);
  command_register(cmd,"a","i",   "needs just an integer");
  command_register(cmd,"b","ir",  "needs integer and real");
  command_register(cmd,"c","s",   "just a string");
  command_register(cmd,"d","i.",  "integer, and optional remaining");
  command_register(cmd,"e","",    "needs nothing");
  command_register(cmd,"quit","", "a long quit");

  while((argv=command_get(cmd))) {                  /* loop getting arguments */
    na = xstrlen(argv,sizeof(string))-1;
    dprintf(0,"Command: %s has argc=%d \n",argv[0],na);
    for (i=0; i<na; i++)
      dprintf(1,"argv[%d] = %s\n",i,argv[i]);

    if (streq(argv[0],"?")) {
      warning("command not understood");
    } else if (streq(argv[0],"quit")) {
      break;
    } else if (streq(argv[0],"a")) {
      do_a(natoi(argv[1]));
    } else if (streq(argv[0],"b")) {
      do_b(natoi(argv[1]), natof(argv[2]));
    } else
      warning("command %s not understood yet",argv[0]);

    freestrings(argv);
  }
  dprintf(0,"All done with \"%s\"\n",name);

}

#endif
