/*
 * COMMAND:   interactive command parser
 *
 *
 */


#include <nemo.h>
#include <command.h>

extern string *burststring(string,string);
extern void freestrings(string *);


command *command_init(string name) 
{
  command *c;
  c = (command *) allocate(sizeof(command));
  c->name = strdup(name);
  c->ncmd = 0;
  return c;
}

#define VALID_TYPES "irs."

void command_register(command *c, string cmd, string type)
{
  int i, clen, n = c->ncmd;
  char *cp;

  for (i=0; i<n; i++) {     /* check if not already known */
    if (streq(c->cmd[i],cmd))
      error("Command %s was already registered",cmd);
  }
  if (n==MAXCMD) error("Too many commands (%d)",MAXCMD);
  c->cmd[n] = strdup(cmd);
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
      if (cp[i] == ".") {
	*cp = 0;
	break;
      }
    }      
    c->nargs[n] = strlen(c->type[n]);
  }
}

string *command_get(command *c)
{
  char line[1024];
  string *sp;
  int i, na;
  
 again:
  printf("%s> ",c->name);
  fflush(stdout);
  clearerr(stdin);
  if (fgets(line,1024,stdin) == NULL)
    return NULL;
  if (line[0] == '?') {
    dprintf(0,"Valid commands: ");
    for (i=0; i<c->ncmd; i++) 
      dprintf(0,"%s ",c->cmd[i]);
    dprintf(0,"\n");
    goto again;
  }

  sp = burststring(line," \n");
  na = xstrlen(sp,sizeof(string))-2;
  for (i=0; i<c->ncmd; i++) {
    if (streq(c->cmd[i],sp[0])) {
      dprintf(0,"Found matching command %s, needs %d, got %d\n",
	      sp[0],c->nargs[i],na);
      if (na < c->nargs[i]) {
	warning("Not enough arguments for %s (need %d)",sp[0],c->nargs[i]);
	freestrings(sp);
	goto again;
      }
      return sp;
    }
  }
  warning("%s: not a valid command",sp[0]);
  freestrings(sp);
  goto again;
}

void command_close(command *c)
{
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
  command_register(cmd,"a","i");
  command_register(cmd,"b","ir");
  command_register(cmd,"c","s");
  command_register(cmd,"d","i.");
  command_register(cmd,"e","");
  command_register(cmd,"quit","");

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
