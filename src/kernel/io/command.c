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
  return c;
}

void command_register(command *c, string cmd, string argtypes)
{

}

string *command_get(command *c)
{
  char line[1024];
  string *sp;
  
  printf("%s> ",c->name);
  fflush(stdout);
  clearerr(stdin);
  if (fgets(line,1024,stdin) == NULL)
    return NULL;

  sp = burststring(line," ");
  return sp;
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
  int n = getiparam("n");
  int na;
  string name = getparam("name");
  command *cmd;
  string *argv;
  
  cmd = command_init(name);
  command_register(cmd,"a","i");
  command_register(cmd,"b","id");
  command_register(cmd,"c","s");
  command_register(cmd,"d","i.");
  command_register(cmd,"e","");

  while((argv=command_get(cmd))) {                  /* loop getting arguments */
    na = xstrlen(argv,sizeof(string))-1;
    dprintf(0,"Command: %s has %d argc\n",argv[0],na);

    if (streq(argv[0],"?")) {
      warning("command not understood");
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
