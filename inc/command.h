/*
 * command.h 
 *  
 * 
 *
 */

#define MAXCMD 1024

typedef struct _command {
  string name;          /* name of the command parser */
  int ncmd;             /* number of commands stored so far */
  string cmd[MAXCMD];   /* name of the valid commands */
  char *type[MAXCMD];   /* i=int d=double s=string .=optional */
  int nargs[MAXCMD];    /* number of required arguments of each command (0=void) */
  string help[MAXCMD];  /* help string */
  
} command;



/* command.c */
command *command_init(string name);
void command_register(command *c, string cmd, string argtypes, string help);
string *command_get(command *c);
void command_read(string filename);
void command_close(command *c);

