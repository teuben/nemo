/*
 * command.h 
 *  
 * 
 *
 */

#define MAXCMD 128

typedef struct _command {
  string name;          /* name of the command parser */
  string cmd[MAXCMD];   /* name of the valid commands */
  int nargs[MAXCMD];    /* number of arguments of each command */
  char *type[MAXCMD];   /* i=int d=double s=string .=optional */
} command;



/* command.c */
command *command_init(string name);
void command_register(command *c, string cmd, string argtypes);
string *command_get(command *c);
void command_close(command *c);
