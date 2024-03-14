/*
 *  NEMOVAR:   display a NEMO variable from the $NEMOVAR file (or future shared memory)
 */

#include <nemo.h>
#include <table.h>
#include <extstring.h>


string defv[] = {
  "var=\n             Name of  variable(s)",
  "value=\n           Optional new value",
  "show=f\n           Show all variables",
  "VERSION=0.2\n      7-mar-2024 PJT",
  NULL,
};

string usage="display NEMO variable(s)";


void nemo_main()
{
  string var = getparam("var");
  bool hasVar = hasvalue("var");
  string value = getparam("value");  
  bool hasVal = hasvalue("value");
  string nemovar = getenv("NEMOVAR");
  stream tstr = stropen(nemovar,"r");
  table *t = table_open(tstr, -1);      // read table as a dumb set of lines
  int i, nrows = t->nr;
  string s;

  dprintf(1,"nemovar: %s\n", nemovar);
  

  if (!hasVar && !hasVal) { 
    for (i=0; i<nrows; i++)
	    printf("%s\n", table_row(t, i));
    return;
  }

  //@TODO - Should use hash in the future?
  if (hasVar && !hasVal) { 
    // writing code might need flock() to lock file
    // loop over lines and find the keyword, then display it
    bool found_var = FALSE;
    for (i=0; i<nrows; i++) {
      s = table_row(t, i);
      if (strncmp(s,var,strlen(var))==0) {
          printf("%s\n",s+strlen(var)+1);      
          found_var = TRUE;  
      }    
    }

    if (!found_var) {
      error("Variable %s not found\n",var);
    }
  }


  if (hasVar && hasVal) {
    stream write_stream = stropen(nemovar,"w!");
    bool found_var = FALSE;
    
    for (i=0; i<nrows; i++) {
      s = table_row(t, i);
      if (strncmp(s,var,strlen(var))==0) {
        fprintf(write_stream,"%s=\"%s\"\n",var,value);    
        found_var = TRUE;  
      } else {
        fprintf(write_stream,"%s\n",s);
      }
    }

    if (!found_var) {
      fprintf(write_stream,"%s=\"%s\"\n",var,value);    
    }

    

    
  }

  // loop over lines and find the keyword, then display it
  
  // printf("%s= TBD\n",var);

}
