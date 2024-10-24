/*
 *  TABGEN:     generate a table with (random) numbers, purely for benchmarking
 */

#include <nemo.h>

string defv[] = {
    "out=???\n     Output table",
    "nr=10\n       Number of rows",
    "nc=5\n        Number of columns",
    "mode=1\n      Mode (1=uniform, 2=normal 3=constant 4=linear, 5=linear)",
    "seed=123\n    Random seed",
    "fmt=%g\n      Format statement for output",
    "sep=s\n       column separator (s=space t=tab c=comma v=vertical bar)",
    "addrow=f\n    Add row number (1=first) as first column?",
    "VERSION=0.8\n 14-feb-2024 PJT",
    "header=None\n     Add a dummy header in a given style [ecsv, ipac]",
    NULL,
};

string usage="Create a table with (random) numbers";


typedef real (*my_real_proc)(real,real);

local real constant(real c1, real c2)
{
  return c2;
}

local real linear(real c1, real c2)
{
  static real count = 0;
  count = count + 1;
  return count;
}

string add_header(string header)
{
  char *token;
  char *result = malloc(100);
  result[0] = '\0';
  int first = 1;
  int arr_len = 0;

  token = strtok(header, ",");

  strcat(result, "| ");
  while (token != NULL) {
        if (!first) {
            strcat(result, " | ");
        }
        strcat(result, token);
        first = 0;
        arr_len++;
        token = strtok(NULL, ",");
  }

  strcat(result, " |");
  strcat(result, "\n");
  strcat(result, "| double");

  for (int i = 1; i < arr_len; i++) {
        strcat(result, " | double");
  }

  strcat(result, " |");

  return result;
}

void nemo_main()
{
  stream ostr = stropen(getparam("out"),"w");
  int i, nr = getiparam("nr");
  int j, nc = getiparam("nc");
  int mode = getiparam("mode");
  int amode = ABS(mode);
  int seed = init_xrandom(getparam("seed"));
  string fmt = getparam("fmt");
  my_real_proc my_random = NULL;
  string seps = getparam("sep");
  bool Qrow = getbparam("addrow");
  char sep[8];
  string header = getparam("header");

  fprintf(ostr, "%s\n", add_header(header));

  if (seps[0] == 'c') strcpy(sep,",");
  else if (seps[0] == 's') strcpy(sep," ");
  else if (seps[0] == 't') strcpy(sep,"\t");
  else if (seps[0] == 'v') strcpy(sep,"|");
  else strcpy(sep,seps);

  dprintf(1,"seed=%d\n",seed);

  if (amode == 1) {
    dprintf(1,"Uniform values between 0 and 1\n");
    my_random = xrandom;
  } else if (amode == 2) {
    dprintf(1,"Uniform values between 0 and 1\n");
    my_random = grandom;
  } else if (amode == 3) {
    dprintf(1,"Uniform values between 0 and 1\n");
    my_random = constant;
  } else if (amode == 4) {
    dprintf(1,"Uniform values between 0 and 1\n");
    my_random = linear;
  } else
    error("Illegal mode=%d",mode);

  if (mode < 0) {
    // only produce random numbers, testing production rate
    real *x = (real *) allocate(nc*nr*sizeof(real));
    int k=0;
    for (i=0; i<nr; i++)
      for (j=0; j<nc; j++)
        x[k++] = my_random(0.0,1.0);
  } else
    // output to file as well
    for (i=0; i<nr; i++) {
      if (Qrow) fprintf(ostr,"%d ", i+1);
      for (j=0; j<nc; j++) {
        if (j>0) fprintf(ostr, "%s", sep);
        fprintf(ostr,fmt, my_random(0.0, 1.0));
      }
      fprintf(ostr,"\n");
    }

}