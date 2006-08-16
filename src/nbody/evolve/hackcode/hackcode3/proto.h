/* code.c */
void startrun(void);
void make_testdata(bool cencon);
void stepsystem(void);

/* code_io.c */
void inputdata(string file);
void initoutput(void);
void stopoutput(void);
void output(void);
void savestate(string file);
void restorestate(string file);

/* grav.c */
void hackgrav(bodyptr p);
void hackwalk(proc sub);

/* hackforce.c */
int  input_data(void);
int  read_snapshot(bodyptr *btab_ptr, int *nobj_ptr, stream instr);
void force_calc(void);
void out_result(void);
void write_snapshot(void);

/* load.c */
void maketree(bodyptr btab, int nbody);

/* util.c */
void pickvec(vector x, bool cf);
