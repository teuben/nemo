/*
** tr2latex - troff to LaTeX converter
** $Id$
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
** 
** This module contains the function prototypes.
*/

int		CAP_GREEK (char *w);
void	GR_to_Greek (char *w, char *ww);
int		main (int argc, char *argv[]);
char *	alternate (char *pin, char *pout, char *w);
char *	do_table (char *pin, char *pout, int *offset);
char *	end_env (char *pout);
void	envoke_stat (int par);
void	errexit (int exitcode);
char *	flip (char *pout, char *w);
char *	flip_twice (char *pout, char *w, char *ww);
int		get_N_lines (char *pin, char *w, int N);
int		get_allargs (char *pin, char ***ppw, int rec);
int		get_arg (char *pin, char *w, int rec);
void	get_brace_arg (char *buf, char *w);
int		get_defword (char *pin, char *w, int *illegal);
int		get_line (char *pin, char *w, int rec);
int		get_multi_line (char *pin, char *w);
int		get_mydef (char *pin, char *w);
int		get_no_math (char *pin, char *w);
char *	get_over_arg (char *pin, char *ww);
int		get_ref (char *pin, char *w);
void	get_size (char *ww, struct measure *PARAMETER);
int		get_string (char *pin, char *w, int rec);
int		get_sub_arg (char *pin, char *w);
int		get_table_entry (char *pin, char *w, int tab);
int		get_till_space (char *pin, char *w);
int		getdef (char *pin, char *ww);
void	getopts (int *p_argc, char *argv[]);
int		getword (char *pin, char *w);
int		is_def (char *w);
int		is_flip (char *w);
int		is_forbid (char *w);
int		is_mathcom (char *w, char *ww);
int		is_mydef (char *w);
int		is_troff_mac (char *w, char *ww, int *arg, int *par);
void	parse_units (char *ww, int *sign, int *units, float *value);
void	process (FILE *in_file, char *f_name, char *pin, char *pout);
int		similar (char *w);
char *	skip_line (char *pin);
int		skip_white (char *pin);
char *	strapp (char *s, char *tail);
char *	strsave (char *s);
void	tmpbuf (FILE *in, char *buffer);
void	troff_tex (char *pin, char *pout, int mid, int rec);
void	usage (int exitcode);
