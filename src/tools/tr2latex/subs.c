/*
** tr2latex - troff to LaTeX converter
** $Id$
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
** 
** Module: subs.c
** 
** These subroutines do (in general) small things for the translator.
** They appear in alphabetical order and their names are unique in the
** first six characters.
*/

#include	"setups.h"
#include	"protos.h"
#include	"simil.h"
#include	"greek.h"
#include	"flip.h"
#include	"forbid.h"
#include	"maths.h"
#include	"macros.h"

extern def_count;
extern mydef_count;

/* compile-time counting of elements */
int GRK_count = (sizeof(GRK_list)/sizeof(GRK_list[0]));
int sim_count = (sizeof(sim_list)/sizeof(sim_list[0]));
int flip_count = (sizeof(flip_list)/sizeof(flip_list[0]));
int forbd_count = (sizeof(forbid)/sizeof(forbid[0]));
int mathcom_count = (sizeof(math)/sizeof(struct math_equiv));
int macro_count = (sizeof(macro)/sizeof(struct macro_table));
	
/*
** alternate fonts (manual macro)
*/
char *alternate (char *pin, char *pout, char *w)
{
	int which = 1;
	char font[MAXWORD], font1[MAXWORD], font2[MAXWORD],
	     ww[MAXWORD], tmp[MAXWORD];
	
	tmp[0] = EOS;
	switch (w[1]) {
	case 'R':	strcpy (font1,"\\rm"); break;
	case 'I':	strcpy (font1,"\\it"); break;
	case 'B':	strcpy (font1,"\\bf"); break;
	}
	switch (w[2]) {
	case 'R':	strcpy (font2,"\\rm"); break;
	case 'I':	strcpy (font2,"\\it"); break;
	case 'B':	strcpy (font2,"\\bf"); break;
	}
	
	strcpy (font, font1);
	while (*pin != '\n' && *pin != EOS) {
		pin += get_arg (pin, ww, 1);
		if (which == 1) {
			sprintf(tmp,"{%s %s}", font1, ww);
			which = 2;
		}
		else {
			sprintf (tmp,"{%s %s}", font2, ww);
			which = 1;
		}
		pout = strapp (pout, tmp);
		while (*pin == ' ' || *pin == '\t')
			pin++;
	}
	return (pout);
}



/*
** check if w is in the GREEK list
*/
int CAP_GREEK (char *w)
{
	int i;
	
	for (i = 0; i < GRK_count; i++) {
		if (strcmp (GRK_list[i], w) == 0)
			return (1);
	}
	return (-1);
}


/*
** translate table
** arguments:
**  pin
**  pout
**	offset    amount to offset pin
*/
char *do_table (char *pin, char *pout, int *offset)
{
	char w [MAXWORD], ww [MAXWORD], format [MAXWORD], tmp [MAXWORD];
	char *ptr;
	int i, len, columns=0;
	int tab = '\t';						/* default tab */
	
	tmp[0] = EOS;
	ptr = pin;							/* remember where we started */
	len = get_line (pin, w, 0);
	if (w [strlen (w) - 1] == ';') {	/* options */
		pin += len;
		if (strncmp (w, "tab", 3) == 0)	/* get the tab charecter */
			tab = w[4];		/* expect something like tab(&); */
		pin = skip_line(pin);
	}
	while (*pin != EOS) {				/* get the LAST format line */
		len = get_line (pin, w, 0);
		if (w[strlen (w) - 1] != '.')	/* not a format line */
			break;
		pin += len;
		for (i = 0; i < len - 1; i++) {
			if (isspace (w[i]))
				continue;
			if (w[i] == 'l')
				format[columns] = 'l';
			else if (w[i] == 'r')
				format[columns] = 'r';
			else
				format[columns] = 'c';
			columns++;
		}
	}
	if (columns == 0) {
		fprintf (stderr, "Sorry, I cannot do tables without a format line\n\
Doing plain translation of table, lines will be commented\n\
You need to fix it yourself\n");
		while (*pin != EOS) {
			(void) getword (pin, w);
			if (strcmp (w,".TE") ==  0) {
				pin += 4;
				break;
			}
			pin += get_line (pin, w, 1);
			*pout++ = '%';
			pout = strapp (pout, w);
			pout = strapp (pout, "\n");
			pin++;		/* skip the \n */
		}
		*offset = pin - ptr;
		return (pout);
	}
	format[columns] = EOS;
	sprintf (tmp, "\\par\n\\begin{tabular}{%s}\n",format);
	pout = strapp (pout, tmp);
	
	while (*pin != EOS) {
		for (i = 0; i < columns - 1; i++) {
			(void) getword (pin, w);
			if (i == 0 && (strcmp (w, "\n") == 0 || strcmp (w, "_") == 0)) {
				pin++;
				i--;
				continue;
			}
			if (strcmp (w, ".TE") == 0) {
				pin += 4;
				if (i == 0) {
					pout -= 3;	/* take back the \\ and the \n */
					*pout = EOS;
				}
				pout = strapp(pout,"\n\\end{tabular}\n\\par\n");
				*offset = pin - ptr;
				return (pout);
			}
			pin += get_table_entry(pin,w,tab);
			pin ++;		/* skip tab */
			troff_tex(w,ww,0,1);
			sprintf(tmp,"%s & ",ww);
			pout = strapp(pout,tmp);
		}
		(void) getword (pin,w);
		if (strcmp(w,".TE") == 0) {
			fprintf(stderr,"Oops! I goofed. I told I you I am not very good at tables\nI've encountered an unexpected end for the table\n\
You need to fix it yourself\n");
			pin += 4;
			pout = strapp(pout,"\\end{tabular}\n\\par\n");
			*offset = pin - ptr;
			return(pout);
		}
		pin += get_table_entry(pin,w,'\n');
		pin++;		/* skip tab */
		troff_tex(w,ww,0,1);
		pout = strapp (pout, ww);
		pout = strapp (pout, "\\\\\n");
	}
	fprintf (stderr, "Oops! I goofed. I told I you I am not very good at tables\n\
File ended and I haven't finished the table!\n\
You need to fix it yourself\n");
	*offset = pin - ptr;
	pout = strapp (pout, "\\end{tabular}\n\\par\n");
	return (pout);
}



/*
** end current environment
*/
char *end_env (char *pout)
{
	if (IP_stat) {
		IP_stat = 0;
		pout = strapp (pout, "\\end{IPlist}");
	}
	if (QP_stat) {
		QP_stat = 0;
		pout = strapp (pout, "\\end{quotation}");
	}
	if (TP_stat) {
		TP_stat = 0;
		pout = strapp (pout, "\\end{TPlist}");
	}
	return(pout);
}


/*
** set flag for current environment
*/
void envoke_stat (int par)
{
	
	switch(par) {
	case 2:
		IP_stat = 1;
		break;
	case 3:
		TP_stat = 1;
		break;
	case 4:
		QP_stat = 1;
		break;
	default:
		break;
	}
}



/*
** do the flipping
*/
char * flip (char *pout, char *w)
{
	int lb=0, rb=0;
	char ww[MAXWORD], tmp[MAXWORD];
	
	ww[0] = EOS;
	tmp[0] = EOS;
	pout--;
	while (isspace (*pout))
		pout--;
	while (1) {
		if (*pout == '{') {
			lb++;
			if (lb > rb)
				break;
		}
		if (*pout == '}')
			rb++;
		if (rb == 0) {
			if (! isspace (*pout) && *pout != '$') {
				pout--;
				continue;
			}
			else
				break;
		}
		pout--;
		if (lb == rb && lb != 0)
			break;
	}
	pout++;
	if (*pout == '\\') {
		pout++;
		(void) getword (pout, tmp);
		sprintf (ww, "\\%s", tmp);
		pout--;
	}
	else if (*pout == '{')
		(void) get_brace_arg (pout, ww);
	else
		(void) getword (pout, ww);
	*pout = EOS;
	sprintf (tmp,"\\%s %s", w, ww);
	pout = strapp (pout, tmp);
	return (pout);
}



/*
** take care of things like x hat under
*/
char * flip_twice (char *pout, char *w, char *ww)
{
	int lb=0, rb=0;
	char tmp1[MAXWORD], tmp2[MAXWORD];
	
	tmp1[0] = EOS;		tmp2[0] = EOS;
	pout--;
	while (*pout == ' ' || *pout == '\t' || *pout == '\n')
		pout--;
	while (1) {
		if (*pout == '{') {
			lb++;
			if (lb > rb)
				break;
		}
		if (*pout == '}')
			rb++;
		if (rb == 0) {
			if (! isspace (*pout) && *pout != '$') {
				pout--;
				continue;
			}
			else
				break;
		}
		pout--;
		if (lb == rb && lb != 0)
			break;
	}
	pout++;
	if (*pout == '\\') {
		pout++;
		(void) getword(pout,tmp2);
		sprintf(tmp1,"\\%s",tmp2);
		pout--;
	}
	else if (*pout == '{')
		(void) get_brace_arg(pout,tmp1);
	else
		(void) getword(pout,tmp1);
	*pout = EOS;
	sprintf(tmp2,"\\%s{\\%s %s}",w,ww,tmp1);
	pout = strapp(pout,tmp2);
	return(pout);
}



/*
** get argument
** arguments:
**  rec=1 means recursive 
*/
int	get_arg (register char *pin, char *w, int rec)
{
	int c,len,i;
	char ww[MAXWORD];
	int delim;
	
	len=0;
	while ((c = *pin) == ' ' || c == '\t') {	/* skip spaces and tabs */
		pin++;
		len++;
	}
	i=0;
	if (*pin == '{' || *pin == '\"') {
		if (*pin == '{')
			delim = '}';
		else
			delim = '\"';
		pin++;
		len++;
		while ((c = *pin++) != EOS && c != delim && i < MAXWORD) {
			if (c == ' ' && delim == '\"')
				ww[i++] = '\\';
			ww[i++] = (char)c;
			len++;
		}
		len++;
	}
	else {
		while ((c = *pin++) != EOS && !isspace(c)
			   /* && c != '$' && c != '}' */ && i < MAXWORD) {
			if (math_mode && c == '~')
				break;
			ww[i++] = (char)c;
			len++;
		}
	}
	ww[i] = EOS;
	if (rec == 1)				/* check if recursion is required */
		troff_tex(ww,w,1,1);
	else
		strcpy(w,ww);
	return(len);
}


/*
** get all arguments
** arguments:
**  rec=1 means recursive
*/
int	get_allargs (register char *pin, char ***ppw, int rec)
{
	int c, i;
	static char *ww [MAXARGS];
	char w [MAXWORD], *instart;
	int nww;
	int delim;

	instart = pin;
	for (nww = 0; ; nww++) {
		while ((c = *pin) == ' ' || c == '\t')	/* skip spaces and tabs */
			pin++;
		if (c == '\n') {
			pin++;
			ww [nww] = EOS;
			break;
		}
		ww [nww] = pin;
		i=0;
		if (*pin == '{' || *pin == '\"') {
			if (*pin == '{')
				delim = '}';
			else
				delim = '\"';
			ww [nww] = ++pin;
			while ((c = *pin++) != EOS && c != delim && i < MAXWORD)
				/* EMPTY */
				;
			pin [-1] = EOS;
		}
		else {
			while ((c = *pin++) != EOS && !isspace(c)
				   /* && c != '$' && c != '}' */ && i < MAXWORD) {
				if (math_mode && c == '~')
					break;
			}
			pin [-1] = EOS;
			if (c == '\n') {
				ww [nww + 1] = EOS;
				break;
			}
		}
	}
	if (rec == 1) {				/* check if recursion is required */
		for (i = 0; ww [i]; i++) {
			if (ww [i] && *ww [i]) {
				troff_tex (ww [i], w, 1, 1);
				if (strcmp (ww [i], w) != 0)
					ww [i] = strsave (w);
			}
		}
	}
	*ppw = ww;
	return (pin - instart);
}




/*
** get argument surrounded by braces
*/
void get_brace_arg (char *buf, char *w)
{
	int c,i, lb=0, rb=0;
	
	i=0;
	while ((c = *buf++) != EOS) {
		w[i++] = (char)c;
		if (c == '{')	lb++;
		if (c == '}')	rb++;
		if (lb == rb)	break;
	}
	w[i] = EOS;
}

/*
** get "define" or .de word
** arguments:
**  pin    delimited by space only
**  w      delimited by space only
*/
int get_defword (char *pin, char *w, int *illegal)
{
	int c,i;
	
	*illegal = 0;
	for (i=0; (c = *pin++) != EOS && !isspace (c) && i < MAXWORD; i++) {
		w[i] = (char)c;
		if (isalpha(c) == 0)
			*illegal = 1;	/* illegal TeX macro */ 
	}
	w[i] = EOS;
	if (*illegal == 0 && is_forbid(w) >= 0)
		*illegal=1;
	return(i);
}


/*
** get the rest of the line
** arguments:
**  rec=1 means recursion is required
*/
int get_line (char *pin, char *w, int rec)
{
	int c,i,len;
	char ww[MAXLINE];
	
	i=0;
	len=0;
	while ((c = *pin++) != EOS && c != '\n' && len < MAXLINE) {
		ww[i++] = (char)c;
		len++;
	}
	ww[i] = EOS;
	if (rec == 1)
		troff_tex(ww,w,0,1);
	else
		strcpy(w,ww);
	return(len);
}


/*
** get multi-line argument
*/
int get_multi_line (char *pin, char *w)
{
	int len=0,l=0,lines=0;
	char tmp[MAXWORD];
	int c1,c2;
	
	w[0] = EOS;	tmp[0] = EOS;
	while (*pin != EOS) {
		c1 = *pin;
		c2 = *++pin;
		--pin;
		if (c1 == '.' && isupper(c2))
			break; 
		lines++;
		if (lines > 1)
			strcat(w," \\\\\n");
		l = get_line(pin,tmp,1);
		strcat(w,tmp);
		len += l+1;
		pin += l+1;
	}
	len--;
	pin--;
	return(len);
}


/*
** get the macro substitution
*/
int get_mydef (char *pin, char *w)
{
	int c1,c2,l,len;
	char tmp[MAXWORD];
	
	tmp[0] = EOS;
	len=1;
	while (*pin != EOS) {
		c1 = *pin;
		c2 = *++pin;
		--pin;
		if (c1 == '.' && c2 == '.')
			break; 
		l = get_line(pin,tmp,1);
		strcat(w,tmp);
		len += l+1;
		pin += l+1;
	}
	return(len);
}


/*
** get N lines
*/
int get_N_lines (char *pin, char *w, int N)
{
	int len=0,l=0,lines=0;
	char tmp[MAXWORD];
	
	w[0] = EOS;	tmp[0] = EOS;
	while (*pin != EOS && lines < N) {
		lines++;
		if (lines > 1)
			strcat(w," \\\\\n");
		l = get_line(pin,tmp,1);
		strcat(w,tmp);
		len += l+1;
		pin += l+1;
	}
	len--;
	pin--;
	return(len);
}


/*
** get text surrounded by quotes in math mode
*/
int get_no_math (char *pin, char *w)
{
	int c,i,len;
	
	len = 0;
	for (i=0; (c = *pin++) != EOS && c != '\"' && i < MAXWORD; i++) {
		if (c == '{' || c == '}') {
			w[i] = '\\';
			w[++i] = (char)c;
		}
		else
			w[i] = (char)c;
		len++;
	}
	w[i] = EOS;
	return(len);
}


/*
** get the denominator of over
*/
char *get_over_arg (char *pin, char *ww)
{
	char w[MAXWORD], tmp1[MAXWORD], tmp2[MAXWORD];
	int len;
	
	w[0] = EOS;
	tmp1[0] = EOS;
	tmp2[0] = EOS;
	pin += getword (pin,tmp1);		/* read first word */
	pin += skip_white (pin);
	len = getword (pin, tmp2);		/* read second word */
	strcat(w,tmp1);	strcat(w," ");
	
	/* as long as there is a sup or sub read the next two words */
	while (strcmp (tmp2, "sub") == 0 || strcmp (tmp2, "sup") == 0) {
		pin += len;
		strcat (w, tmp2);
		strcat (w, " ");
		pin += skip_white (pin);
		pin += getword (pin, tmp1);
		strcat (w, tmp1);
		strcat (w, " ");
		pin += skip_white (pin);
		len = getword (pin, tmp2);
	}
	troff_tex (w, ww, 0, 1);
	return (pin);
}



/*
** get reference
*/
int get_ref (char *pin, char *w)
{
	int len=0, l=0, lines=0;
	char tmp[MAXWORD];
	
	w[0] = EOS;	tmp[0] = EOS;
	while (*pin != EOS) {
		if (*pin == '\n')
			break;
		(void) getword(pin,tmp);
		if (tmp[0] == '.' && isupper(tmp[1])) {
			/* these commands don't cause a break in reference */
			if (strcmp (tmp, ".R") != 0 && strcmp (tmp, ".I") != 0
				&& strcmp(tmp,".B") != 0)
				break; 
		}
		else if (tmp[0] == '.' && !(isupper(tmp[1]))) {
			/* these commands don't cause a break in reference */
			if (strcmp (tmp, ".br") != 0 && strcmp (tmp, ".bp") != 0)
				break; 
		}
		l = get_line (pin, tmp, 1);
		lines++;
		if (lines > 1)
			strcat (w, " ");
		strcat (w, tmp);
		len += l+1;
		pin += l+1;
	}
	len--;
	pin--;
	return (len);
}


/*
**
*/
void get_size (char *ww, struct measure *PARAMETER)
{
	int sign=0, units=0;
	float value;
	
	if (ww[0] == EOS) {
		if (PARAMETER->def_value == 0) {
			PARAMETER->value = PARAMETER->old_value;
			strcpy(PARAMETER->units,PARAMETER->old_units);
		}
		else {
			PARAMETER->value = PARAMETER->def_value;
			strcpy(PARAMETER->units,PARAMETER->def_units);
		}
	}
	else {
		PARAMETER->old_value = PARAMETER->value;
		strcpy (PARAMETER->old_units, PARAMETER->units);

		parse_units (ww, &sign, &units, &value);

		if (units == 'p')
			strcpy (PARAMETER->units, "pt");
		else if (units == 'i')
			strcpy (PARAMETER->units, "in");
		else if (units == 'c')
			strcpy (PARAMETER->units, "cm");
		else if (units == 'm')
			strcpy (PARAMETER->units, "em");
		else if (units == 'n') {
			value = .5*value;	/* n is about half the width of m */
			strcpy (PARAMETER->units, "em");
		}
		else if (units == 'v')
			strcpy(PARAMETER->units,"ex");
		else if (units == 0) {
			if (sign == 0 || PARAMETER->old_units[0] == EOS)
				strcpy(PARAMETER->units,PARAMETER->def_units);
			else
				strcpy(PARAMETER->units,PARAMETER->old_units);
		}
		else {
			fprintf(stderr,"unknown units %c, using default units\n", units);
			strcpy(PARAMETER->units,PARAMETER->def_units);
		}
		if (sign == 0)
			PARAMETER->value = value;
		else
			PARAMETER->value = PARAMETER->old_value + sign*value;
	}
}



/*
** get the rest of the line -- Nelson Beebe
** arguments:
**  rec=1 means recursion is required
*/
int get_string (char *pin, char *w, int rec)
{
	register int c,i,len;
	char ww[MAXLINE];
	register char *start;
	
	if (*pin != '\"')
		return(get_line(pin,w,rec));
	start = pin;				/* remember start so we can find len */
	i=0;
	pin++;						/* point past initial quote */
	while ((c = *pin++) != EOS && c != '\"' && c != '\n' && i < MAXLINE)
		ww[i++] = (char)c;
	ww[i] = EOS;
	if (c != '\n')				/* flush remainder of line */
		while ((c = *pin++) != '\n')
		/* EMPTY */
			;
	len = pin - start - 1;		/* count only up to NL, not past */
	if (rec == 1)
		troff_tex(ww,w,0,1);
	else
		strcpy(w,ww);
	return(len);
}



/*
** get the argument for sub and sup
*/
int get_sub_arg (char *pin, char *w)
{
	int c,len,i;
	char ww[MAXWORD], tmp[MAXWORD];
	
	len=0;	tmp[0] = EOS;
	while ((c = *pin) == ' ' || c == '\t') {
		pin++;
		len++;
	}
	i=0;
	while ((c = *pin++) != EOS && c != ' ' && c != '\t' && c != '\n'
		   && c != '$' && c != '}' && c != '~' && i < MAXWORD) {
		ww[i++] = (char)c;
		len++;
	}
	ww[i] = EOS;
	if (strcmp(ww,"roman") == 0  || strcmp(ww,"bold") == 0
		|| strcmp(w,"italic") == 0) {
		(void) get_arg(pin,tmp,0);
		sprintf(ww,"%s%c%s",ww,c,tmp);
		len += strlen(tmp)+1;
	}
	troff_tex(ww,w,0,1);		/* recursive */
	return(len);
}


/*
**
*/
int	get_table_entry (char *pin, char *w, int tab)
{
	int c, i=0;
	
	for (i=0; (c = *pin++) != EOS && c != tab && i < MAXWORD; i++)
		w[i] = (char)c;
	w[i] = EOS;
	
	return(i);
}



/*
** get characters till the next space
*/
int get_till_space (char *pin, char *w)
{
	int c,i;
	
	for (i=0; (c = *pin++) != EOS && c != ' ' && c != '\n'
		 && c != '\t' && i < MAXWORD; i++)
		w[i] = (char)c;
	w[i] = EOS;
	return(i);
}



/*
** get the define substitution
*/
int getdef (char *pin, char *ww)
{
	int c,i,len;
	int def_delim;
	char w[MAXWORD];
	
	def_delim = *pin++;		/* take first character as delimiter */
	len=1;
	i=0;
	while ((c = *pin++) != EOS && c != def_delim && i < MAXWORD) {
		len++;
		w[i++] = (char)c;
	}
	w[i] = EOS;
	len++;
	if (c != def_delim) {
		fprintf(stderr,
				"WARNING: missing right delimiter in define, define=%s\n",w);
		len--;
	}
	troff_tex(w,ww,0,1);		/* now translate the substitution */
	return(len);
}



/*
** get an alphanumeric word (dot also)
*/
int getword (char *pin, char *w)
{
	int c,i;
	
	for (i=0; (c = *pin++) != EOS
		 && (isalpha(c) || isdigit(c) || c == '.') && i < MAXWORD; i++)
		w[i] = (char)c;
	if (i == 0 && c != EOS)
		w[i++] = (char)c;
	w[i] = EOS;
	return(i);
}



/*
** change GREEK to Greek
*/
void GR_to_Greek (char *w, char *ww)
{
	*ww++ = '\\';
	*ww++ = *w;
	while(*++w != EOS)
		*ww++ = tolower(*w);
	*ww = EOS;
}



/*
** check if w was defined by the user
*/
int is_def (char *w)
{
	int i;
	
	for (i=0; i < def_count; i++) {
		if (strcmp(def[i].def_macro,w) == 0)
			return(i);
	}
	return(-1);
}



/*
** check if w is in the flip list
*/
int is_flip (char *w)
{
	int i;
	
	for (i=0; i < flip_count; i++) {
		if (strcmp(flip_list[i],w) == 0)
			return(i);
	}
	return(-1);
}



/*
** check if w is one of those sacred macros
*/
int is_forbid (char *w)
{
	int i;
	
	for (i=0; i < forbd_count; i++) {
		if (strcmp(forbid[i],w) == 0)
			return(i);
	}
	return(-1);
}



/*
** check if w has a simple correspondence in TeX
*/
int is_mathcom (char *w, char *ww)
{
	int i;
	
	for (i = 0; i < mathcom_count; i++)
		if (strcmp (math[i].troff_symb, w) == 0) {
			strcpy (ww, math[i].tex_symb);
			return (i);
		}
	return (-1);
}



/*
** check if w is user-defined macro
*/
int is_mydef (char *w)
{
	int i;
	
	for (i=0; i < mydef_count; i++) {
		if (strcmp(mydef[i].def_macro,w) == 0)
			return(i);
	}
	return(-1);
}



/*
** check if w is a macro or plain troff command
*/
int is_troff_mac (char *w, char *ww, int *arg, int *par)
{
	int i;
	
	for (i=0; i < macro_count; i++)
		if (strcmp(macro[i].troff_mac,w) == 0) {
			strcpy(ww,macro[i].tex_mac);
			*arg = macro[i].arg;
			*par = macro[i].macpar;
			return(i);
		}
	return(-1);
}


/*
**
*/
void parse_units (char *ww, int *sign, int *units, float *value)
{
	int len, k=0, i;
	char tmp[MAXWORD];
	
	len = strlen(ww);
	if (ww[0] == '-')
		*sign = -1;
	else if (ww[0] == '+')
		*sign = 1;
	if (*sign != 0)
		k++;
	
	i=0;
	while (k < len) {
		if (isdigit(ww[k]) || ww[k] == '.')
			tmp[i++] = ww[k++];
		else
			break;
	}
	tmp[i] = EOS;
	sscanf(tmp,"%f",value);
	i=0;
	if (k < len) {
		*units = ww[k++];
		if (k < len)
			fprintf(stderr, "Suspect problem in parsing %s, unit used is %c\n",
					ww, *units);
	}
}



/*
** check if w is in the similar list
*/
int similar (char *w)
{
	int i;
	
	for (i=0; i < sim_count ; i++) {
		if (strcmp(sim_list[i],w) == 0)
			return(1);
	}
	return(-1);
}



/*
** ignore the rest of the line
*/
char * skip_line (char *pin)
{
	while (*pin != '\n' && *pin != EOS)
		pin++;
	if (*pin == EOS)
		return(pin);
	else
		return(++pin);
}



/*
** skip white space
*/
int skip_white (char *pin)
{
	int c,len=0;
	
	while ((c = *pin++) == ' ' || c == '\t' || c == '\n')
		len++;	
	return(len);
}



/*
** copy tail[] to s[], return ptr to terminal EOS in s[]
*/
char * strapp (register char *s, register char *tail)
{
	while ((*s++ = *tail++) != 0)
#ifdef DEBUG
		if (debug_o)
			putchar (s[-1]);
#else
		/*EMPTY*/
		;
#endif
	return (s-1);			/* pointer to EOS at end of s[] */
}



/*
** copy input to buffer, buffer holds only MAXLEN characters
*/
void tmpbuf (FILE *in, char *buffer)
{
	int c;
	unsigned int l=0;
	
	while (l++ < MAXLEN && (c = getc(in)) != EOF)
		*buffer++ = (char) c;
	if (l >= MAXLEN) {
		fprintf (stderr,"Sorry: document is too large\n");
		exit(-1);
	}
	*buffer = EOS;
}



/*
** save a string by allocating space
*/
char *strsave (char *s)
{
	char *res;
	
	if ((res = malloc (strlen (s) + 1)) == NULL) {
		fprintf (stderr, "strsave: Can't allocate %d bytes\n", strlen(s) + 1);
		errexit (errno);
	}
	strcpy (res, s);
	return (res);
}

