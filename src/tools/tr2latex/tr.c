/*
** tr2latex - troff to LaTeX converter
** $Id$
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
** 
** Module: tr.c
**
** This module contains the HARD-WIRED rules of the translator.
** It should be handled with care.
*/

#include	"setups.h"
#include	"protos.h"

int def_count = 0;
int mydef_count = 0;

extern bool man;		/* man flag */

void troff_tex (char *pin, char *pout, int mid, int rec)
{
	char eqn_no[MAXWORD], w[MAXWORD], ww[MAXLINE], tmp[MAXWORD], tmp2[MAXWORD];
	char *p, **pw;
	int len,c,c1,c2,i,j;
	int ref = 0;
	int put_brace = 0;
	int first_word = 1;
	int no_word = 1;
	int arg = 0;
	int par = 0;
	int illegal = 0;
	int floating = 0;
	static int delim_defd = 0;	/* whether math delimiter has been defined */
	static char DELIM [] = "$";
	float flen;
	int N;
	int RSRE = 0;			/* block indentation */
	int EXEE = 0;			/* example indentation */
	int NTNE = 0;			/* note indentation */
	int thisfont = 1;		/* default font is roman */
	int lastfont = 1;		/* default last font is roman */
	int offset = 0;			/* amount to offset pin */
	
	*pout = EOS;
	w[0] = EOS;
	ww[0] = EOS;
	tmp[0] = EOS;
	tmp2[0] = EOS;
	while (*pin != EOS) {
		len = getword(pin,w);
		c1 = (first_word? '\n': pin[-1]);
		c2 = pin[0];
		pin += len;
		if (!isspace(w[0]))
			no_word = 0;
		/* first check if we are in math mode */
		if (math_mode) {
			len = get_till_space(pin,ww);
			sprintf(tmp,"%s%s",w,ww);
			if (strcmp(w,"delim") == 0) {
				delim_defd = 1;
				pin += skip_white(pin);
				DELIM[0] = *pin;
				pin = skip_line(pin);
			}
			/* check if it is a math delimiter; switch to non-math mode if so*/
			else if (delim_defd && strcmp(w,DELIM) == 0) {
				math_mode = 0;
				*pout++ = '$';
			}
			/* check for illegal macros here */
			else if (len > 0 && def_count > 0 && (i=is_def(tmp)) >= 0) {
				pin += len;
				pout = strapp(pout,def[i].replace);
			}
			/* See if it is a (legally) defined macro */
			else if (def_count > 0 && (i=is_def(w)) >= 0) {
				if (def[i].illegal)
					pout = strapp(pout,def[i].replace);
				else {
					pout = strapp(pout,"\\");
					pout = strapp(pout,w);
				}
			}
			/* Search for commands in some order;
			   start with non-alphanumeric symbols */
			else if (strcmp(w,"#") == 0 || strcmp(w,"&") == 0
					 || strcmp(w,"%") == 0 || strcmp(w,"_") == 0) {
				pout = strapp(pout,"\\");
				pout = strapp(pout,w);
			}
			else if (strcmp(w,"=") == 0) {
				if (*pin == '=') {
					pin++;
					pout = strapp(pout,"\\equiv");
				}
				else
					pout = strapp(pout,"=");
			}
			else if (strcmp(w,"<") == 0 || strcmp(w,">") == 0) {
				if (*pin == '=') {
					pin++;
					if (strcmp(w,"<") == 0)
						pout = strapp(pout,"\\le");
					else
						pout = strapp(pout,"\\ge");
				}
			}
			else if (strcmp(w,"-") == 0) {
				if (*pin == '>') {
					pin++;
					pout = strapp(pout,"\\to");
				}
				else if (*pin == '+') {
					pin++;
					pout = strapp(pout,"\\mp");
				}
				else
					*pout++ = '-';
			}
			else if (strcmp(w,"+") == 0) {
				if (*pin == '-') {
					pin++;
					pout = strapp(pout,"\\pm");
				}
				else
					*pout++ = '+';
			}
			else if (strcmp(w,"\"") == 0) {
				len = get_no_math(pin,ww);
				pin += len+1;
				if (len > 1) {
					sprintf(tmp,"\\ \\it\\hbox{%s}",ww);
					pout = strapp(pout,tmp);
				}
				else if (len == 1)
					*pout++ = ww[0];
			}
			/* Now search for symbols that start with a captial */
			else if (strcmp(w,".EN") == 0) {
				math_mode = 0;
				if ((len=strlen(eqn_no)) > 0)
				{
					sprintf(tmp,"\\eqno %s",eqn_no);
					pout = strapp(pout,tmp);
				}
				eqn_no[0] = EOS;
				c1 = *--pout;
				c2 = *--pout;
				if (c1 == '\n' && c2 == '$')
					*--pout = EOS;
				else {
					pout += 2;
					pout = strapp(pout,"$$");
				}
			}
			/* Now search for symbols that start with a small letter */
			else if (strcmp(w,"bold") == 0 || strcmp(w,"roman") == 0 ||
					 strcmp(w,"italic") == 0) {
				pin += get_arg(pin,ww,1);
				if (strcmp(w,"bold") == 0) {
					sprintf(tmp,"{\\bf %s}",ww);
					pout = strapp(pout,tmp);
				}
				else if (strcmp(w,"roman") == 0) {
					sprintf(tmp,"{\\rm %s}",ww);
					pout = strapp(pout,tmp);
				}
				else {
					sprintf(tmp,"{\\it %s}",ww);
					pout = strapp(pout,tmp);
				}
			}
			else if (strcmp(w,"define") == 0) {
				if (def_count >= MAXDEF) {
					fprintf(stderr, "Too many defines. MAXDEF=%d\n",MAXDEF);
					exit(-1);
				}
				for (i=0; *--pout != '$' && (unsigned)i < MAXLEN; i++)
					tmp[i] = *pout;
				tmp[i] = EOS;
				strcat(tmp,"$$");
				*--pout = EOS;
				pin += skip_white(pin);
				pin += get_defword(pin,w,&illegal);
				pin += skip_white(pin);
				pin += getdef(pin,ww);
				if (illegal) {
					def[def_count].illegal = 1;
					fprintf(stderr, "illegal TeX macro, %s, replacing it\n",w);
					p = (char *)malloc((unsigned)(strlen(ww)+1)*sizeof(char));
					strcpy(p,ww);
					def[def_count].replace = p;
				}
				else {
					def[def_count].illegal = 0;
					sprintf(tmp2,"\\def\\%s{%s}\n",w,ww);
					pout = strapp(pout,tmp2);
				}
				p = (char *)malloc((unsigned)(strlen(w)+1)*sizeof(char));
				strcpy(p,w);
				def[def_count++].def_macro = p;
				pin += skip_white(pin);
				for (j=i+1; j >= 0; j--)
					*pout++ = tmp[j];
				tmp[0] = EOS;
			}
			else if (strcmp(w,"gsize") == 0 || strcmp(w,"gfont") == 0)
				pin = skip_line(pin);
			else if (strcmp(w,"left") == 0 || strcmp(w,"right") == 0) {
				sprintf(tmp,"\\%s",w);
				pout = strapp(pout,tmp);
				pin += skip_white(pin);
				len = getword(pin,ww);
				if (strcmp(ww,"floor") == 0) {
					pin += len;
					if (strcmp(w,"left") == 0)
						pout = strapp(pout,"\\lfloor");
					else
						pout = strapp(pout,"\\rfloor");
				}
				else if (strcmp(ww,"nothing") == 0 || ww[0] == '\"') {
					pin += len;
					*pout++ = '.';
					if (ww[0] == '\"')	pin++;
				}
				else if (*pin == '{' || *pin == '}')
					*pout++ = '\\';
			}
			else if (strcmp(w,"over") == 0) {
				if (!first_word) {
					pout--;
					for (i=0; isspace (*pout); i++)
						tmp[i] = *pout--;
					if (*pout == '}' && put_brace == 0)
						*pout = ' ';
					else {
						for (; !isspace(*pout) && *pout != '$'; i++)
							tmp[i] = *pout--;
						put_brace = 0;
						*++pout = '{';
					}
					for (j=i-1; j >= 0; j--)
						*++pout = tmp[j];
					*++pout = EOS;
				}
				pout = strapp(pout,"\\over");
				pin += skip_white(pin);
				*pout++ = ' ';
				if (*pin == '{')
					pin++;
				else {
					pin = get_over_arg(pin,ww);
					pout = strapp(pout,ww);
					if (*pin != EOS || !first_word)
						*pout++ = '}';
				}
			}
			else if (strcmp(w,"size") == 0)
				pin += get_arg(pin,ww,0);
			else if (strcmp(w,"sup") == 0 || strcmp(w,"to") == 0 ||
					 strcmp(w,"sub") == 0 || strcmp(w,"from") == 0) {
				while ((c = *--pout) == ' ' || c == '\t' || c == '\n')
					/* EMPTY */
					;
				*++pout = EOS;
				if (strcmp(w,"sup") == 0 || strcmp(w,"to") == 0)
					pout = strapp(pout,"^");
				else
					pout = strapp(pout,"_");
				pin += skip_white(pin);
				len = get_sub_arg(pin,ww);
				pin += len;
				if (len > 1) {
					sprintf(tmp,"{%s}",ww);
					pout = strapp(pout,tmp);
					len = skip_white(pin);
					pin += len;
					(void) getword(pin,ww);
					if (strcmp(ww,"over") == 0)
						put_brace = 1;
					pin -= len;
				}
				else
					pout = strapp(pout,ww);
			}
			else if (strcmp(w,"up") == 0 || strcmp(w,"down") == 0
					 || strcmp(w,"fwd") == 0 || strcmp(w,"back") == 0) {
				if (strcmp(w,"up") == 0) {
					pout = strapp(pout,"\\raise");
					strcpy(tmp,"ex");
				}
				else if (strcmp(w,"down") == 0) {
					pout = strapp(pout,"\\lower");
					strcpy(tmp,"ex");
				}
				else if (strcmp(w,"fwd") == 0) {
					pout = strapp(pout,"\\kern");
					strcpy(tmp,"em");
				}
				else if (strcmp(w,"back") == 0) {
					pout = strapp(pout,"\\kern-");
					strcpy(tmp,"em");
				}
				pin += skip_white(pin);
				pin += getword(pin,ww);
				len = atoi(ww);		flen = len/100.;
				ww[0] = EOS;
				sprintf(tmp2,"%4.2f%s",flen,tmp);
				pout = strapp(pout,tmp2);
			}
			/* Now check if the word is a member of a group */
			else if (CAP_GREEK(w) > 0) {
				GR_to_Greek(w,ww);
				pout = strapp(pout,ww);
			}
			else if (is_flip(w) >= 0) {
				if (!first_word) {
					len = skip_white(pin);
					pin += len;
					(void) getword(pin,ww);
					if (is_flip(ww) >= 0) {
						pin += strlen(ww);
						pout = flip_twice(pout,w,ww);
					}
					else {
						pin -= len;
						pout = flip(pout,w);
					}
				}
				else {
					pout = strapp(pout,"\\");
					pout = strapp(pout,w);
				}
			}
			else if (is_mathcom (w, ww) >= 0)
				pout = strapp (pout, ww);
			else if (similar(w) > 0) {
				pout = strapp (pout, "\\");
				pout = strapp (pout, w);
			}
			
			/* if none of the above math commands matched, it is an ordinary
			   symbol; just copy it */
			
			else
				pout = strapp(pout,w);
		}
		
		/* check if it is a math delimiter; switch to math mode if so */
		
		else if (strcmp(w,"$") == 0 && de_arg > 0) {
			de_arg++;
			*pout++ = '#';
		}
		else if (delim_defd && strcmp(w,DELIM) == 0) {
			math_mode = 1;
			*pout++ = '$';
		}
		else if (strcmp(w,"$") == 0)
			pout = strapp(pout,"\\$");
		
		/* check if it is a non-math troff command */
		
		else if (c2=='.' && !mid && (c1=='\n' || c1==EOS || first_word)) {
			/* Search in some order; start with non-alphanumeric characters */
			if (strcmp(w,".") == 0) {
				c1 = *pin;
				c2 = *++pin;
				if (c1 == '\\' && c2 == '\"') {
					++pin;
					pin += get_line(pin,ww,0);
					pout = strapp(pout,"%");
					pout = strapp(pout,ww);
				}
				else {
					fprintf(stderr,
							"I cannot translate troff macro .%c%c\n",c1,c2);
					pin += get_line(pin,ww,0);
					sprintf(tmp,"%%.%c%c",c1,c2);
					pout = strapp(pout,tmp);
					pout = strapp(pout,ww);
					if (*pin == EOS)	*pout++ = '\n';
				}
			}
			/* Now search for commads that start with a capital */
			else if (strcmp(w,".AB") == 0) {
				pin += get_arg(pin,ww,0);
				if (strcmp(ww,"no") == 0)
					pout = strapp(pout,"\\bigskip");
				else
					pout = strapp(pout,"\\begin{abstract}");
			}
			else if (strcmp(w,".B") == 0 || strcmp(w,".bf") == 0 ||
					 strcmp(w,".I") == 0 || strcmp(w,".it") == 0 ||
					 strcmp(w,".R") == 0 || strcmp(w,".rm") == 0 ||
					 strcmp(w,".P") == 0) {
				if (strcmp(w,".R") == 0 || strcmp(w,".rm") == 0)
					strcpy(w,"rm");
				else if (strcmp(w,".B") == 0 || strcmp(w,".bf") == 0)
					strcpy(w,"bf");
				else if (strcmp(w,".I") == 0 || strcmp(w,".it") == 0)
					strcpy(w,"it");
				else {
					switch(lastfont) {
					case 1:
						strcpy (w, "rm");
						thisfont = 1;
						break;
					case 2:
						strcpy (w, "it");
						thisfont = 2;
						break;
					case 3:
						strcpy (w, "bf");
						thisfont = 3;
						break;
					default:
						strcpy (w, "rm");
						thisfont = 1;
						break;
					}
				}
				pin += get_arg(pin,ww,1);
				if (ww[0] == EOS) {
					pout = strapp(pout,"\\");
					pout = strapp(pout,w);
				}
				else {
					sprintf(tmp,"{\\%s %s",w,ww);
					pout = strapp(pout,tmp);
					while (1) {
						pin += get_arg(pin,ww,1);
						if (ww[0] == EOS)
							break;
						sprintf(tmp," %s",ww);
						pout = strapp(pout,tmp);
					}
					pout = strapp(pout,"}");
				}
			}
			else if (man && (strcmp(w,".BR") == 0 || strcmp(w,".BI") == 0
							 || strcmp(w,".IR") == 0 || strcmp(w,".IB") == 0
							 || strcmp(w,".RI") == 0 || strcmp(w,".RB") == 0)){
				pout = alternate(pin,pout,w);
				pin = skip_line(pin);
				*pout++ = '\n';
			}
			else if (strcmp(w,".BX") == 0) {
				pin += get_arg(pin,ww,1);
				sprintf(tmp,"\\fbox{%s}",ww);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".EQ") == 0) {
				math_mode = 1;
				put_brace = 0;
				pout = strapp(pout,"$$");
				len = get_arg(pin,eqn_no,0);
				if (strcmp(eqn_no,"I") == 0 || strcmp(eqn_no,"L") == 0) {
					fprintf(stderr,"lineups are ignored\n");
					pin += len;
					len = get_arg(pin,eqn_no,0);
				}
				if ((strlen(eqn_no)) > 0)
					pin += len;
				len = get_arg(pin,tmp,0);
				if (strcmp(tmp,"I") == 0 || strcmp(tmp,"L") == 0) {
					fprintf(stderr,"lineups are ignored\n");
					pin += len;
				}
			}
			else if (strcmp (w,".IP") == 0) {
				pin += get_arg (pin, ww, 1);
				pin = skip_line (pin);
				if (IP_stat == 0)
					pout = strapp(pout,"\\begin{IPlist}\n");
				sprintf(tmp,"\\IPitem{{%s}}\n",ww);
				pout = strapp(pout,tmp);
				if (de_arg > 0)		mydef[mydef_count].par = 2;
				else			IP_stat = 1;
			}
			else if (strcmp(w,".KE") == 0) {
				if (floating)
					pout = strapp(pout,"\\end{figure}");
				else
					pout = strapp(pout,"}");
				floating = 0;
			}
			else if (strcmp(w,".KF") == 0) {
				floating = 1;
				pout = strapp(pout,"\\begin{figure}");
			}
			else if (strcmp(w,".QP") == 0) {
				if (de_arg > 0)
					mydef[mydef_count].par = 4;
				else
					QP_stat = 1;
				pout = strapp(pout,"\\begin{quotation}");
			}
			else if (strcmp(w,".RE") == 0) {
				if (RSRE == 0)
					fprintf(stderr,".RE with no matching .RS\n");
				else {
					sprintf(tmp,"\\ind{%d\\parindent}", RSRE--);
					pout = strapp(pout,tmp);
				}
			}
			else if (strcmp(w,".RS") == 0) {
				pin += get_arg (pin, ww, 1);
				pin = skip_line (pin);
				sprintf(tmp,"\\ind{%d\\parindent}", ++RSRE);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".Re") == 0) {
				if (ref == 0)
					pout = strapp(pout,"\\REF\n");
				ref++;
				pin = skip_line(pin);
				pin += get_ref(pin,ww);
				sprintf(tmp,"\\reference{%s}",ww);
				pout = strapp(pout,tmp);
			}
			else if (man && (strcmp(w,".TP") == 0 || strcmp(w,".HP") == 0)) {
				if (IP_stat && TP_stat) {
					pout = strapp(pout,"\\end{IPlist}%\n");
					IP_stat = 0;
				}
				if (QP_stat && TP_stat) {
					pout = strapp(pout,"\\end{quotation}%\n");
					QP_stat = 0;
				}
				pin = skip_line(pin);
				pin += get_line(pin,ww,1);
				if (TP_stat == 0) {
					sprintf(tmp,"\\begin{TPlist}{%s}\n",ww);
					pout = strapp(pout,tmp);
				}
				sprintf(tmp,"\\item[{%s}]",ww);
				pout = strapp(pout,tmp);
				if (de_arg > 0)
					mydef[mydef_count].par = 3;
				else
					TP_stat = 1;
			}
			else if (man && (strcmp(w,".TH") == 0)) {
				/* expect something like .TH LS 1 "September 4, 1985"*/
				pin += get_allargs (pin, &pw, 1);
				if (pw [0] == NULL || pw [1] == NULL)
					fprintf (stderr, "Missing argument in troff macro .TH\n");
				else {
					for (j = 2; j <= 5; j++)
						if (pw [j] == NULL)
							break;
					for (; j <= 5; j++)
						pw[j] = "";
					sprintf (tmp, "\\phead{%s}{%s}{%s}{%s}{%s}",
							 pw[0], pw[1], pw[2], pw[3], pw[4]);
					pout = strapp (pout, tmp);
				}
			}
			else if (man && strcmp (w, ".PN") == 0) {
				pin += get_allargs (pin, &pw, 1);
				if (pw [0] == NULL)
					fprintf (stderr, "Missing argument in troff macro .PN\n");
				else {
					sprintf (tmp, "{\\tt{}%s}%s ",
							 pw [0], pw [1]? pw [1]: "");
					pout = strapp (pout, tmp);
				}
			}
			else if (man && strcmp (w, ".CW") == 0)
				pout = strapp (pout, "\\tt{}");
			else if (man && strcmp (w, ".MS") == 0) {
				pin += get_allargs (pin, &pw, 1);
				if (pw [0] == NULL || pw [1] == NULL)
					fprintf (stderr, "Missing arguments in troff macro .MS\n");
				else {
					sprintf (tmp, "{\\tt{}%s}(%s)%s ",
							 pw [0], pw [1], pw [2]? pw [2]: "");
					pout = strapp (pout, tmp);
				}
			}
			else if (man && strcmp (w, ".EX") == 0) {
				pin += get_arg (pin, ww, 1);
				pin = skip_line (pin);
				sprintf(tmp, "\\nofill\\ind{%d\\parindent}\\tt{}", ++EXEE);
				pout = strapp(pout,tmp);
			}
			else if (man && strcmp (w, ".EE") == 0) {
				if (EXEE == 0)
					fprintf(stderr,".EE with no matching .EX\n");
				else {
					--EXEE;
					pout = strapp(pout, "\\fill ");
				}
			}
			else if (man && strcmp (w, ".NT") == 0) {
				++NTNE;
				pin += get_allargs (pin, &pw, 1);
				if (pw [0] == NULL)
					pout = strapp (pout, "\\beginnotec{NOTE}");
				else {
					if (strcmp (pw [0], "C") == 0) {
						if (pw [1])
							pout = strapp (pout, "\\beginnotec{NOTE}");
						else {
							sprintf (tmp, "\\beginnotec{%s}", pw [1]);
							pout = strapp (pout, tmp);
						}
					}
					else {
						sprintf (tmp, "\\beginnote{%s}", pw [0]);
						pout = strapp (pout, tmp);
					}
				}
			}
			else if (man && strcmp (w, ".NE") == 0) {
				if (NTNE == 0)
					fprintf(stderr,".NE with no matching .NT\n");
				else {
					--NTNE;
					pout = strapp(pout, "\n\\end{quote}\n");
				}
			}
			else if (strcmp(w,".TS") == 0)
			{
				fprintf(stderr,"I am not very good at tables\n\
I can only do very simple ones. You may need to check what I've done\n");
				pin = skip_line(pin);
				pout = do_table(pin,pout,&offset);
				pin += offset;
				offset = 0;		/* reset */
			}
			/* Now search for commands that start with small letters */
			else if (strcmp(w,".TE") == 0) {
				fprintf(stderr,"Oops! I goofed. I told you I am not very good at tables.\nI have encountered a table end but I am not in table mode\n");
			}
			else if (strcmp(w,".de") == 0) {
				de_arg = 1;
				if (mydef_count >= MAXDEF) {
					fprintf(stderr,
							"Too many .de's. MAXDEF=%d\n",MAXDEF);
					exit(-1);
				}
				pin += skip_white(pin);
				pin += get_defword(pin,w,&illegal);
				pin += skip_white(pin);
				pin += get_mydef(pin,ww);
				mydef[mydef_count].arg_no = de_arg;
				if (illegal) {
					mydef[mydef_count].illegal = 1;
					fprintf(stderr,
							"illegal TeX macro, %s, replacing it\n",w);
					p = (char *)malloc((unsigned)(strlen(ww)+2)*
									   sizeof(char));
					sprintf(p,"%s",ww);
					mydef[mydef_count].replace = p;
				}
				else {
					mydef[mydef_count].illegal = 0;
					sprintf(tmp,"\\def\\%s",w);
					pout = strapp(pout,tmp);
					for (j=1; j<de_arg; j++) {
						sprintf(tmp,"#%d",j);
						pout = strapp(pout,tmp);
					}
					sprintf(tmp,"{%s}\n",ww);
					pout = strapp(pout,tmp);
				}
				p = (char *)malloc((unsigned)(strlen(w)+2)*sizeof(char));
				sprintf(p,".%s",w);
				mydef[mydef_count++].def_macro = p;
				pin = skip_line(pin);
				de_arg = 0;
			}
			else if (strcmp(w,".ds") == 0) {
				pin += get_arg(pin,w,0);
				pin += skip_white(pin);
				pin += get_line(pin,ww,1);
				if (strcmp(w,"LH") == 0) {
					sprintf(tmp,"\\lefthead{%s}",ww);
					pout = strapp(pout,tmp);
				}
				else if (strcmp(w,"RH") == 0) {
					sprintf(tmp,"\\righthead{%s}",ww);
					pout = strapp(pout,tmp);
				}
				else if (strcmp(w,"CF") == 0) {
					if (index(ww,'%') == 0) {
						sprintf(tmp,"\\footer{%s}",ww);
						pout = strapp(pout,tmp);
					}
					else
						pout = strapp(pout, "\\footer{\\rm\\thepage}");
				}
				else {
					fprintf(stderr,"I do not understand .ds %s\n",w);
					sprintf(tmp,"%%.ds %s %s",w,ww);
					pout = strapp(pout,tmp);
				}
			}
			else if (strcmp(w,".sp") == 0) {
				pin += get_arg(pin,ww,0);
				(void) get_size(ww,&space);
				sprintf(tmp,"\\par\\vspace{%3.1f%s}",space.value,space.units);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".in") == 0) {
				pin += get_arg(pin,ww,0);
				(void) get_size(ww,&indent);
				sprintf(tmp,"\\ind{%3.1f%s}",indent.value,indent.units);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".ls") == 0) {
				pin += get_arg(pin,ww,0);
				(void) get_size(ww,&linespacing);
				sprintf(tmp,"\\baselineskip=%3.1f%s",linespacing.value,
						linespacing.units);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".so") == 0) {
				pin += get_arg(pin,ww,0);
				sprintf(tmp,"\\input %s",ww);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".ti") == 0) {
				pin += get_arg(pin,ww,0);
				tmpind.value = indent.value;
				strcpy(tmpind.units,indent.units);
				(void) get_size(ww,&tmpind);
				sprintf(tmp,"\\tmpind{%3.1f%s}",
						tmpind.value,tmpind.units);
				pout = strapp(pout,tmp);
			}
			else if (strcmp(w,".vs") == 0) {
				pin += get_arg(pin,ww,0);
				(void) get_size(ww,&vspace);
				sprintf(tmp,"\\par\\vspace{%3.1f%s}",
						vspace.value,vspace.units);
				pout = strapp(pout,tmp);
			}
			/* check if it is a member of a group */
			else if (mydef_count > 0 && (i=is_mydef(w)) >= 0) {
				if (mydef[i].par > 0) {
					if (de_arg > 0)
						mydef[mydef_count].par = mydef[i].par;
					else {
						pout = end_env(pout);
						pout = strapp(pout,"\n");
					}
				}
				if (mydef[i].illegal)
					pout = strapp(pout,mydef[i].replace);
				else {
					w[0] = '\\';	/* replace dot by backslash */
					pout = strapp(pout,w);
				}
				for (j=1; j <mydef[i].arg_no; j++) {
					pin += get_arg(pin,ww,1);
					sprintf(tmp,"{%s}",ww);
					pout = strapp(pout,tmp);
				}
				if (de_arg == 0)
					envoke_stat(mydef[i].par);
			}
			else if ((i=is_troff_mac(w,ww,&arg,&par)) >= 0) {
				if (par > 0) {
					if (de_arg > 0)
						mydef[mydef_count].par = par;
					else
						pout = end_env(pout);
				}
				pout = strapp(pout,ww);
				if (ww[0] == EOS)
					pin = skip_line(pin);
				if (ww[0] != EOS && arg == 0) {
					pin = skip_line(pin);
					*pout++ = '\n';
				}
				if (arg > 0) {
					if (arg == 1) {
						pin += skip_white(pin);
						pin += get_string(pin,ww,1);
					}
					else {
						if (isupper(w[1])) {
							pin = skip_line(pin);
							pin += get_multi_line(pin,ww);
						}
						else {
							pin += get_arg(pin,tmp,0);
							pin = skip_line(pin);
							if (tmp[0] == EOS)
								N = 1;
							else
								N = atoi(tmp);
							pin += get_N_lines(pin,ww,N);
						}
					}
					sprintf(tmp2,"{%s}",ww);
					pout = strapp(pout,tmp2);
				}
			}
			/* if none of the above commands matched, it is either
			   an illegal macro or an unknown command */
			else {
				len = get_till_space(pin,ww);
				sprintf(tmp,"%s%s",w,ww);
				if (mydef_count > 0 && (i=is_mydef(tmp)) >= 0) {
					pin += len;
					if (mydef[i].par > 0) {
						if (de_arg > 0)
							mydef[mydef_count].par=mydef[i].par;
						else {
							pout = end_env(pout);
							pout = strapp(pout,"\n");
						}
					}
					pout = strapp(pout,mydef[i].replace);
					for (j=1; j <mydef[i].arg_no; j++) {
						pin += get_arg(pin,ww,1);
						sprintf(tmp,"{%s}",ww);
						pout = strapp(pout,tmp);
					}
					if (de_arg == 0)
						envoke_stat(mydef[i].par);
				}
				else {
					fprintf(stderr, "I cannot translate troff macro %s\n",w);
					pin += get_line(pin,ww,0);
					pout = strapp(pout,"%");
					pout = strapp(pout,w);
					pout = strapp(pout,ww);
					if (*pin == EOS)
						*pout++ = '\n';
				}
			}
		}
		
		/* some manuals have commented lines beginning with ''' */
		else if (c2=='\'' && !mid && (c1=='\n' || c1==EOS || first_word)) {
			if (*pin == '\'') {
				pin++;
				if (*pin == '\'') {
					pin++;
					pout = strapp(pout,"%");
				}
				else
					pout = strapp(pout,"''");
			}
			else
				pout = strapp(pout,"'");
		}
		
		/* See if it is one of these symbols */
		
		else if (strcmp(w,"#") == 0 || strcmp(w,"&") == 0 ||
				 strcmp(w,"{") == 0 || strcmp(w,"}") == 0 ||
				 strcmp(w,"%") == 0 || strcmp(w,"_") == 0 ||
				 strcmp(w,"~") == 0 || strcmp(w,"^") == 0 ) {
			pout = strapp(pout,"\\");
			pout = strapp(pout,w);
			if (strcmp(w,"~") == 0 || strcmp(w,"^") == 0)
				pout = strapp(pout,"{}");
		}
		
		else if (strcmp(w,">") == 0 || strcmp(w,"<") == 0
				 || strcmp(w,"|") == 0) {
			sprintf(tmp,"$%s$",w);
			pout = strapp(pout,tmp);
		}
		
		/* check for backslash commands */
		
		else if (strcmp(w,"\\") == 0) {
			switch (*pin) {
			case ' ':
			case '\t':
			case '\n':
				pout = strapp(pout,"\\");
				*pout++ = *pin++;
				break;
			case EOS:
				break;
			case '-':
				pin++;
				pout = strapp(pout,"--");
				break;
			case '~':
			case '^':
				pin++;
				pout = strapp(pout,"\\/");
				break;
			case '0':
				pin++;
				pout = strapp(pout,"\\ ");
				break;
			case 'e':
				pin++;
				pout = strapp(pout,"\\bs ");
				break;
			case '\\':
				pin++;
				if (*pin == '$' && de_arg > 0) {
					pin++;
					de_arg++;
					*pout++ = '#';
				}
				else
					pout = strapp(pout,"\\bs ");
				break;
			case '`':
			case '\'':
				break;				/* do nothing */
			case '"':
				pin++;
				pin += get_line(pin,ww,0);
				pout = strapp(pout,"%");
				pout = strapp(pout,ww);
				break;
			case '|':
				pin++;
				pout = strapp(pout,"\\,");
				break;
			case '&':
				pin++;
				break;
			case '(':
				c1 = *++pin;
				c2 = *++pin;
				pin++;
				if (c1 == 'e' && c2 == 'm')
					pout = strapp(pout,"---");
				else if (c1 == 'd' && c2 == 'e')
					pout = strapp(pout,"$^\\circ$");
				else
					fprintf(stderr,
							 "I am not prepared to handle \\(%c%c\n",c1,c2);
				break;
			case 's':
				pin +=3;
				break;
			case '*':
				c1 = *++pin;
				pin++;
				switch (c1) {
				case ':':	pout = strapp(pout,"\\\""); break;
				case 'C':	pout = strapp(pout,"\\v"); break;
				case ',':	pout = strapp(pout,"\\c"); break;
				case '(':
					sprintf(tmp,"\\%c",c1);
					pout = strapp(pout,tmp);
					break;
				default:
					fprintf(stderr,"I am not prepared to handle \\*( cases\n");
					pin += 2;
					break;
				}
				if (c1 != '(') {
					c1 = *pin++;
					sprintf(tmp,"{%c}",c1);
					pout = strapp(pout,tmp);
				}
				break;
			case 'f':
				c1 = *++pin;
				pin++;
				switch (c1) {
				case '1':
				case 'R':
					lastfont = thisfont;
					thisfont = 1;
					if (isspace (*pin)) {
						*pout++ = ' ';
						pin++;
					}
					pout = strapp(pout,"%\n\\rm ");
					break;
				case '2':
				case 'I':
					lastfont = thisfont;
					thisfont = 2;
					if (isspace (*pin)) {
						*pout++ = ' ';
						pin++;
					}
					pout = strapp(pout,"%\n\\it ");
					break;
				case '3':
				case 'B':
					lastfont = thisfont;
					thisfont = 3;
					if (*pin == ' ' || *pin == '\t' ||
						*pin == '\n' || *pin == '\f') {
						*pout++ = ' ';	pin++;
					}
					pout = strapp(pout,"%\n\\bf ");
					break;
				case 'P':
					/* preserve white space - from Nelson Beebe  */
					if (isspace (*pin)) {
						*pout++ = ' ';
						pin++;
					}
					switch(lastfont) {
					case 1:
						pout = strapp(pout,"\\rm%\n");
						thisfont = 1;
						break;
					case 2:
						pout = strapp(pout,"\\it%\n");
						thisfont = 2;
						break;
					case 3:
						pout = strapp(pout,"\\bf%\n");
						thisfont = 3;
						break;
					default:
						pout = strapp(pout,"\\rm%\n");
						thisfont = 1;
						break;
					}
					break;
				default:
					fprintf(stderr, "I do not understand \\f%c yet\n",c1);
					break;
				}
				break;
			default:
				fprintf(stderr,"I am not prepared to handle \\%c\n",*pin++);
				break;
			}
		}
		
		/* if non of the above checks, its a dull word; copy it */
		
		else
			pout = strapp(pout,w);
		*pout = EOS;
		ww[0] = EOS;
		tmp[0] = EOS;
		tmp2[0] = EOS;
		if (!no_word)
			first_word = 0;
	}
	/* if file end, close opened environments and delimitters */
	if (rec == 0) {
		if (IP_stat)	pout = strapp(pout,"\\end{IPlist}\n");
		if (QP_stat)	pout = strapp(pout,"\\end{quotation}\n");
		if (TP_stat)	pout = strapp(pout,"\\end{TPlist}\n");
	}
	
	*pout = EOS;
}
