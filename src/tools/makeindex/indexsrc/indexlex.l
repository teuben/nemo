/* $Header */
/* Lex Lexical Analyzer for LaTeX index processor */
/* Roman numeral code written by John Renner (adobe!renner@decwrl.dec.com) */
%{
#undef input
#define input()	(((yytchar=yysptr>yysbuf?U(*--yysptr):getc(InputFile)) \
    ==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
%}
letter		[A-Za-z]
integer		-?[0-9]{1,9}
roman		[mdclxvi]+
whitespace	[ \t\n]
%%
\\indexentry		{
			    LineNo = EndLineNo;
			    return(INDEXENTRY);
			}
{roman}			{			
			    LineNo = EndLineNo;
			    strcpy(yylval.value, yytext);
			    return(ROMAN);
			}
\\{letter}+		{	/* Control word */
			    LineNo = EndLineNo;
			    strcpy(yylval.value, yytext);
			    return(CONTROLSEQUENCE);
			}
\\[^{letter}]		{	/* Control symbol */
			    LineNo = EndLineNo;
			    strcpy(yylval.value, yytext);
			    return(CONTROLSEQUENCE);
			}
\%.*\n			LineNo = EndLineNo++;	/* Comment. Don't pass to parser */
{integer}		{
			    LineNo = EndLineNo;
			    strcpy(yylval.value, yytext);
			    return(INTEGER);
			}
\{ |
\} |
\,			{
			    LineNo = EndLineNo;
			    return(yytext[0]);
			}
{whitespace}+		{
			    char	*p;	/* Utility variable */

			    LineNo = EndLineNo;
			    for (p = yytext; *p; p++)
				if (*p == '\n')
				    EndLineNo++;
			    strcpy(yylval.value, yytext);
			    return(WHITESPACE);
			}
[^\\\{\}\%\, \t\n0-9]+	{
			    LineNo = EndLineNo;
			    strcpy(yylval.value, yytext);
			    return(STRING);
			}
%%
