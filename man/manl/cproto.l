.\" $Id$
.\"
.de EX          \"Begin example
.ne 5
.if n .sp 1
.if t .sp .5
.nf
.in +.5i
..
.de EE		\"End example
.fi
.in -.5i
.if n .sp 1
.if t .sp .5
..
.TH CPROTO 1 "February 28, 1991"
.SH NAME
cproto \- generate C function prototypes from C source code
.SH SYNOPSIS
.B cproto
[ 
.I option \fP...\fI
] [
.I file \fP...\fI
]
.SH DESCRIPTION
.B Cproto
reads C source code files and outputs function prototypes for external
functions defined in the source to standard output.
The function definitions may be in the old style or ANSI style.
Optionally,
.B cproto
also outputs declarations for any external variables defined in the file.
If no
.I file
argument is given,
.B cproto
takes its input from the standard input.
.SH OPTIONS
.TP
.B \-e
Output the keyword
.B extern
in front of each declaration having global scope.
.TP
.BI \-f n
Set the style of function prototype where
.I n
is a number from 0 to 4.
For example, consider the function definition
.EX
main (argc, argv)
int argc;
char *argv[];
{
 ...
}
.EE
If the value is 0, then no prototypes are generated.
When set to 1, the output is:
.EX
int main(/*int argc, char *argv[]*/);
.EE
For a value of 2, the output has the form:
.EX
int main(int /*argc*/, char */*argv*/[]);
.EE
The default value is 3.
It produces the full function prototype:
.EX
int main(int argc, char *argv[]);
.EE
A value of 4 produces prototypes guarded by a macro:
.EX
int main P_((int argc, char *argv[]));
.EE
.TP
.BI \-m name
Set the name of the macro used to guard prototypes when option -f4 is selected.
By default it is "P_".
.TP
.B \-d
Omit the definition of the prototype macro named by the -m option.
.TP
.B \-p
Disable promotion of formal parameters in function prototypes.
By default, parameters of type
.B char
or
.B short
in traditional style function definitions are promoted to type
.B int
in the function prototype.
Parameters of type
.B float
get promoted to 
.B double
as well.
.TP
.B \-s
Also output
.B static
declarations.
.TP
.B \-v
Also output declarations for variables defined in the file.
.TP
.BI \-F string
Set the format used to output each prototype.
The string is a template in the form
.EX
" int main ( a, b )"
.EE
where each space in the string may be replaced with whitespace characters.
For example, the option
.EX
-F"int main(\\n\\ta,\\n\\tb\\n\\t)"
.EE
will produce prototypes in the format
.EX
int main(
        int argc,
        char *argv[]
        );
.EE
.TP
.BI \-D name\[=value\]
This option is passed through to the preprocessor and is used to define 
symbols for use with conditionals such as
.I #ifdef.
.TP
.BI \-U name
This option is passed through to the preprocessor and is used to remove
any definitions of this symbol.
.TP
.BI \-I directory
This option is passed through to the preprocessor and is used to specify
a directory to search for files that are referenced with 
.I #include.
.TP
.B \-V
Print version information.
.SH AUTHOR
Chin Huang
cthuang@contact.uucp
chin.huang@canrem.uucp
.SH BUGS
When cproto finds an error, it usually outputs the not very descriptive
message "syntax error".
.TP
Options that take string arguments only interpret the following
character escape sequences:
.EX
\\n	newline
\\t	tab
.EE
.SH "SEE ALSO"
cc(1),
cpp(1)
.SH FILES
NEMO/usr/lib/cproto - source code location within NEMO
