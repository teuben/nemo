.\" Mcstrip Cyclo
.\" Tools for analysing cyclomatic complexity
.TH "cyclo" 1  "Programmer's tools"
.SH NAME
mcstrip, cyclo \- tools for analysing cyclomatic complexity
.SH SYNOPSIS
.B mcstrip
[-l] [-c] [-r] [-a] [inputfile]
.br
.B cyclo 
[-l] [-t] [-f|-F] [-n func#] [-p] [-s scale] [-c]
.in +6
[-d] [-i]
.in -6
.br
.SH DESCRIPTION
These two tools are used to measure the cyclomatic complexity of
a piece of ANSI C or C++ code.  A postscript flowgraph of the 
functions can also be produced.
.SH COMMANDS
.B mcstrip
.br
This command removes comments, string and character constants from
the input file, or standard input if none is specified.  If any of 
the following options are given, it will give details of the file
instead.
.br

.in +3
[-l] print number of lines in file
.br
[-c] print number of lines that contain comments
.br
[-r] print ratio of lines against lines with comments
.br
[-a] print all of the above
.br
.in -3

.B cyclo
.br
This command takes as its standard input a file that has had the comments,
string and character constants removed (by
.I mcstrip
).  It will then produce either flowcharts
or the cyclomatic complexity number of one or more files in the input.
.br

Information
.in +3
[-l] prints results of tokenisation
.br
[-t] prints results of flow generation
.br
.in -3
Functions
.in +3
[-i] ignores functions declared within a
.br
.in +5
struct/class/union
.in -5
.br 
[-n func#] limits output of -p,-f,-F,-c to the given
.br
.in +11
function number, rather than the default of all functions.
.in -11
.br
[-f] prints function name(s)
.br
[-F] prints functions called by each source function
.br
[-c] print cyclomatic complexity of function(s)
.br
[-p] produce postscript flow graph
.br
[-s scale] scales postscript output (float)
.br
[-d] prints function names to postscript standard
.in +5
output
.in -5
.br
.in -3
.SH RETURN CODES
.B mcstrip
always returns zero.
.B cyclo
returns one in the event of an error, otherwise it returns zero.
.SH EXAMPLES
.B Information about file
.in +3
$ 
.I mcstrip -a main.C
.br
63 comment lines
.br
1199 lines
.br
5 % comments to lines ratio
.br
.in -3
.B Function names
.in +3
$
.I mstrip main.C | cyclo -f
.br
Token::Token
.br
Token::~Token
.br
print_funcs
.br
display_arc
.br
display
.br
.in -3
.B Function calls
.in +3
$
.I mcstrip main.C | cyclo -F -n 3
.br
print_funcs
    setw
    setiosflags
    print_token
    resetiosflags
.in -3
.B Flow charts
.in +3
$
.I mcstrip main.C | cyclo -p -n 3 > flow.ps
.br
$
.I gs flow.ps
.br
.SH MORE DETAILED INFORMATION
The cyclomatic complexity measure counts the number of decisions made within
a piece of code.  This number is equal to the number of linearly independent
paths through the code.  Its main application is in establishing test cases.
Be aware that the measure may not be useful.  It is recommended that further 
reading is done before using it.

The following papers contain further relevant information:
.br

.in +3
\'A Complexity Measure\' McCabe T.J., IEEE Transactions on Software Engineering
vol SE-2, #4 1976, pp 308-320.
.br

\'A Critique of Cyclomatic Complexity as a Software Metric\' Shepperd M., Software
Engineering Journal, March 1988, pp 30-36.
.br

\'Software Defect Prevention Using McCabe\'s Complexity Metric\' Ward W.T., Hewlett-Packard
Journal, April 1989, pp 64-69.
.br

.in -3
.SH SEE ALSO
cflow(1) wc(1)
.SH LIMITATIONS
Only the ANSI C method of declaring functions is recognised.
.br

Constructors that have an initialiser list have their names incorrectly
detected.
.br

Lexical elements that span more than one line are not recognised.  This includes
strings terminating with a continuation slash, struct and the following tag, and
function names and the following opening parenthesis.
.SH COPYRIGHT
(c) 1993 Roger Binns

These tools were produced by Roger Binns for a fourth year project as part of
a computer science degree, for the Computer Science department, Brunel
University, Uxbridge, Middlesex UB8 3PH, United Kingdom.

This software is provided in good faith, having been developed by Brunel
University students as part of their normal course work.  It should not be
assumed that Brunel has any rights of ownership, and the University cannot
accept any liability for its subsequent use.  It is a condition of any such
use that the user idemnifies the University against any claim (including
third party claims) arising therefrom.
