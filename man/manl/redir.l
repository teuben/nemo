.TH REDIR 1 "5 February 1989"
.SH NAME
redir - redirect standard output and/or standard error (diagnostic output)
.SH SYNOPSIS
redir [-[a][oeOE] file] ... command arg1 arg2 ...
.SH DESCRIPTION
Redir executes the named command with the given arguments,
redirecting the standard output and/or standard input according
to the options given.  This program exists to remedy a serious
deficiency in csh(1), namely that csh does not allow stdout
and stderr to be redirected separately.  The options are:
.TP 20
.B -o file        
Redirect the standard output of the named command
to the file "file", only if "file" does not yet exist.
.TP
.B -O file        
Redirect the standard output of the named command
to the file "file", clobbering the existing
version of "file" if it exists.
.TP
.B -e     
Redirect the standard error (diagnostic output) of the
named command to file "file", only if "file" does not yet exist.
.TP
.B -E     
Redirect the standard error (diagnostic output) of the
named command to the file "file", clobbering the existing
version of "file" if it exists.
.TP
.B -ao file       
Append the standard output of the named command
to the file "file", whether or not "file" exists.
.TP
.B -ae file       
Append the standard error (diagnostic output) of the
named command to the file "file", whether or not "file" exists.


-aO and -aE are synonyms for -ao and -ae.
.SH EXAMPLE
The following example demonstrates a means for piping the output
of sed through ctags and then into awk; errors from ctags will be
placed in the file "errors":
.nf

% ... | sed -f sed-script | redir -e errors ctags -x | awk -f awk-script | ...

.fi
When you want to pipe command A through command B,
collecting the stderr of A in A.err and the stderr of B in B.err, while
displaying the stdout of B, you could do the following:
.nf

%          redir -e A.err A | redir -e B.err B
.fi
.SH AUTHOR
Bob Glickstein
.nf
       Information Technology Center
       Carnegie Mellon University
       Pittsburgh, PA
       22-Jan-89
       ARPAnet: bobg@andrew.cmu.edu
.fi
.SH UDPATES
.nf
.ta +1i +4i
5-feb-88	argc==1 check added	PJT
.fi
