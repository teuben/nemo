.TH UPDATE 1NEMO "22 May 1987" 
.SH NAME
update \- combine multiple fortran source files for CYBER 205's UPDATE
.SH SYNOPSIS
.PP
\fBupdate [hfile1  hfile2 ... hfilen] [\fB-d\fP file1 file2 ... filen]
.SH DESCRIPTION
When using multiple source files, and using the standard 'include' facility
that  many compilers offer (unix, VMS), this utility converts these files
into one source deck that the Cyber UPDATE facility understands.
.PP
The outputfile name is the filename of the last file (\fIfilen\fP) with
an extension ".UPD". For insanity reasons the output filename is in
upper case.
.PP
This program has only been tested with fortran source files, since
it thinks comment lines start with a '\fBC\fP'. The C-definition of comments
(\fB/*\fP through \fB*/\fP) is not implemented here. As long as the word
\fBinclude\fP does not appear in those comments, a C-source deck should
also work.
.SH OPTIONS
.TP 20
\fB-c\fP
leave out the comments [default: they are left in].
.TP
\fB-a\fIext\fP
the default extension of the input filenames,
in case all filenames have the same extension
[default: none, user must supply them]
.SH BUGS
See remarks on C-source files under DESCRIPTION.
.SH "SEE ALSO"
\fIusage(0)\fP
.SH FILES
.nf
.ta +2.0i
anchor:/usr/pjt/205/doc/update*   	sources and executables
.fi
.SH AUTHOR
Peter Teuben
.SH "UPDATE HISTORY"
.nf
.ta +2.0i +2.0i
31-Mar-87	created    	PJT
22-May-87	include files bug	PJT
.fi
