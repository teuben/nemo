.TH TR2TEX 1 "14 June 1991"
.UC 4
.SH NAME
tr2latex \- convert a document from troff to LaTeX
.SH SYNOPSIS
.B tr2latex
[
.B \-m
] [
.B \-t
] [
.BI \- n
] [
.B \-s
.I style
] [
.B \-o
.I outfile
]
.IR filename ...
.SH DESCRIPTION
.B Tr2latex
converts a document typeset in
.B troff
to a
.B LaTeX
format.
It is intended to do the first pass of the conversion. The user
should then finish up the rest of the conversion and customize the
converted manuscript to his/her liking.
It can also serve as a tutor for those who want to convert from
troff to LaTeX.
.PP
Most of the converted document will be in LaTeX
but some of it may
be in plain
.B TeX.
It will also use some macros in
.B troffms.sty
or
.B troffman.sty
which are included in the package and must be available to the document
when processed with LaTeX.
.PP
If there is more than one input file, they will all be converted into
one LaTeX document.
.PP
.B Tr2latex
understands most of the
.B -ms
and
.B -man
macros and
.B eqn
preprocessor symbols. It also understands several plain
.B troff
commands. Few
.B tbl
preprocessor commands are understood to help convert very simple tables.
.PP
When converting manuals, use the
.B -m
flag. Otherwise
.B tr2latex
assumes to translate a
.B -me
input file.
.PP
If a troff command cannot be converted, the line that contain that
command will be commented out.
.PP
NOTE: if you have
.B eqn
symbols, you must have the in-line mathematics delimiter defined by
.B delim
in the file you are converting. If it is defined in another
setup file, that setup file has to be concatenated with the
file to be converted, otherwise
.B tr2latex
will regard the in-line math as ordinary text.
.SH OPTIONS
.TP 12
.B \-m
Convert manual pages. This makes
.B tr2latex
understand most of the
.B -man
macros. It uses the style file
.BI troffman.sty .
.TP
.B \-t
Produce twoside page style.
.TP
.B \-\fIn\fP
Use a font size of
.IR n pt.
The default font size is 12pt for
.B \-man
and 11pt otherwise.
.TP
.B "\-s \fIstyle\fP"
Use the style file
.I style
instead of the default
.BR article.sty .
.TP
.B "\-o \fIoutfile\fP"
Write output to file
.IR outfile .
.SH BUGS
Many of these bugs are harmless. Most of them cause local errors
that can be fixed in the converted manuscript.
.PP
\- Some macros and macro arguments are not recognized.
.PP
\- Commands that are not separated from their argument by a space are
not properly parsed (e.g .sp3i).
.PP
\- When some operators (notably over, sub and sup) are renamed (via define),
then they are encountered in the text,
.B tr2latex
will treat them as
ordinary macros and will not apply their rules.
.PP
\- rpile, lpile and cpile are treated the same as cpile.
.PP
\- rcol, lcol are treated the same as ccol.
.PP
\- Math-mode size, gsize, fat, and gfont are ignored.
.PP
\- lineup and mark are ignored. The rules are so different.
.PP
\- Some troff commands are translated to commands that require
delimiters that have to be explicitly put. Since they are
sometimes not put in troff, they can create problems.
Example: .nf not closed by .fi.
.PP
\- When local motions are converted to \\raise or \\lower, an \\hbox
is needed, which has to be put manually after the conversion.
.PP
\- 'a sub i sub j' is converted to 'a_i_j' which TeX
parses as 'a_i{}_j}' with a complaint that it is vague. 'a sub {i sub j}'
is parsed correctly and converted to 'a_{i_j}'.
.PP
\- Line spacing is not changed within a paragraph in TeX
(which is a bad practice anyway).
TeX uses the last line spacing in effect in that paragraph.
.SH TODO
Access registers via
.B .nr
command.
.SH FILES
$(TEXLIB)/macros/troffman.sty
.br
$(TEXLIB)/macros/troffms.sty
.SH AUTHOR
Kamal Al-Yahya, Stanford University
.br
Christian Engel, Aachen University of Technology