.TH UUEXPLODE 1L "11 January 1994"
.IX uuexplode
.SH NAME
uuexplode - decode multiple and/or multipart uuencoded files to their
target files
.SH SYNOPSIS
\fBuuexplode\fR [-hqslcfw] [-t dirname] [filename, filename...]
.SH DESCRIPTION
\fBuuexplode\fR is a "smart", fast uudecoder that will decode one or
more files consisting of complete uuencoded files. The target files can be
written in any user specified directory or in current. It is upward compatible
with uudecode(1) provided that the encoder program used to create the input
to \fBuuexplode\fR has more or less followed the format described in
uuencode(5). See Technical Notes below.
.PP
If no filename(s) is supplied, input is read from stdin. It can also write to
stdout and hence act as a proper filter.
.PP
\fBuuexplode \fR automatically ignores lines that do not appear to be part of
the uuencoded file. This is particularly useful in uudecoding files that have
been concatenated from several email or news postings \fIwithout\fR stripping
off the headers or otherwise editing extraneous content. A lot of effort was
put into ensuring that erroneous target files are not written, several error
and warning messages are written to stdout when \fBuuexplode\fR finds
something in the input it doesn't accept or just finds strange. The messages
are meant to be self explanatory and contain line numbers in the input file.
Note that when warnings are issued, the output is still most likely correct.
.SH OPTIONS
.TP 8
\fB-h\fR
Print usage message stdout.
.TP
\fB-q\fR
Quiet. Do not print any diagnostics, errors or warning messages - no matter
what goes wrong.
.TP
\fB-s\fR
Write data to stdout instead of to the target file. Only the very first
target file found in the input will be written, then \fBuuexplode\fR exits.
By the "first" file is meant the first file actually containing any decodable
uucode at all. This file may not be in the first input file if there are more
than one. All warning/error messages go to stderr so there is no need to use
-q with -s.
.TP
\fB-l\fR
Lenient mode. This allows for the short data line, the one not beginning with
an "M"at the end of a uucode block, to have a non-standard length.
Use this as a second attempt if a target file seems to come out corrupted
even if no warnings were printed during decoding. Note: -l causes the security
of uucode recognition to drop quite considerably, do not use casually.
.TP
\fB-c\fR
Careless mode. This stops the checks for illegal characters in the short
data line. Its main use is in conjunction with -f when it allows decoding
of uucode in the character range [>..}] that some non-standard encoders
produce. This kind of code is transparent to a "dumb" uudecode program.
Never use this option unless you know you have non-standard uucode.
.TP
\fB-f\fR
Fast mode. This stops the checks for illegal characters in the normal uucode
data lines. It makes the decoding go faster, but not significantly (in fact
probably not even noticably). Normally this would only be used with -c
if you know you have uucode with a non-standard character range. The
sequrity of uucode recognition with -f -c is \fImuch\fR lower than normal.
.TP
\fB-w\fR
Overwite. Without this, \fBuuexplode\fR will refuse to overwrite a target file
that already exists.
.TP
\fB-t dirname\fR
Target directory. The named directory does not have to exist, it will be
created. Should this fail, \fBuuexplode\fR exits and no target files will be
written. Without -t, the target files always end up in current directory.
.SH TECHNICAL NOTES
This is for advanced users who know the uucode format well and want to
know more about what \fBuuexplode\fR accepts as valid uucode.
A number of constraints are put on it:
.TP 3
1
The first line after the 'begin'-line is a uucode line if all characters
are valid, it begins with [SPC..M,'] and the length is not totally
absurd with respect to what the first char was.
.TP 3
2
There is at least a \\n behind each uucode line. There may be a \\r also.
.TP 3
3
If there are extra chars at the end of the data lines it is present in
all lines, even the last `-line. And always the same number of chars.
.TP 3
4
There are no extra chars behind the 'end'-line except \\n or \\r\\n
.TP 3
5
The normal uucode lines are always 45 bytes (60 chars) and have 'M'
as length indicator. A "standard" uudecode program does not care about
this at all.
.TP 3
6
No valid uucode contains ASCII > 'a'. Note that one can also use
the range [>..}] since bit 6 of the chars are unused. Again a "standard"
uudecode program will happily decode such data, uuexplode will only do
this with -f -cc. There are uuencoders (e.g. on the Macintosh) that produce
such uucode.
.TP 3
7
The end of a uucode block resides in the same input file (or stdin)
as the corresponding 'begin'-line.
.TP 3
8
All parts belonging to the same target file(s) are in correct order in the
input.
.TP 3
9
The input contains no characers with ASCII code less than 32.
.PP
A small list of capabilities concerning garbage lines:
.br
No assumptions whatsoever are made about the format of a garbage line.
.br
Handles garbage between all the last three lines of a uucode block.
.br
Handles bogus SPACE/` lines between parts (and at the end) of the same block.
.br
Handles bogus short data lines in the middle of a uucode block.
.br
Handles short or even empty uucode blocks correctly.
.br
Tolerates optional checksum/"guard" character(s) at the end of data lines.
.br
Tolerates SPACE or ` as zero-code. The empty-line can be SPACE, ` or ``.
.SH RETURNS
If all decoding went well, the returncode is 0. Otherwise it is the
logical OR of the following values:
.br
1 \ \ At least one warning was issued
.br
2 \ \ At least one error occurred, possibly casuing a target file not to be 
created
.br
4 \ \ A target file already existed (only if -w was not actice)
.SH BUGS
Some of the warnings and errors are not quite as self explanatory as one would
wish. A user that does not know the uucode format and does not want to know,
may not be able to understand why the file he tried to decode failed.
.PP
The options -c and -l require at least rudimentary knowledge of the uucode
format to be really useful.
.PP
There is no proper standard formally specifying the format of uucode or what
a uudecode algorithm shall accept as valid input. This program tries to
be lenient enough to accept all the common variants while still being as
pedantic as possible.
There are variants of uuencoders out there producing (at least capable of
producing) uucode that uuexplode does not accept. A present no such
\fIcommonly\fR used uuencoders are known to me.
.PP
The algorithms used to detect valid uucode are not based on proper statistical
analysis of uuencoded data with garabage lines in it. Instead, some fairly
rigid theoretical analysis, partly based on what kind of garbage normally
occur in Email and News postings (but not limited to those) was done.
It \fImay\fR not work as well with totally random garbage interleaved with the
uucode.
.SH SEE ALSO
uuencode(5), uuencode(1), uudecode(1)
.SH AUTHOR
Michael Bergman (euambn@eua.ericsson.se)

Uuexplode was originally based on a program called "kiss", written by
Kevin Young. This program was part of a collection of uucode tools called
\fBuuxfer\fR. The actual decoding function that decodes one line of
uucode is also written by Kevin and altough all the other code has been
almost completely rewritten this remains unaltered.

