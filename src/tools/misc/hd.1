.TH HD 1 "January 28, 1986"
.UC
.SH NAME
hd \- dump the contents of a file in hex
.SH SYNOPSIS
.B hd
[ name... ]
.br
.SH DESCRIPTION
For each argument,
.I hd
dumps to standard output the contents of the file in hex
and ASCII, and the offsets into the file for each byte.
If no file is specified, the standard input is used.
If one of the files cannot be opened, a message is printed
to that effect and the next file is tried.
.PP
Dumps from
.I hd
are meant to go to a printer, so form feed characters
are inserted to cause page ejects at the end of the
page.  Pages are formatted as follows:  At the top of
the page is the name of the file being dumped and the
column headings.  The column headings label the offset
of the byte in that column from the address given at
the beginning of the line containing the byte.
.PP
After the column headings, the dump of the file is presented.
Each line contains 16 bytes of data from the file.  The first
item on the line is the address within the file of the
first byte of the line.  If the address is preceded by a "+",
then it is the beginning of a 512 byte block.
.PP
The next item on the line is the hexadecimal representation of the bytes
in the file.  For machines with the so\-called "little\-endian"
architecture (eg. VAXes, IBM PCs), the bytes are arranged right to
left, that is, the byte whose address is given at the beginning of the
line is at the rightmost part of the hex dump.  For machines with
"big\-endian" architectures (eg. Sun Workstations), the bytes are
arranged left to right. 
.PP
Following the hex dump is the ASCII dump.  It is always arranged
left to right.  Characters which are not printable (eg. control
characters) are represented by a ".".
.PP
When the end of file is reached,
.I hd
starts a new page.
.SH BUGS
On MS-DOS machines, most C standard I/O libraries do newline conversion,
that is, a \\r\\n sequence is turned into a single \\n.  This is
corrected when reading from files which are given on the command line,
but not when reading from stdin.
.PP
This manual page uses the term ASCII when it really shouldn't.  Hd will
compile and run properly on an EBCDIC machine, since <ctype.h> and
isprint() are used.
.SH FILES
.nf
.ta +2.0i
~src/nemo/util  	hd.c hd.1
.SH UPDATE HISTORY
.nf
.ta +2.0i
xx-xxx-xx	Version 1.0, created (Kenneth J. Montgomery)
28-Jan-87	Version 2.0, various mods (Donald L. Nash)
01-Feb-87	Version 2.0a, nextfile() bug, hexdump left to right (PJT)
.fi
