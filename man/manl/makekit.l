.TH MAKEKIT 1l
.\" $Header$
.SH NAME
makekit \- split files up into shell archive packages
.SH SYNOPSIS
.B makekit
[
.B -1
] [
[
.B -e
] [
.BI -h\| #
] [
.BI -i\| name
] [
.BI -k\| #
] [
.B -m
] [
.BI -n\| name
] [
.BI -o\| name
] [
.B -p
] [
.BI -s\| #[k]
] [
.BI -t\| text
] [ file... ]
.SH DESCRIPTION
.I Makekit
reads a list of files and directories, determines their sizes,
and parcels them up into a series of shell archives such that all the
archives are of reasonable size.
It then invokes
.IR shar (1l)
to actually create the archives.
.PP
By default, no archive will be larger than about 50,000 bytes; this may be
changed by using the ``\-s'' option.
If the number given with the ``\-s'' option ends with the letter ``k''
then the size is multiplied by 1024, otherwise it is taken to be the
desired maximum size, in bytes.
Each archive will have a name that looks like
.IR Part nn,
where ``nn'' represents the two-digit sequence number (with leading zero
if needed).
The leader part of the archive name may be changed with the ``\-n'' option.
The ``\-n'' is also useful when write permission to the directory being
archive is denied; e.g., ``\-n/tmp/KERNEL.''
.PP
.I Makekit
reads its list of files on the command line, or standard input
if none are given.
It is also possible to specify an input filename with the ``\-i'' option.
The input should contain a list of files, one to a line, to separate.
In addition, if each input line looks like this:
.RS
filename\ \ \ whitespace\ \ \ optional-digits\ \ \ whitespace\ \ \ text
.RE
then
.I makekit
will ignore the spaces and digits, but remember the text associated with
each file, and output it with the filename when generating the ``shipping
manifest.''
Further, the ``\-h'' option may be given to have the program skip the
indicated number of lines in the input; this option is provided so that
.I makekit
can more easily re-parse the manifests it has generated.
.PP
The generated manifest will be sent to the standard output.
An alternate output file may be given by using the ``\-o'' option; if
the output file exists,
.I makekit
will try to rename it with an extension of
.IR \&.BAK \&.
If the ``\-o'' option is used,
.I makekit
will add that name to the list of files to be archived; the ``\-e''
option may be given to exclude the manifest from the list.
.PP
The ``\-m'' option is the same as giving the options,
\&``-iMANIFEST -oMANIFEST -h2.''
This is a common way to regenerate a set of archives after the first
use of
.I makekit
in a directory.
.PP
If a large number of kits has to be generated, you may need to give
the ``\-k'' option to increase the maximum number of kits to be
generated.
.PP
After partitioning the files and directories,
.I makekit
calls
.I shar
with the proper options to generate archives in a series.
Each resultant archive will, when executed, check to see if all the parts
are present.
If the ``\-1'' option is used, then
.I makekit
will not instruct
.I shar
to generate the checks (by not passing on the ``\-n'' and ``\-e'' options).
By using the ``\-t'' option, you can specify a line of starting instructions
to display to the recipient when all pieces have been unpacked.
This is useful when resending part of a series that has probably already
been unpacked by the recipient.
See
.I shar
for more information on multi-part archives.
If the ``\-x'' option is used,
.I shar
is not called, but the manifest is still created.
.PP
.I Makekit
normally reorders its input so that the archives are as ``dense'' as
possible, with the exception that directories are given priority over
files, and a file named
.I README
is the first of all.
The manifest is also sorted in alphabetical order; this makes it easy
to locate ``missing'' files when the distribution is a large one.
The ``\-p'' option may be used to override both sortings, however,
and preserve the original order of the input list in generating
both the manifest, and the shell archives.
.SH NOTES
.I Makekit
tries to partition the files so that all directories are in the first archive.
This usually means the first archive must be the first one to be unpacked.
.PP
Saying ``the `\-k' option is to help prevent runaway packaging'' is probably
.I "post hoc propter hoc"
reasoning.
.SH "SEE ALSO"
findsrc(1l), shar(1l)
