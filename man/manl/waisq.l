.TH WAISQ 1 "Thu Feb 20 1992" "Thinking Machines"
.SH NAME
waisq
.SH DESCRIPTION
.PP
This is a shell level interface to WAIS questions that have been built
using xwais.  It allows a user to submit questions to sources for queries,
and to retrieve documents from questions (and the appropriate source).  It
can build or modify a question, including the modification performed during
a search.
.SH USAGE
.PP
.B waisq \\
.br
  [\-f question_file] /* defaults to stdin, stdout */
.br
  [\-s sourcedir] /* or ~/wais-sources/ */
.br
  [\-S sourcename]
.br
  [\-c common_sourcedir] /* or it is ignored */
.br
  [-m max_Result_Docs] /* defaults to 40 */
.br
  [\-v document_number]
.br
  [\-g] /* to do a search */
.br
  [\-h] /* this message */
.br
  [\-t] /* test a source */
.br
  [\-V] /* print version */
.br
  [word [word [...]]]
.br
or
.PP
\fB waisq \fR\- [other options] /* for interactive mode */
.SH OPTIONS
.TP 10
\fB\-f \fIquestion_filename
the pathname for a question structure as stored
in a file.  It uses standard UNIX conventions.  If no filename is provided,
stdin is used.  If the filename is a dash ( - ), no question will be read,
so you must specify a sourcename and the words to search.  This will create
a question structure and write it to stdout.

.TP 10
\fB\-s \fIsourcedir
a directory where source structures are stored.  It must end with a '/'.
If it is not specified, ~/wais-sources/ is used.  If none of those exist,
the application will probably fail, unless the common_sourcedir exists.
Only files that end in .src in the directory are considered sources.
.TP 10
\fB\-S \fIsourcename
the pathname of a source structure file.  It must exist in sourcedir or
common_sourcedir, and end in .src
.TP 10
\fB\-c \fIcommon_sourcedir
Similar to sourcedir, but usually common for a site. It need not exist.
.TP 10
\fB\-m \fImaxDocs
The limit on the number of results returned by a server.  Default value is
40.
.TP 10
\fB\-v \fIdocument_number
The number of the document to retrieve from a valid query.  This can be
from a previously stored question, or the current question if the -g option
is used.
.TP 10
\fB\-t
Test a source to see if it is servicing wais requests.  The first word will
be taken as the name of a source.
.TP 10
\fB\-g
Submit the query.  You must specify the -g option in order to actually
perform a query.  You need not specify it if you have already performed a
query from which you wish to retrieve a document.
.TP 10
\fB\-V
Print the version number and date of the program.
.TP 10
\fBword\fR...
optional search keywords to replace any search keywords in the question.
.TP 10
\fB\-
If the first argument is a dash ( \- ) then waisq operates in an
interactive mode.  It then reads question structures from standard input,
and returns results to standard output.  If the question has search-words
or relevant documents, it performs the query, and returns the resulting
question structure to standard output.  If there are only result documents,
it retrieves them and prints them to standard output.  It continues to read
it's input until an end of file is encountered, and then it exits.  All
other command line switches will work, except \-f.  This makes waisq useful
as a tool for user interfaces to access the WAIS system.
.SH EXAMPLES
.TP
to create a new question and submit it for a query:
.br
waisq -f - -S sourcename -g word1 word2 > my_question
.TP
to submit an existing question for query:
.br
waisq -f question_name -g
.TP
To retrieve a response:
.br
waisq -f question_name -v 1

.SH SEE ALSO
waisserver(1), waisindex(1), waissearch(1), xwais(1), xwaisq(1)
.PP
Wide Area Information Servers Concepts by Brewster Kahle.
.br
brewster@think.com

.SH DIAGNOSTICS
The diagnostics produced by the application are meant to be
self\-explanatory.

.SH BUGS

.SH AUTHOR
Jonny Goldman <Jonathan@Think.COM>
