.TH WAISINDEX 1 "Sun May 10 1992" "Thinking Machines"
.SH NAME
waisindex \- Indexes files
.SH SYNOPSIS
.B waisindex
[\ \-d\ index_filename\ ] [\ \-a\ ] [\ \-r\ ] [\ \-mem\ mbytes\ ] [\ \-register\ ] 
[\ \-export\ ] [\ \-e\ [\ file\ ]\ ] [\ \-l\ log_level\ ]
[\ \-pos\ |\ -nopos\ ] [\ \-nopairs\ |\ -pairs\ ] [\ \-nocat\ ] [\ \-T\ type\ ]
[\ \-t\ type\ ] [\ \-contents | \-nocontents\ ] filename filename ...
.SH DESCRIPTION
.PP
.B waisindex
creates an index of the words in files so that they can be searched quickly
(see waissearch).  The index takes about as much disk space as the original
text.  It also creates a new source structure named index_filename.src if none exists.
.SH OPTIONS
.PP
.TP 10
\fB\-d\fI index_filename
This is the base filename for the index files.  Therefore if /usr/local/foo
is specified, then the index files will be called /usr/local/foo.dct etc.
.br
The index should be stored on the local file system of the machine running
waisindex.  It works over NFS, but it is much slower.
.TP 10
.B \-a
Append this index to an existing one.  Useful for incremental additions or
updates.  This will only add onto an index, so that if a file has changed,
it will get reindexed, but the old entries will not be purged.  Therefore,
to save space, it is a good idea to reindex the whole set of files
periodically.
.TP 10
.B \-r
Recursively index subdirectories.
.TP 10
.B\-mem
How much main memory to use during indexing.  This variable will have a
large effect on how fast indexing is done.
.TP 10
.B \-register
Register this database with the directory of servers.  You are encouraged
to register databases, but only ones that will be consistently running.
The directory of servers is available to anyone that is on the internet or
can phone in.
.TP 10
\fB\-export
This causes the resulting source description file to include the host-name
and tcp-port for use by the clients.  Otherwise the file contains no
connection information, and is expected to be used only for local searches.
.TP 10
\fB\-e\fI [ filename ]
Redirect error output to pathname, if supplied, or to /dev/null.
Error output defaults to stderr, unless -s is selected,
in which case it defaults to /dev/null.
.TP 10
\fB\-l\fI log_level 
set logging level.  Currently only levels 0, 1, 5 and 10 are meaningful:
Level 0 means log nothing (silent).  Level 1 logs only errors and warnings
(messages of HIGH priority), level 5 logs messages of MEDIUM priority (like
indexing filename info).  Level 10 logs everything.
.TP 10
\fB\-pos (\-nopos)
Include (don't include - default) word position information in the index.
This will increase the index size, but will allow search engines to do
proximity.
.TP 10
\fB\-nopairs (\-pairs)
Don't build (build - the default) word pairs from consecutive capitalized
words.
.TP 10
\fB\-nocat
Inhibits the creation of a catalog.  This is useful for databases with a
large number of documents, as the catalog contains 3 lines per document.
.TP 10
\fB\-contents (\-nocontents)
Include (exclude) the contents of the file from the index.  The filename
and header will still be indexed.  Default is type depedant.
.TP 10
\fB\-T type
Sets the TYPE of the document to "type".
.TP 10
\fB\-t\fI type
This is the format of files that are handled by waisindex.  It is easy to
parse a different format, but that has to be done by changing the source
(ircfiles.c).  To find out the list of currently known types, execute the
waisindex command with no arguments and it will list them.
.TP 10
\fBfilename filename...
These are the files that will be indexed according to the arguments above.
To insure the files are registered in the filename table correctly, it is 
advised that these be full paths (beginning with a /).  If the database is
to be used from a machine other than the machine on which the index is
created, this should be a machine-independant path.

.SH SEE ALSO
waissearch(1), waisserver(1), waissearch-gmacs(1), xwais(1), xwaisq(1)
.PP
Wide Area Information Servers Concepts by Brewster Kahle.
.br
Brewster@think.com

.SH DIAGNOSTICS
.PP
The diagnostics produced by the waisindex are meant to be self-explanatory.

.SH BUGS
.PP
It temporarily takes twice the space it needs for an index.

Due to some compile time constants the document table is limited to 16
Megabytes.  This limits the indexer to databases with headlines that add up
to less than 16 megabytes (since thats the principal component of the
table).  This is typically a problem for database types where a record is
essentially a headline (one_line, archie).

See the note in ir/README in the wais distribution for more detail.
