.TH XWAISQ 1 "Thu Feb 20 1992" "Thinking Machines"
.SH NAME
xwaisq \- X interface to the WAIS questions
.SH SYNOPSIS
.B xwaisq
[ standard X toolkit arguments ] [ \fIQuestionName\fR ]
.SH DESCRIPTION
.PP
.B xwaisq
is an X Windows based client for the WAIS system.  It is the application
used for modifying and submitting questions to a WAIS server.
.PP
This application is heavily dependant on the X Toolkit resource mechanism.
Questions are stored to the user's "questionDirectory" (usually
~/wais-questions) and can be retrieved using this application.
.PP
A question is posed in three parts:
.PP
1.  Search words, labeled "Tell me about"
.br
2.  Sources, labeled "In Sources"
.br
3.  Relevance Feedback, labeled "Similar to"
.PP
Search words can be entered into the "Tell me about" frame by simply typing
into the frame.  Pressing return will initiate the search (however, if you
have not yet specified a source, you'll be warned, and no search will
actually occur).
.PP
To add a source  to  a question, press  "Add Source", and select the source
from the popup  menu.  To remove  a source  from  the question,  select the
source, and click on the "Delete Source" button.
.PP
To add a document to a question, select the document in the Result Document
list, and press "Add Document".   To remove  a document  from the question,
select the document, and click on the "Delete Document" button.
.PP
When  you have you question completely  defined, you can search the sources
by pressing the "Search" button, or <RETURN> in the "Tell me about" window.
.PP
To view  a resulting  document, select the  document and  press  the "View"
button.  A  new Document  window  will appear.    If you'd like to see  the
keywords matched, press the "Find Key" button.   To save the document  to a
file, press  the "Save to  File"  button.  When  you are finished  with the
document, press "Done" button in the Document window.
.PP
When you are finished with the question close its window using  the "Done"
button.  Any changes to the question  will be written  back to the question
directory.
.SH FILES
.ta 3i
.PP
Xwais	application Default file.
.br
$(HOME)/wais-questions	location of saved questions.
.Pbr
$(HOME)/wais-sources	location for personal sources.
.SH RESOURCES
.ta 3i
.PP
The following resources may be changed to modify the default behavior of
the application.  I recommend you look at the resource file Xwais to get a
better idea of all the resources available.
.PP
Xwais*questionDirectory:		~/wais-questions/
.br
Xwais*documentDirectory:		~/wais-documents/
.br
Xwais*userSourceDirectory:		~/wais-sources/
.br
Xwais*commonSourceDirectory
.br
Xwais*seedWords
.br
Xwais*initialSource
.br
Xwais*questionName
.br
Xwais*DoSearch
.PP
The last four resources are usually unspecified.  Under special
circumstances they may be used to fully specify a question without a file.
.SH BUGS
.PP
\fBxwaisq\fR is heavily dependant on it's resources.  If it is unable to
find the proper resource file, it will most likey become unusable.
.PP
\fBxwaisq\fR has no limit to resizing, so if the window is resized too
small, some buttons may be lost.
.SH SEE ALSO
.PP
xwais(1), waisindex(1), waisserver(1)
.SH AUTHOR
Jonny Goldman <Jonathan@Think.COM>
