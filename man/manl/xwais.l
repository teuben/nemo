.TH XWAIS 1 "Thu Feb 20 1992" "Thinking Machines"
.SH NAME
xwais \- X interface to the WAIS system
.SH SYNOPSIS
.B xwais
[ standard X toolkit arguments ]
.SH DESCRIPTION
.PP
.B xwais
is a simple interface that allows the selection and modification
of WAIS questions and sources.

There are basically three kinds of windows you will see in Xwais:
.PP
1. the main Xwais window.
.br
2. a number of Question windows.
.br
3. a Source window.
.PP
In the main Xwais window you'll see two lists - Questions and Sources.
.PP
To open a question, select the question in the Question List and press
the "Open" button in the  Question frame.  To create a  new  question,
press the "New" button in  the  Question frame.   Press "Help" in  the
resulting Question window    for  information  on how  to   manipulate
questions.
.PP
To add a new source, press the "New"  button in the  Source Frame.  To
edit an old source,  press  the "Open"  button.
.PP
To exit the application, click twice on the "Quit"  button.  This only
closes the Xwais "shell".  Questions must be closed explicitly.
.PP
Remember, your window  manager may be managing  some of these windows,
so if you  click on a  button and see  no window,  check to see if the
window has actually been iconified.
.SH FILES
.ta 3i
.br
Xwais	application Default file.
.br
$(HOME)/wais-questions	location of saved questions.
.br
$(HOME)/wais-sources	location for personal sources.
.SH RESOURCES
.ta 3i
.PP
The following resources may be changed to modify the default behavior of
the application.  I recommend you look at the resource file Xwais to get a
better idea of all the resources available.
.br
Xwais*questionDirectory:		~/wais-questions/
.br
Xwais*userSourceDirectory:		~/wais-sources/
.br
Xwais*commonSourceDirectory:		/proj/wais/wais/wais-sources/
.SH BUGS
.PP
\fBxwais\fR currently has a limit of twenty questions. It will arbitrarily
trunctate the list to the first twenty files it sees.
.PP
\fBxwais\fR does not sort the question or source lists.  It really should.
.PP
\fBxwais\fR is heavily dependant on it's resources.  If it is unable to
find the proper resource file, it will most likey become unusable.
.PP
\fBxwais\fR has no limit to resizing, so if the window is resized too
small, some buttons may be lost.
.SH SEE ALSO
.PP
xwaisq(1), waisindex(1), waisserver(1)
.SH AUTHOR
Jonny Goldman <Jonathan@Think.COM>
