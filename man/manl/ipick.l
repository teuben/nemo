.\ " $Id$
.\ "
.TH IPICK 1 "28 February 1993"
.\ "
.\ " $Log$
.\ " Revision 1.2  2010/12/04 04:59:45  pteuben
.\ " long forgotten update
.\ "
.\ " Revision 1.1.1.1  2004-09-22 02:33:14  pteuben
.\ " import from 1993 comp.sources.misc
.\ "
.\ " Revision 1.1  1993/02/27  21:53:58  markd
.\ " Initial revision
.\ "
.\" Had to fix examples using -i to use -v (new invert option)
.\" I have no *roff manual so this manpage may well be faulty.
.\" Having said that, the good bits of this manpage come care of
.\" DaviD W. Sanderson <dws@ssec.wisc.edu>
.SH NAME
ipick \- A screen-based filter to
.IR i nteractively
.I pick
lines
.SH SYNTAX
.B ipick
.RB [ \-abdhrRvV ]
.RB [ \-m
.IR minimum ]
.RB [ \-M
.IR maximum ]
.if n .ti +.6i
.RB [ \-t
.IR fixed-title-text ]
.RB [ \-T
.IR stdin-title-lines ]
.if n .ti +.6i
.RB [ \-X
.IR xterm-name-substring ]
.RI [ filename ]
.SH DESCRIPTION
.PP
.I ipick
reads lines of text from the standard input or the optional
.I filename
and uses
.IR curses (3)
to present them as a full-screen selection list.
.I ipick
provides numerous commands
to select, navigate, scroll and search
through this list.
On quitting,
.I ipick
writes the selected lines to the standard output.
.PP
Typically you would use
.I ipick
as a final-filter to glue your neat, pre-stored
pipelines and scripts together in a friendly way
so that people other than Unix-o-philes can use them.
.\"
.\" Make the description brief so that the OPTIONS show up early.
.\" For those interested, there is a more details MOTIVATION section.
.\"
.PP
See
.BR MOTIVATION ,
towards the end, for a more detailed discussion.
.SH OPTIONS
.PP
Options can appear in any order so long as they
.I precede
the
.IR filename .
.TP
.B \-a
Automatically exit when the number of lines selected is within the
minimum and maximum values allowed
(see
.B \-m
and
.B \-M
options).
.TP
.B \-b
Activate the audible alarm on invalid keystroke commands.
Not normally needed as
.I ipick
always generates an error message.
Normally you would set this option only for inexperienced users.
.TP
.B \-d
Drain the standard input on exit.
This avoids the possibility of upstream processes
receiving a SIGPIPE.
This is more a nicety than a necessity.
Furthermore,
using this option could prove expensive if the upstream process
is a long way from finishing!
.TP
.BR \-h , " \-?"
Print an extended help message describing these options.
.TP
.B \-m
Minimum number of lines that the user must select before
.I ipick
will exit.
.TP
.B \-M
maximum number of lines that the user may select select before
.I ipick
will exit.
.TP
\&
You would typically use the minimum and maximum settings
in conjunction with the
.B \-a
and/or
.B \-r
options to ensure an orderly and
predictable outcome of the picking process.
.TP
.B \-r
Restricted mode.
In this mode the user cannot escape the clutches of
.I ipick
except as defined by the other command-line options.
This disables the shell and pipe commands.
Additionally, this disables keyboard signals by
setting the terminal to raw mode instead of cbreak mode.
Note that raw mode has other side-effects; see
.IR stty (1)
for more details.
.TP
.B \-R
Do not set terminal input to raw mode.
The main effects of this option
is that the
.IR stty (1)
control character sequences remain active and are thus capable of
generating
keyboard signals such as SIGINT and SIGQUIT.
Note that
.I ipick
has an abort command that can be bound to your normal ``intr'' or
``quit'' characters; see
.IR stty (1)
for more details on raw mode.
.TP
.BI \-t " fixed-title-text"
Text to use as a title on the screen.
.I ipick
creates a ``fixed title'' containing the given text
starting at the top line of the screen.
.I ipick
does not scroll this part of the title horizontally with the data.
If the
.I fixed-title-text
contains embedded newlines,
.I ipick
handles them correctly.
.TP
.BI \-T " stdin-title-lines"
.I ipick
will read
.I stdin-title-lines
lines from the standard input and use these
as the ``variable title'' which follows the ``fixed title''.
This option is especially useful when the upstream process
generates a title line,
such as
.IR ps (1),
or
.IR w (1).
.I ipick
scrolls this part of the title horizontally with the data.
.TP
\&
If the standard input is less than
.I stdin-title-lines
long, then
.I ipick
terminates silently.
.TP
\&
If you define a title with either
.B \-t
or
.BR \-T ,
then
.I ipick
separates the title from the data with a line of hyphens.
.TP
.B \-v
Invert the selections.
With this option,
.I ipick
writes the lines the user
.I did not
select to stdout.
.TP
.B \-V
Print a message containing the version,
compilation options, location of the system-wide
customization file and the obligatory copyright message.
.TP
.BI \-X " xterm-name-substring"
.I ipick
compares the
.I xterm-name-substring
with the terminal type defined in $TERM.
If it matches,
.I ipick
sends the escape sequences needed to enable and disable mouse-tracking.
The default
.I xterm-name-substring
is ``xterm'', so real
.IR xterm (1)
users need do nothing.
The comparison is case-sensitive.
See
.B "XTERM MOUSE TRACKING"
for more information.
.PP
The
.BR \-m ,
.B \-M
and
.B \-r
options allow you to carefully control the actions
and predict the results of a user.
For example, if you use
.TP
\&
ipick \-m 1 \-M 1 \-r
.PP
then you can write the rest of the pipeline assuming that
.I ipick
will produce one and only one line of output.
.SH "KEYBOARD COMMANDS"
.PP
.I ipick
accepts numerous keyboard commands,
perhaps the most important being
.I ?
which provides online help.
.PP
.I ipick
is ecumenical regarding keybindings
as it implements a reasonable set of
.IR emacs ,
.I vi
and
.I more
keybindings concurrently!
Additionally,
.I ipick
binds all of the commonly used commands to function keys as defined
in the TERMINFO (or the corresponding termcap) definition file.
This means that users can avoid learning the idiosyncratic nature of
.I vi
and
.I emacs
keystrokes.
.PP
Finally,
.I ipick
can be customized
with both a system wide startup file and
and a user startup file.
See
.BR CUSTOMIZATION ,
for a detailed discussion of these files.
For the purposes of this section it is assumed that
no customization is in effect.
.PP
The list of keyboard commands uses the following symbols:
.nr n \w'^V, <Space>, K_npage, K_index\0\0'+1n-1/1n
.TP \nn
^A \(em control-A
Press Control and ``A'' together
E-A \(em Escape-A
Press Escape followed by ``A''
.TP
.RI K_ tttt
identifies
the TERMINFO (or termcap)
.RI key_ tttt
used.
This is the ``Variable name'',
.I not
the ``Capname'' in the man page.
.PP
.I Selection Commands
.TP \nn
<CR>, K_ent
Toggle the selection state of the
current line, then move to the next line if present
.PD 0
.TP
0\-9
Select the specific line number
.TP
S
Select all lines
.TP
C
Clear all lines
.TP
+
Toggle the state of unread lines
.TP
s, K_select
Select the ``line-range'' nominated
.TP
c, K_clear
Clear the ``line-range'' nominated
.TP
t
Toggle the selection state of the ``line-range'' nominated
.PD
.PP
The ``line-range'' given to the ``s'', ``c'' and ``t'' commands may
consist of any one of:
.PD
.TP 4
\(bu
A number, series of numbers or range of numbers
.nf
(e.g. 1 5 7\-14 21\-19)
.fi
.TP
\(bu
The string ``visible'' meaning all lines currently visible on the
screen.
.TP
\(bu
The string ``all'' meaning all lines.
.PD
.PP
You may shorten both ``visible'' and ``all'' to just ``v'' and ``a''
respectively.
.I ipick
ignores the case of these strings.
.PD
.PP
.I Positioning commands
.TP \nn
.PD 0
T, K_beg
Top of File
.TP
B, K_end
Bottom of File
.TP
H, K_home
Top of screen
.TP
L, K_ll
Bottom of screen
.TP
^N, j, K_down
Next line
.TP
^P, k, K_up
Previous line
.PD
.PP
.I Vertical scrolling
.TP \nn
^U
Up half the screen
.PD 0
.TP
^D
Down half the screen
.TP
E-v, b, K_ppage, K_rindex
Up full screen
.TP
^V, <Space>, K_npage, K_index
Down full screen
.PD
.PP
.I Horizontal scrolling
.TP \nn
^B, h, K_left
Scroll left one character
.PD 0
.TP
^F, l, K_right
Scroll right one character
.TP
^I (TAB)
Scroll right one tabstop
.TP
E-i, K_cbt
Scroll left one tabstop
.TP
^A, ^ (circumflex)
Scroll to beginning of line
.TP
^E, $
Scroll to end of line
.TP
<
Scroll left half screen
.TP
>
Scroll right half screen
.PD
.PP
.I Searching
.TP \nn
^S, E-s, /, K_find
Forward search
.PD 0
.TP
^R, E-r, \e
Reverse search
.TP
n
Redo forward search
.TP
N
Redo reverse search
.TP
*, K_next
Find next selected line
.TP
&, K_prev
Find previous selected line
.PD
.PP
.I Miscellaneous
.TP \nn
\&., K_redo
Redo last command
.PD 0
.TP
g, K_move
Go to line
.TP
q, Q, ^X^C, ZZ, K_exit
Quit
.TP
^C
Abort
.TP
^L, K_refresh
Refresh
.TP
?, K_help
Provide online help
.TP
!
Shell command
.TP
|
Pipe command.
Pipe the current line into the command
.PD
.SH CUSTOMIZATION
.PP
When
.I ipick
is invoked, it reads initialization commands from the
site customization
file (normally /usr/local/lib/ipickrc but see the
.B \-V
option)
and
then the user's customization file .ipickrc in the
user's home directory.
Both files are optional.
The file can contain blank lines and comments where comments are
denoted by the first occurence of the ``#'' character.
.PP
A command must be completely contained on a single
line as there is no continuation character.
Each token can be separated by one or more
white-space characters.
A token can be quoted if it requires embedded spaces with either a
single or double quote.
.PP
.PD
.I Substitution sequences
.PP
The following special substitution sequences are recognized:
.ft I
.nr n \w'^char\0\0\0'+1n/1n
.ft
.TP \nn
.B \eb
A backspace character
.TP
.B \ef
A form-feed
.TP
.B \ee
An escape
.TP
.B \en
The newline character
.TP
.B \er
A carriage-return
.TP
.B \es
A space
.TP
.B \et
A TAB
.TP
.B \e\e
The back-slash character
.TP
.B \eNNN
The character represented by the  octal value of
.I NNN
.TP
.B ^char
A control character.
.I char
can be ``A-Z'' or ``@''.
.PP
.I Customization commands
.\" .nr n \w'bind-terminfo'+1n/1n
.\" .TP \nn
.PP
.B bind-key
.I function
.I keysequence
.RI [ helptext ]
.TP
\&
Bind a set of keystrokes to a specific function.
.PP
.B bind-terminfo
.I function
.I terminfo-capability
.RI [ helptext ]
.TP
\&
Bind a terminfo capability to a specific
function.
If the terminfo functionality has not been built into
.I ipick
then this effectively becomes a no-op.
.PP
.B bind-termcap
.I function
.I termcap-capability
.RI [ helptext ]
.TP
\&
Bind a termcap capability to a specific
function.
If the termcap functionality has not been
built into
.I ipick
then this effectively becomes a no-op.
.PP
.B include
.I includefile
.PP
Where:
.PP
.ft I
.nr n \w'terminfo-capability\0\0\0'+1n/1n
.ft
.TP \nn
.I function
is a valid
.I ipick
function; see
.BR FUNCTIONS ,
for a complete list.
.TP
.I keysequence
is any sequence of keystrokes.
The
.I substitution
character will typically be useful for this parameter.
.TP
.I helptext
is an optional text string that will be displayed as part of the online
help screen.
If not present then the keysequence or capability name will be used.
.TP
.I terminfo-capability
A valid terminfo keystroke capability name; see
.IR terminfo (5)
for more details.
.TP
.I termcap-capability
A valid termcap keystroke capability name; see
.IR termcap (5)
for more details.
.TP
.I includefile
a path of the file to include.
.TP
\&
If
.I includefile
has a leading tilde ``~'' then
.I ipick
performs the usual tilde expansion.
That is the string between
the tilde and the first ``/'' character
is treated as a login name whose home directory is substituted.
An empty login name implies $HOME.
.TP
\&
If
.I includefile
contains any shell characters then
.I ipick
assumes that it's a pipeline command
and uses
.IR popen (3S)
rather than
.IR fopen (3V)
to open the command.
.TP
\&
If
.I includefile
is treated as a pipeline command and there is no trailing pipe character ``|''
then
.I ipick
preceeds the pipeline string with
the
.IR cat (1)
command.
.PD
.PP
The follow sample
.I .ipickrc
file demostrates most of the functionality
discuss in this section.
.PP
.nf
#
# This file has some comments in it
#

include "~fred/.ipickrc"         # Fred has a cool startup file
include '/usr/local/lib/ipick_${EDITOR:-EMACS}'

include "grep terminfo /usr/local/lib/ipick_wyse |"

include "$HOME/env/test_bindings"

bind-key QUIT "\ee\es\e023" "NT250 Exit key"

# Bind the begining of line terminal definition

bind-terminfo beginning-of-line kbeg
bind-termcap beginning-of-line @1
.fi
.SH FUNCTIONS
The follow is a list of available functions.
.PP
.ft I
.nr n \w'scroll-right-screen\0\0\0'+1n/1n
.ft
.TP \nn
.I select-range
Enter a range of line numbers to pick
.TP
.I select-all
Mark all lines as picked
.TP
.I clear-range
Enter a range of line numbers which are
to have their pick state cleared
.TP
.I clear-all
Mark all lines as not picked
.TP
.I toggle-current
Toggle the picked state of the current
line and move the cursor to the next line
.TP
.I toggle-range
Enter a range of line numbers which are
to have their pick state toggled
.TP
.I toggle-unread
Toggle the picked state of all unread lines
.TP
.I top-of-screen
Move the cursor to the top of the screen
.TP
.I bottom-of-screen
Move the cursor to the bottom of the screen
.TP
.I previous-line
Move the cursor up one line
.TP
.I next-line
Move the cursor down one line
.TP
.I quit
Exit and write the pick lines to stdout
.TP
.I abort
Exit without writing anything to stdout
.TP
.I help
Display the online help screen
.TP
.I refresh
Clear and re-display the current screen
.TP
.I scroll-left-char
Move the cursor left one position
.TP
.I scroll-right-char
Move the cursor right one position
.TP
.I beginning-of-line
Move the cursor to the first character position
.TP
.I end-of-line
Move the cursor to the last character position of the current line
.TP
.I scroll-tab
Move the cursor forward one tabstop
.TP
.I scroll-backtab
Move the cursor back one tabstop
.TP
.I scroll-left-screen
Move the cursor left by one half the width of the screen
.TP
.I scroll-right-screen
Move the cursor right by one half the width of the screen
.TP
.I scroll-up-half
Move the cursor up by one half of the height of the screen
.TP
.I scroll-down-half
Move the cursor down by one half of the height of the screen
.TP
.I scroll-up-full
Move the cursor up a full screen
.TP
.I scroll-down-full
Move the cursor down a full screen
.TP
.I beginning-of-file
Move the cursor to the first line
.TP
.I end-of-file
Move the cursor to the end of all
the input lines.
If data is from an upstream process, ipick will read until
this process writes end of file
.TP
.I search-forward
Enter the search forward string
.TP
.I search-backward
Enter the reverse search string
.TP
.I re-search-forward
Continue the search forwards
.TP
.I re-search-backward
Continue the search backwards
.TP
.I next-selected
Move the cursor to the next line selected
.TP
.I previous-selected
Move the cursor to the previous line selected
.TP
.I goto-line
Enter a line number to position to
.TP
.I shell
Enter a shell command
.TP
.I pipe
Enter a  command which will have the
current line written to it's stdin
.TP
.I redo-command
Redo the last command
.TP
.I xterm-mouse
Define the xterm escape prefix.
Only to be used if you really know what you're
doing (Which I don't when it comes to
.IR nroff (1))
.TP
.I invalid-command
Treat as an invalid command and generate an error message
.SH "XTERM MOUSE TRACKING"
.PP
.I ipick
has limited support for
.IR xterm 's
``mouse tracking'' capability (the X11
``normal tracking mode'' not the X10 compatibility mode).
.PP
To enable this facility in a particular
.I xterm
you have to send it a special escape sequence.
To quote from
.IR xterm (1),
``[mouse tracking] is enabled
by specifying parameter 1000 to DECSET''.
.I ipick
does this automatically when it detects a terminal type of
.IR xterm .
If your
.I xterm
clone uses a different name, you can use the
.B \-X
option to tell ipick what it is.
.PP
When you enable this facility, your
.I xterm
will pass any mouse events to
.I ipick
in a form that
.I ipick
recognizes.
In all cases, the down event defines the start of a range of lines and
the up event defines the end of the range \(em so dragging is useful.
.PP
Each mouse button has the following function:
.nr n \w'Button\0\0\0'+1n-1/1n
.TP \nn
.B Button
.B Description
.TP
1
Set the selected status of the range
.PD 0
.TP
2
Toggle the selected status of the range
.TP
3
Clear the selected status of the range
.PD
.PP
Note that
.I ipick
ignores augmentation of mouse events
with the Shift, Control and Meta keys.
.PP
The
.I xterm
facility is limited in that chording the mouse buttons seems
to be undefined.
Furthermore, the release (or mouse-up) doesn't specify the button.
Accordingly
.I ipick
takes a conservative approach to mouse-events and tends to discard
anything unexpected.
.ne 11
.SH DIAGNOSTICS
All error messages should be self-explanatory.
.PP
.I Exit codes.
.TP
0
Normal termination \(em at least one line selected
.TP
1
Normal termination \(em no lines selected
.TP
2
Abnormal termination
.SH NOTES
.PP
If no input exists, or if the
.B \-T
option causes
.I ipick
to consume all its input, then
.I ipick
terminates silently with an exit code of 1.
.PP
Each Newline terminated string in
.I fixed-title-text
(see the
.B \-t
option)
becomes a separate line on the fixed-title section of the screen.
.PP
.I ipick
only reads lines from the standard input as needed
(and after each keyboard command)
rather than reading all input lines on starting.
This is especially useful if the upstream process generates output lines
slowly as
.I ipick
is able to display lines as soon as they become available \(em within
the constraints of any buffering.
It is also useful if the upstream process has the potential to generate
an enormous number of lines of output prior to completion.
.PP
.I ipick
processes binary files correctly, but it's hard to imagine that
this capability is especially useful.
.PP
When constructing pipelines, be aware of the fact that many
commands don't take multiple parameters.
In these cases, use xargs \-l1 if your system has it.
.SH RESTRICTIONS
.PP
.I ipick
is designed to process modest amounts of data.
The data is held in memory and all functions are coded simplistically.
If you ask
.I ipick
to handle large amounts of data, it does so sluggishly and
consumes excessive system resources.
.PP
.I ipick
takes a passive approach to ambiguous function key
definitions \(em later definitions override earlier definitions.
.PP
The search function does not handle regular expressions.
.PP
Because of the way in which
.I ipick
reads the standard input, using
.I ipick
with data from the keyboard does not work
as you would want (In fact
.I ipick
should probably insist on a pipe or a file as input,
just as
.I more
does).
The workaround is to use the ``here document'' capability
of the shell (the ``<<'' redirection).
.PP
.I ipick
arbitrarily expands tab characters to 8-column tabstops.
.SH BUGS
.PP
The online help does not display properly if the screen is less than 80
columns wide.
.PP
Directions implied by movement and scrolling commands apply to
the cursor, not the data.
.PP
.I ipick
does not know when include files opened with
.IR popen (1)
fail.
This can sometimes cause
.I ipick
to wait for ever depending on the state of the child process.
.SH MOTIVATION
.PP
It is the very essence of Unix to make useful
commands with the generic construct:
.IP
generate_listing |
.I FILTER
| do_something
.PD 0
.PP
or
.IP
do_something `generate_listing |
.IR FILTER `
.PD
.PP
The problem is
.IR FILTER .
Getting it correct for the simplistic case is easy,
making it perfect and bullet-proof is not.
This is especially true if the pipeline is being
developed for the user community.
.PP
Think about how you typically build a pipeline for
the following requests:
.IP
``Kill my awk process, it's hung my terminal.''
.br
``Remove my print job, I've run the wrong report.''
.br
``Remove that message queue, then re-run the daemon.''
.PP
Either you construct a nice obscure
.I FILTER
using some
combination of
.IR grep (1),
.IR awk (1),
.IR perl (1),
or
.IR sed (1)
\(em usually
after
you've had a look at the
.I generate_listing
output a couple of times to make sure
you don't zap the wrong thing!
Alternatively you run the
.I generate_listing
program a couple of times until you've memorized the
relevant identifier (such as pid, job number, queue id), then you
run the
.I do_something
program and re-type the relevant identifier
trying as best you can to avoid mis-typing and mis-remembering.
.PP
In other words, tedious and error-prone.
.PP
Of course, when the
time comes to give your neat pipeline or script
to the user community, how do you give it an
easy to use, safe, bullet-proof interface?
Do you knock up a quick shell wrapper with
token prompts, perfunctory checking and an
interface that's
.I almost
the same as your last shell wrapper?
.PP
If any of these situations sound familiar then
.I ipick
may well be your pipeline panacea!
(Well, at least marginally useful.)
.PP
This is because
.I ipick
makes the user the final part of the
.I FILTER
in a safe, friendly manner, often obviating the need for
shell wrappers and such.
.PP
The example of selecting and removing a set of print jobs
is the easiest way to demonstrate
.IR ipick .
With the pipeline:
.nf
.IP
lpq | grep `whoami` | ipick | awk '{ print $3 }' | xargs lprm
.fi
.PP
You use
.I ipick
to select the print jobs to be removed and the
pipeline does the rest, removing only those print jobs
selected with
.IR ipick .
Another example:
.nf
.IP
kill \-9 `ps | ipick \-T1 | awk '{ print $1 }'`
.fi
.PP
you simply select the lines with the relevant pids then
exit from
.I ipick
\(em the pipeline does the rest by killing only those processes
you selected with
.IR ipick .
.SH EXAMPLES
.PP
A few more examples to get you started.
Normally you would define each of these
as a function or alias in your shell.
(Of course, I present these examples as ideas.
They are not complete, bullet-proof functions.)
.PP
.I "Pick source files to edit."
.nf
.IP
vi `ls \-1 *.c | ipick \-m1`
.fi
.PP
.I "Pick a directory to cd to."
.nf
.IP
cd `ls \-l | grep '^d' | ipick \-m1 \-M1 \-a | cut \-c46\-`
.fi
.PP
.ne 6
.I "Pick files to extract from a tar file in the default tape drive."
.nf
.IP
tar t | ipick | xargs tar xv
.fi
.PP
I find this especially useful when the archive
contains l-o-o-o-ng filenames.
.PP
Actually, in its current form
the above example has a number of limitations, so a more complete
solution to an interactive tar is:
.nf
.IP
.ta \w'\0\0\0\0\0\0\0\0'u
tar tvf ${1:\-/dev/rmt8} | sed \-e 's./$..' | ipick |
.br
	cut \-c42\- | xargs \-v \-t tar xvf ${1:\-/dev/rmt8} {}
.fi
.DT
.PP
.I "Clean up a directory containing many junk files."
.PP
The following pipeline does this by letting you select
the files you want to keep, and deleting the rest.
.nf
.IP
.ta \w'\0\0\0\0\0\0\0\0'u
rm `ls \-l | grep '^\-' | ipick \-v \-t "Which files to keep ?" |
.br
	cut \-c46\-`
.fi
.DT
.PP
.ne 19
.I "Part of the login script to set the terminal type"
.nf
.IP
.ta \w'xterm\0\0\0\0'u
#! /bin/sh
.sp
\&...
.sp
export TERM
TERM=`ipick \-m1 \-M1 \-r \-a \-T3 <<EOD | cut \-f1 \-d' '
.sp
	Pick the terminal type that you are logged into
.sp
vt100	The old grey terminals in the conference room
xterm	The new fancy terminals in the bosses office
sun	One of the workstations in the sysadmin's office!
EOD`
.sp
\&...
.fi
.DT
.PP
.I "Site conventions for predefined files"
.PP
If you get really keen, you can have a site convention for
all pre-defined
.I ipick
files, such as:
.TP
\(bu
first two lines are always header
.TP
\(bu
first space separated field is always the output
selection value
.PP
Then you can define a generic function or script
(let's call it ipickf) such as:
.nf
.IP
ipick \-a \-m1 \-M1 \-r \-T2 $HOME/pickfiles/$1 | cut \-f1 \-d' '
.fi
.PP
then use it around the traps as:
.nf
.IP
TERM=`ipickf TERM`
.fi
.PP
.SH VERSION
.PP
.I ipick
version 1.1, dated 28 February, 1993.
.SH AUTHOR
.PP
Copyright (c) 1993, Mark Delany <markd@werple.apana.org.au>
All rights reserved.
.PP
Substantial man page improvements by DaviD W. Sanderson <dws@ssec.wisc.edu>
.PP
.I ipick
may only be copied under the terms of either the Artistic License
or
the GNU General Public Licence,
which may be found in the ipick source kit.
.SH "SEE ALSO"
.IR stty (1),
.IR xargs (1V),
.IR xterm (1L),
.IR curses (3V),
.IR terminfo (5),
.IR termcap (5)
