.TH IPICK 1 "1 August 1992"
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
.RB [ \-abdhrsv ]
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
.B \-s
Catch keyboard signals.
Normally
.I ipick
terminates on receiving a keyboard signal (SIGINT or SIGTERM).
With this option
.I ipick
catches SIGINT and aborts the current command.
This is only useful for long-running
commands and even then the problem is that all processes in the
same process group get the signal and typically terminate.
Hence, this option is marginally useful at best.
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
Print a message containing the version,
compilation options and the obligatory copyright message.
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
is able to display lines as soon as they become available - within the
constraints of any buffering.
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
There is no
.I .ipickrc
file to add or override the default keybindings \(em yet!
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
Arguably
.I ipick
could use
.IR select ()
and look ahead for keyboard signal characters rather
than rely on cbreak mode as this method would protect other
processes in the same process group.
However, this is less portable and messy.
.PP
Directions implied by movement and scrolling commands apply to
the cursor, not the data.
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
.I sed (1)
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
tar tvf ${1:-/dev/rmt8} | sed -e 's./$..' | ipick | \e
.br
	cut -c42- | xargs -i -t tar xvf ${1:-/dev/rmt8} {}
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
version 1.0 Gamma1, dated 1 August, 1992.
.SH AUTHOR
.PP
Copyright (c) 1992, Mark Delany <markd@werple.pub.uu.oz.au>
All rights reserved.
.PP
Substantial man page improvements by DaviD W. Sanderson <dws@ssec.wisc.edu>
.PP
.I ipick
may only be copied and used under the terms and conditions
of the COPYRIGHT notice found in the source kit.
.SH "SEE ALSO"
.IR stty (1),
.IR xargs (1V),
.IR xterm (1L),
.IR curses (3V),
.IR terminfo (5),
.IR termcap (5)
