.TH CONV123 1 "6 Aug 1993" "123 Utilities"
.SH NAME
conv123 \- The ASCII to Lotus(TM) spreadsheet converter

.SH SYNOPSIS
.B conv123
[inputfile] [outputfile]

.SH DESCRIPTION

.I conv123
is a utility for reading a specially-formatted flat
ASCII file and converting it into a Lotus/123(TM) ".wk1" file.
Each specially-formatted ASCII file creates one entire spreadsheet.
This spreadsheet file is in its final, directly-readable 
binary form, and need undergo no further translation before being
read directly by Lotus/123(TM).  This feature makes
.I conv123
especially
suited to networked environments, where report writers or programs
running on transaction processing machines can easily and directly
create spreadsheet output, instantly ready for copying, downloading,
or direct access in the MS-DOS(TM) or Macintosh(TM) environments.
In addition to individual data values, certain common aggregate
functions (formulas) may be specified (sum,min,max,variance,std
deviation, and average), and five binary functions (addition,
subtraction, multiplication, division, and percentage) may also
be specified.

The input file format has been optimized for generation by a
simple top-down, left-to-right report generator,
and is human-readable.  The input logic has been made fairly intelligent,
automatically discerning numeric and non-numeric values,
calculating formula ranges, and calculating the maximum number
of decimal places in a row or column.

.SH OPTIONS
There are currently no command-line options for
.I conv123.


.SH PARAMETERS
The only parameters are the input and output file names, both
of which are optional;
.I conv123
can act as a filter, and by
default accepts its input from stdin and sends its output to
stdout.

.SH INPUT FILE FORMAT
The major feature of the input file format is a "definition"
line, consisting of all dashes (-) and spaces.  Each string of
contiguous dashes on the "definition" line will become a column
in Lotus/123(TM), the width of that column being determined by the
number of dashes.  Characters in the input file that fall outside
these "column definitions" are ignored.  In addition, up to two
lines above the "definition" line can be used for header information,
though this is not required.  There may not be any white space before
the first dash on the definition line; in other words, the definition
line MUST begin in the first column of the input file.

Following is an example input file, followed by a more detailed
explanation of its contents:

                   Jan          Feb         Mar         Qtr
     Product      Sales        Sales       Sales       Total
    ---------  -----------  ----------- ----------- -----------
    widget-1     45000.87       54390.9    75947.12     @RSUM
    widget-2     47000.78       64550.0    73647.77     @RSUM
    widget-3     47800.11       67830.5    54536.37     @RSUM
    widget-4     33003.30       76530.7    98653.83     @RSUM
    widget-5     12000.44       64836.7    27637.91     @RSUM

                     @CSUM        @CSUM       @CSUM     @CSUM

This is a fairly simple spreadsheet with ten rows and five
columns.  The fifth column is actually made of formulas
("row sums"), as is the last row ("column sums").
.I conv123
is
smart enough to know that each data value in the first column
is non-numeric, and creates the formulas in the fifth column
accordingly, i.e., the "row sums" are sums of the values in
the second through the fourth columns.  In addition, since
.I conv123
keeps track of the maximum number of decimal places
used in a given row or column, the number of decimal places
displayed in Lotus/123(TM) will be two (2) for all formulas,
with the exception of the CSUM in the third column, which will
display with one (1) decimal place.  Note also that the position
of the number or formula code in each column is not significant,
as long as the entire values fits under the dashed-line column
definition(s).  Also, most spreadsheet programs will re-position
the values within each cell as the values are displayed.


.SH FORMULA CODES
As has been shown in the above example, cell contents may be
decimal numbers or character values, with the width of the cell
determined by the width of the contiguous dashes in the "definition"
line.  Cell contents may also be formula codes.  A formula code is
a contiguous set of characters, always beginning with a commercial
"at" sign (@).  There are three types: row, column and both.  The
type indicator is a single upper case letter, one of "R" (row), "C"
(column) or "B" (both), immediately following the "@" ("at" sign).
The formula function immediately follows the type indicator, and may
be an aggregate function or a binary function.  The range for an
aggregate function is calculated automatically, being the range from
the earliest available cell -- i.e., the cell farthest "above" for
columns, or farthest to the "left" for rows -- to the cell immediately
adjacent to the formula definition.  Note that this "farthest" cell
location is re-adjusted when a non-numeric cell definition is 
encountered, being set to begin at the next numeric cell (this does
not include non-cells, i.e., positions within the spreadsheet with
no value).  The automatic range start position may be overridden for
"row" functions by placing a numeric value at the end of the formula
code (still contiguous).  A numeric value at the end of a formula code
is treated as an absolute column number, with the first column
defined as column one (1).

The range for a binary function is simply the last two adjacent
cells "above" the formula definition for column functions, and
"left" of the formula definition for row functions.  The adjacent
cell is always used in the calculation; an override number after
the function code in a row function is treated as an absolute
column number overriding the penultimate cell.  For commutative
functions (addition and multiplication), order is not significant.
For non-commutative functions (subtraction and division), the
definition is a - b, or a/b, where a is the penultimate cell, and b
is the immediately adjacent cell.  If there is an override column
number, THIS ORDER IS REVERSED -- i.e., b - a and b/a, where a is
the overriding column number and b is the immediately adjacent cell.
This situation is made clearer in the examples below.

Following is a summary of the available functions, which MUST be in
upper case where applicable:

    SUM -- Performs addition on a range of cells.
    AVG -- Takes the average value on a range of cells.
    MIN -- Extracts the minimum value on a range of cells.
    MAX -- Extracts the minimum value on a range of cells.
    VAR -- Calculates the variance on a range of cells.
    STD -- Calculates the standard deviation on a range of cells.
     +  -- Adds two cells.
     -  -- Subtracts two cells.
     *  -- Multiplies two cells.
     /  -- Divides two cells.
     %  -- Divides two cells, and multiplies the result by 100.
    NUL -- This is actually not a function, but performs the task
           of a "placeholder", preventing the "resetting" of auto-
           matically calculated numeric ranges without creating
           cell contents.  This function is most frequently used to
           extend the row range of aggregate column functions.


.SH FORMULA EXAMPLES
    Here are some sample input lines, with result following.

    input:
    ---------  -----------  ----------- ----------- -----------
    widget-1      100.00       200.00      300.00     @RSUM

    result:
    ---------  -----------  ----------- ----------- -----------
    widget-1        100.00       200.00      300.00      600.00


    input:
    ---------  -----------  ----------- ----------- -----------
    widget-1       100.00       200.00      300.00     @R+

    result:
    ---------  -----------  ----------- ----------- -----------
    widget-1        100.00       200.00      300.00      500.00


    input:
    ---------  -----------  ----------- ----------- -----------
    widget-1       100.00       200.00      300.00     @R/

    result:
    ---------  -----------  ----------- ----------- -----------
    widget-1        100.00       200.00      300.00        0.67


    input:
    ---------  -----------  ----------- ----------- -----------
    widget-1       100.00       200.00      300.00     @R/3

    result:
    ---------  -----------  ----------- ----------- -----------
    widget-1        100.00       200.00      300.00        1.50

    Note in the above example that the order of divisor/dividend as
    regards the immediately adjacent cell is the reverse of the
    previous example.


 input:
 ---------  --------  -------- -------- ---------- ------- ----------
 widget-1    100.00    200.00     @RSUM     300.00   400.0     @RSUM
 widget-2    100.00    100.00     @RSUM     100.00   200.0     @RSUM4
                                   @NUL
              @CSUM     @BSUM     @CSUM

 result:
 ---------  --------  -------- -------- ----------- ------- ----------
 widget-1     100.00    200.00   300.00      300.00   400.0     700.00
 widget-2     100.00    100.00   200.00      100.00   200.0     500.00

              200.00    500.00     0.00



.SH SPECIAL CONSIDERATIONS

The dashed-line definitions must begin at the first character
of the line in which they occur.

Output files are basically in Lotus(TM) ".WK1" format, and output
files should be named with this extension to provide maximum
information to spreadsheet programs.

Binary functions in column formula codes are allowed, but rarely
make much sense.

The program is simple by design, and makes no attempt to incorporate
useful, more complex features such as absolute cell positioning,
general-purpose formulae, macros, etc.  It also has not undergone
a complete test suite.


.SH BUGS

The maximum number of columns is currently hard-coded into the 
program at 124; memory should be dynamically allocated to support
an arbitrary spreadsheet size.

Optimizations could be applied at several points.

Starting/ending points for formulae should be arbitrarily
specifiable.
