.TH CATPS 1NEMO "24 June 1997"
.SH NAME
catps, catpgps - catenate postscript files
.SH SYNOPSIS
catps file1 file2 .. fileN 
.SH DESCRIPTION
Catenates postscript files. You probably should use the version which 
was written for specific type of postscript files, e.g. 
\fIcatps\fP, \fIcatpgps\fP for PGPLOT based postscript etc.
.PP
The \fIpgmerge(1)\fP tool from the psutils package can also
merge multiple postscript files, but as their documentation
says, they are for a very specific case, it does not merge
files in the general case.
.SH SEE ALSO
yapp(3NEMO), layout(1NEMO), psmerge(1)
.SH HISTORY
.nf
.ta +1i +4i
24-jun-97	added catpgps	PJT
.fi
