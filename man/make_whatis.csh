#! /bin/csh -f
# Generate a 'Whatis' html file
#

echo "<HTML>"
echo "<HEAD>"
echo "<TITLE>NEMO WhatIs Page </TITLE>"
echo "</HEAD></BODY>"

echo "<CENTER>"
echo -n "[<A HREF=index1.html>man1</A>,"
echo -n "<A HREF=index3.html>man3</A>,"
echo -n "<A HREF=index5.html>man5</A>,"
echo -n "<A HREF=index6.html>man6</A>,"
echo -n "<A HREF=index8.html>man8</A>,"
echo -n "<A HREF=indexl.html>manl</A>".
echo -n "<A HREF=whatis.html>whatis</A>]".
echo "</CENTER>"

echo "<H1>NEMO Whatis Manual Reference  </H1>"
echo "Created on: `date` by `whoami`@`hostname``domainname`<P>"
echo "<\!-- Created by $0 - do not edit -->"
echo "Note: this page may still have some bad cross-references"
echo "<DIR COMPACT>"
awk -f whatis.awk whatis
echo "</DIR></BODY></HTML>"

