#! /bin/csh -f
# Generate a 'bibcode' html file
#

echo "<HTML>"
echo "<HEAD>"
echo "<TITLE>NEMO BIBCODE Page </TITLE>"
echo "</HEAD></BODY>"

echo "<CENTER>"
echo -n "[<A HREF=index1.html>man1</A>,"
echo -n "<A HREF=index3.html>man3</A>,"
echo -n "<A HREF=index5.html>man5</A>,"
echo -n "<A HREF=index6.html>man6</A>,"
echo -n "<A HREF=index8.html>man8</A>,"
echo -n "<A HREF=indexl.html>manl</A>,"
echo -n "<A HREF=whatis.html>whatis</A>,"
echo -n "<A HREF=bibcode.html>bibcode</A>]"
echo "</CENTER>"

echo "<H1>NEMO BIBCODE Manual Reference  </H1>"
echo "Created on: `date` by `whoami`@`hostname``domainname`<P>"
echo "<\!-- Created by $0 - do not edit -->"
echo "<DIR COMPACT>"
grep ^@ads man?/*.? | awk '{print $1,$2}' | awk -F/ '{print $1,$2}' | sed s/:@ads// | awk -f bibcode.awk

echo "</DIR></BODY></HTML>"

