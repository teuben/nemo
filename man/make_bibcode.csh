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

echo "These are the BIBCODE references within NEMO that refer to papers, so we can"
echo "establish a link between code and paper. For most codes we do"
echo "this via their unix man page (.1, and even .5). In these man pages"
echo "there should be a reference in the FILES section where the code"
echo "is located in NEMO.   For some codes that do not have a manual page,"
echo "the source code file will be mentioned."
echo "<P>"
echo "Created on: `date` by `whoami`@`hostname``domainname`<P>"
echo "<\!-- Created by $0 - do not edit -->"
echo "<DIR COMPACT>"

# @todo need the comment field
grep ^@ads man?/*.? | awk '{print $1,$2}' | awk -F/ '{print $1,$2}' | sed s/:@ads// | awk -f bibcode.awk

# @todo: scan all source code files in usr and src for ^@ads

echo "</DIR></BODY></HTML>"

