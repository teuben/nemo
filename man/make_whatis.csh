#! /bin/csh -f
# Generate a 'Whatis' html file
#

echo "<HTML>"
echo "<HEAD>"
echo "<TITLE>NEMO WhatIs Page </TITLE>"
echo "</HEAD></BODY>"
echo "<H1>NEMO Whatis Manual Reference  </H1>"
echo "Created on: `date` by `whoami`@`hostname``domainname`<P>"
echo "<\!-- Created by $0 - do not edit -->"
echo "Note: this page may still have some bad cross-references"
echo "<DIR COMPACT>"
awk -f whatis.awk man/whatis
echo "</DIR></BODY></HTML>"

