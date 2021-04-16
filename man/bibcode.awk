
# awk script, called by make_bibcode.csh
#   BUG:   it produces <A HREF="xyopen,.y.html">

# <a href=bimaudoc/index.html><img src=/icons/blue.ball.gif alt="*"></a>
#  <a href=bimaudoc/index.html>Bima User Documentation</a>


# from: man1 CGS.1 2005A&A...433...57T
#   to: 
# stem: https://ui.adsabs.harvard.edu/abs/

BEGIN {
    url = "https://ui.adsabs.harvard.edu/abs"
    printf("<TABLE BORDER>\n");
    printf("<td>   <B> NAME </B> </td>\n")
    printf("<td>   <B> BIBCODE </B> </td>\n")
    printf("<td>   <B> COMMENTS </B> </td>\n")
}

{
    manl = $1;
    name = $2;
    bibcode = $3;

    printf("<tr>");
    printf("<td> <A HREF=\"%s.html\"> %s </A> </td>\n",name,name);
    printf("<td> <A HREF=\"%s/%s\">   %s </A> </td>\n",url,bibcode,bibcode);
    printf("</tr>\n");
}

END {
    printf("</TABLE>\n");
}
