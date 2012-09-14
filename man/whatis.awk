
# awk script, called by make_whatis.csh
#   BUG:   it produces <A HREF="xyopen,.y.html">

# <a href=bimaudoc/index.html><img src=/icons/blue.ball.gif alt="*"></a>
#  <a href=bimaudoc/index.html>Bima User Documentation</a>

{
   # printf("<li>");
   cat = substr($2,2,1);        # check if this is ok
   name = $1;
   # proper man_html 
   printf("<LI><A HREF=\"%s.%s.html\"> %s </A> ",name,cat,name);
   # old cat
   # printf("<LI><A HREF=\"man/cat%s/%s.%s\"> %s </A> ",cat,name,cat,name);
   printf("%s\n",$0);
}

