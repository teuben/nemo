

# from: @ads 2005A&A...433...57T
#   to: 
# stem: https://ui.adsabs.harvard.edu/abs/

# @todo  s/\.1NEMO/\.1/g   and also for 3,5,6,8,l


BEGIN {
    url = "https://ui.adsabs.harvard.edu/abs"
}

{
    if ($1 == "@ads") {
	bibcode = $2;
	printf("@ads <A HREF=%s/%s> %s </A>\n", url, bibcode, bibcode);
    } else if ($1 == "<body") {
	print $0;
	printf("HTML automatically generated with <A HREF=http://manpages.ubuntu.com/manpages/bionic/man1/rman.1.html>rman</A><br>\n");
	
    } else
	print $0;
}

