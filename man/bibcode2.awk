

# from: @ads 2005A&A...433...57T
#   to: 
# stem: https://ui.adsabs.harvard.edu/abs/

# @todo  s/\.1NEMO/\.1/g   and also for 3,5,6,8,l


BEGIN {
    url1 = "https://ui.adsabs.harvard.edu/abs"
    url2 = "http://manpages.ubuntu.com/manpages/bionic/man1/rman.1.html"
    url3 = "https://astronemo.readthedocs.io"
}

{
    if ($1 == "@ads") {
	bibcode = $2;
	printf("@ads <A HREF=%s/%s> %s </A>\n", url1, bibcode, bibcode);
    } else if ($1 == "<body") {
	print $0;
	printf("This HTML automatically generated with <A HREF=%s>rman</A> for <A HREF=%s>NEMO</A><br>\n",url2,url3);
	
    } else
	print $0;
}

