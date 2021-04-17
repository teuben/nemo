

# from: @ads 2005A&A...433...57T
#   to: 
# stem: https://ui.adsabs.harvard.edu/abs/

BEGIN {
    url = "https://ui.adsabs.harvard.edu/abs"
}

{
    if ($1 == "@ads") {
	bibcode = $2;
	printf("@ads <A HREF=%s/%s> %s </A>\n", url, bibcode, bibcode);
    } else
	print $0;
}

