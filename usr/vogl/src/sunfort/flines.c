/*
 * setlin_
 */
setlin_(n)
	int	*n;
{
	setlinestyle(*n);
}

/*
 * deflin_
 */
deflin_(n, s)
	int	*n, *s;
{
	deflinestyle(*n, *s);
}
/*
 * setlinestyle_
 */
setlinestyle_(n)
	int	*n;
{
	setlinestyle(*n);
}
/*
 * deflinestyle_
 */
deflinestyle_(n, s)
	int	*n, *s;
{
	deflinestyle(*n, *s);
}

/*
 * linewi_
 */
linewi_(n)
	int	*n;
{
	linewidth(*n);
}

/*
 * linewidth_
 */
linewidth_(n)
	int	*n;
{
	linewidth(*n);
}
