parameter a, b;
observation y;
data x;

main()
{
	while(import())
		export(y - (a + b*x));
}


