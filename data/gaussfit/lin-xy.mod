parameter a, b;
observation x;
data y;

main()
{
	while(import())
		export(y - (a + b*x));
}


