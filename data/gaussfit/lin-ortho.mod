parameter a, b;
observation x,y;

main()
{
	while(import())
		export(y - (a + b*x));
}


