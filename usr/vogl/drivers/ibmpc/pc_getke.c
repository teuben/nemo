int
pc_getkey()
{
	return((int)getch());
}

int
pc_checkkey()
{
	if (kbhit())
		return((int)getch());

	return(0);
}
