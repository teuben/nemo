#include <stdio.h>

main()
{
	int	c1, c2, c3;

	while ((c1 = getchar()) != EOF) {
		if (c1 == '_') {
			if ((c2 = getchar()) == 'c') {
				if ((c3 = getchar()) == '_') {
					putchar(c3);
				} else {
					putchar(c1);
					putchar(c2);
					putchar(c3);
				}
			} else {
				putchar(c1);
				putchar(c2);
			}
		} else {
			putchar(c1);
		}
	}
}
			
