#include <stdio.h>

main(argc, argv)
int argc;
char *argv[];
{
    char string[30];
    int i,j;
    real y=0.9;
    
    init_X(argc, argv);

    x_space(0.0,0.0,1.0,1.0);

    x_move(0.2,y);

    x_font_by_name("-jis-fixed-medium-r-normal--0-0-75-75-c-0-jisx0208.1983-0");

    for(i=32; i<256; i+=16)
    {
	sprintf(string,"%3d:   ",i);
	for(j=i; j<i+16; j++)
	    string[6+(j%16)] = j;
	x_label(string);
	y = y-0.03;
	x_move(0.2,y);
    }

    x_doplot();
    x_buttonwait();

    close_X();
}

