/* testing symbols */

int i_glob;
static int i_loc;

float f_glob=1.0;
static float f_loc=2.0;

extern double sqrt(double);

int i_fun(float x)
{
	return(sqrt(x) + 21.0);
}


