/*
 * minv4
 *
 *	find the inverse of the 4 by 4 matrix b using gausian elimination
 * and return it in a.
 * 
 *	(We don't actually use this yet in VOGL - maybe one day).
 */
minv4(a, b)
	Matrix	a, b;
{
	float	val, val2;
	int	i, j, k, ind;
	Matrix	tmp;

	identmatrix(a);

	copymatrix(tmp, b);

	for (i = 0; i != 4; i++) {

		val = tmp[i][i];		/* find pivot */
		ind = i;
		for (j = i + 1; j != 4; j++) {
			if (fabs(tmp[j][i]) > fabs(val)) {
				ind = j;
				val = tmp[j][i];
			}
		}

		if (ind != i) {			/* swap columns */
			for (j = 0; j != 4; j++) {
				val2 = a[i][j];
				a[i][j] = a[ind][j];
				a[ind][j] = val2;
				val2 = tmp[i][j];
				tmp[i][j] = tmp[ind][j];
				tmp[ind][j] = val2;
			}
		}

		if (val == 0.0)
			fatal("art: singular matrix in minv4.\n");

		for (j = 0; j != 4; j++) {
			tmp[i][j] /= val;
			a[i][j] /= val;
		}

		for (j = 0; j != 4; j++) {	/* eliminate column */
			if (j == i)
				continue;
			val = tmp[j][i];
			for (k = 0; k != 4; k++) {
				tmp[j][k] -= tmp[i][k] * val;
				a[j][k] -= a[i][k] * val;
			}
		}
	}
}
