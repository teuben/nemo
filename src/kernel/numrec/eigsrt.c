void eigsrt(d,v,n)
float d[],**v;
int n;
/*205401*/
{
	int k,j,i;
	float p;

	for (i=1;i<n;i++) {
		k=i;
		p=d[i];
		for (j=i+1;j<=n;j++) {
			if (d[j] >= p) {
				k=j;
				p=d[j];
			}
		}
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}
