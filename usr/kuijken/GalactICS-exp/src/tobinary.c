#include <stdio.h>

typedef struct {
	float r[6];
} phase;

main()
{
	int i, nobj;
	float t, *m;
	phase *r;

	fscanf(stdin,"%d %f\n",&nobj,&t);

	m = (float *) calloc(nobj,sizeof(float));
	r = (phase *) calloc(nobj,sizeof(phase));

	for(i=0; i<nobj; i++) {
		fscanf(stdin,"%f %f %f %f %f %f %f\n",
			m+i, &r[i].r[0],&r[i].r[1],&r[i].r[2],
				&r[i].r[3],&r[i].r[4],&r[i].r[5]);
	}

	fwrite(&nobj,sizeof(int),1,stdout);
	fwrite(m,nobj*sizeof(float),1,stdout);
	fwrite(&t,sizeof(float),1,stdout);
	fwrite(r,nobj*sizeof(phase),1,stdout);
}
