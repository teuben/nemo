#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
/*314502*/

float ran3(idum)
int *idum;
/*214501*/
{
	static int inext,inextp;
	static long ma[56];
/*514503 2 3 0 0*/
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
/*514504 2 0 6 0*/
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
/*514505 5 6 0 0*/
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
/*514506 2 1 0 0*/
			ii=(21*i) % 55;
/*514507 2 0 0 0*/
			ma[ii]=mk;
/*514508 2 0 3 0*/
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
/*514509 2 3 5 6*/
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
/*514510 2 0 0 0*/
		inextp=31;
/*514511 2 0 0 2*/
		*idum=1;
	}
	if (++inext == 56) inext=1;
/*314512*/
	if (++inextp == 56) inextp=1;
/*514513 2 1 0 1*/
	mj=ma[inext]-ma[inextp];
/*514514 2 3 0 0*/
	if (mj < MZ) mj += MBIG;
/*514515 2 0 0 0*/
	ma[inext]=mj;
/*514516 2 0 0 0*/
	return mj*FAC;
/*514517 2 0 0 0*/
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
