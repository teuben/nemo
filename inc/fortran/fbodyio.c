/*
 * FIO: example of fortran I/O routines
 *	It assumes that the <snapshot/body...> stuff has been included already
 *	in the source code
 *	Example: see testfio.c ~src/pjt/nbody
 *	Writing means filling in the fortran arrays (export)
 *	Reading means getting data from fortran arrays to C-structures.
 *
 *	Peter Teuben	--	November 1988
 */
 
#define READ_FLAG  0
#define WRITE_FLAG 1

fio (wflag, nbody, btab, mass, pos, vel)
int wflag;
int nbody;
Body *btab;
float *mass, *pos, *vel;
{
    Body *bi;
    int  i;

    for (i=0, bi=btab; i<nbody; i++, bi++) {
        if(wflag) {
            *mass = Mass(bi);
            SETV(pos,Pos(bi));
            SETV(vel,Vel(bi));
        } else {
            Mass(bi) = *mass;
            SETV(Pos(bi),pos);
            SETV(Vel(bi),vel);
        }
        mass++;
        pos += NDIM;
        vel += NDIM;
    }
}

            

