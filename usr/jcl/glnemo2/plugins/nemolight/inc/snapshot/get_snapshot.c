/*
 * GET_SNAPSHOT.C: read in snapshots in their striped form.
 *		Although the data is now not ijn the form
 *		of the Body structure, it has the advantage
 *		that no intermediate arrays are needed.
 */


/* Similar to a body, we define a whole snapshot with pointers to
 * the homogeneous arrays
 */

typedef struct _snapshot {
    int nbody;
    int bits;           /* what's all present in here */
    real time;
    real *mass;        
    real *phase;
} SS, *SSptr;



get_snapshot(instr, ssptr)
stream instr;			/* input stream, of course */
SSptr ssptr;			/* pointer to body array */
{
    real time, *rbuf;
    int nbody, bits = 0;
    char *realloc();

    if (get_tag_ok(instr, SnapShotTag)) {
	get_set(instr, SnapShotTag);            /* open  the snapshot */

        if (get_tag_ok(instr, ParametersTag)) {
            get_set(instr, ParametersTag);
            get_data(instr, NobjTag, IntType, &nbody, 0);
            ssptr->nbody = nbody;
            if (get_tag_ok(instr, TimeTag)) {       /* time data specified? */
                get_data_coerced(instr, TimeTag, RealType, &time, 0);
                bits |= TimeBit;                  /*   set time bit */
                ssptr->time = time;
            }
            get_tes(instr, ParametersTag);
        } else
            error("get_snapshot: No parameters...");

        if (get_tag_ok(instr, ParticlesTag)) {
            get_set(instr, ParticlesTag);
            rbuf = ssptr->mass;
            if (get_tag_ok(instr, MassTag)) {
                if (rbuf==NULL)
                  rbuf = (real *) malloc(nbody * sizeof(real));
                else
                  rbuf = (real *) realloc(rbuf, nbody * sizeof(real));
                if (rbuf == NULL)
                    error("get_snapshot: not enuf memory for %d masses\n",nbody);
	        get_data_coerced(instr, MassTag, RealType, rbuf, 
                        nbody, 0);
	        bits |= MassBit;
                ssptr->mass = rbuf;
            }
            rbuf = ssptr->phase;
            if (get_tag_ok(instr, PhaseSpaceTag)) {
                if (rbuf==NULL)
                  rbuf = (real *) malloc(NDIM*2*nbody * sizeof(real));
                else
                  rbuf = (real *) realloc(rbuf, NDIM*2*nbody * sizeof(real));
                if (rbuf == NULL)
                    error("get_snapshot: not enuf memory for %d phases\n",nbody);
	        get_data_coerced(instr, PhaseSpaceTag, RealType, rbuf, 
                        nbody, 2, NDIM, 0);
	        bits |= PhaseSpaceBit;
                ssptr->phase = rbuf;
            }
            get_tes(instr, ParticlesTag);
        } else
            error("get_snapshot: No particles...");

        get_tes(instr, SnapShotTag);
	ssptr->bits = bits;
        return(1);
    } else
        return(0);
}
