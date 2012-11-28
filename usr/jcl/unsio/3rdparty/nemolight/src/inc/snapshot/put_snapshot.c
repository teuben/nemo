/*
 * PUT_SNAPSHOT.C: write snapshots in their striped form.
 *		Although the data is now not in the form
 *		of the Body structure, it has the advantage
 *		that no intermediate arrays are needed.
 */

put_snapshot(instr, ssptr)
stream instr;			/* input stream, of course */
SSptr ssptr;			/* pointer to body array */
{
    real time, *rbuf;
    int nbody, bits;
    
    nbody = ssptr->nbody;
    bits = ssptr->bits;
    time = ssptr->time;

    put_set(instr, SnapShotTag);            /* open  the snapshot */

    put_set(instr, ParametersTag);
    put_data(instr, NobjTag, IntType, &nbody, 0);
    put_data(instr, TimeTag, RealType, &time, 0);
    put_tes(instr, ParametersTag);

    put_set(instr, ParticlesTag);
    /*  COORDSYS !! */
    rbuf = ssptr->mass;
    if (rbuf == NULL && bits&MassBit)
        error("put_snapshot: no masses supplied");
    else
        put_data(instr, MassTag, RealType, rbuf, nbody, 0);
    rbuf = ssptr->phase;
    if (rbuf == NULL && bits&PhaseSpaceBit)
        error("put_snapshot: no phases supplied");
    put_data(instr, PhaseSpaceTag, RealType, rbuf, nbody, 2, NDIM, 0);
    put_tes(instr, ParticlesTag);

    put_tes(instr, SnapShotTag);
}
