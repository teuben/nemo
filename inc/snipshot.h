/*
 * SNIPSHOT.H: definitions and accessor macros for a simple internal
 * representation of an N-body system.  Snip-Snap was a Duch comedy
 * revue, hence the name snipshot to complement snapshot.
 */

/*
 * SNIPBODY structure.
 */

typedef struct {
    real snipmass;
    real snipphase[2][NDIM];
} snipbody, *snipbodyptr;

#define SnipMass(p)	((p)->snipmass)
#define SnipPos(p)	((p)->snipphase[0])
#define SnipVel(p)	((p)->snipphase[1])
#define SnipPhase(p)	((real *) (p)->snipphase)

/*
 * SNIPSHOT structure.
 */

typedef struct {
    real sniptime;
    int snipnobj;
    snipbodyptr snipbodies;
} snipshot, *snipptr;

#define SnipTime(s)	((s)->sniptime)
#define SnipNobj(s)	((s)->snipnobj)
#define SnipBodies(s)	((s)->snipbodies)

/*
 * Functions for manipulating snip-shots.
 */

snipptr snap2snip(/* itemptr */);
itemptr snip2snap(/* snipptr */);

snipptr allocsnip(/* int */);
snipptr copysnip(/* snipptr */);

void freesnip(/* snipptr */);

void zerocm4snip(/* snipptr */);

