/*
 *  POTENTIAL: 
 */

#ifndef _potential_h
#define _potential_h

typedef struct {
    string name;
    string pars;
    string file;
    struct a_potential *next;
} a_potential;

proc get_potential    (string, string, string);
proc get_inipotential (void);
real get_pattern      (void);

#endif
