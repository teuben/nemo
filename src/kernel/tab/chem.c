/*
 *  Mol Weights
 */

#include <stdinc.h>

typedef struct {
  char *atom;
  int number;
  int weight;
} atom;

static atom a[] = {
    "H",    1,  1,
    "D",    1,  2,
    "He",   2,  4,
    "Li",   3,  7,
    "Be",   4,  9,
    "B",    5,  11,
    "C",    6,  12,
    "N",    7,  14,
    "O",    8,  16,
    "F",    9,  19,
    "Ne",   10, 20,
    "Na",   11, 23,
    "Mg",   12, 24,
    "Al",   13, 27,
    "Si",   14, 28,
    "P",    15, 31,
    "S",    16, 32,
    "Cl",   17, 35, /* ? */
    "Ar",   18, 40,
    "K",    19, 39,
    "Ca",   20, 40,
    "Sc",   21, 45,
    "Ti",   22, 48,
    "V",    23, 51, 
    "Cr",   24, 52,
    "Mn",   25, 55,
    "Fe",   26, 56,
    NULL,   0,  0,  NULL,
};

/*
 *  return atomic number (number of protons)
 */

int at_number(name)
char *name;
{
    atom *p;

    for (p=a; *p->atom; p++)
        if (streq(name,p->atom)) return p->number;
    return 0;
}

