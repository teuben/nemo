/*
 * MODULE.H: data structure with slots for major operations.
 */

typedef struct {
    string  modname;
    void    (*modinit)(/* itemptr, stream */);
    itemptr (*modwork)(/* itemptr, itemptr, stream */);
    void    (*modsave)(/* stream, stream */);
    void    (*modrestore)(/* stream, stream */);
} module, *modptr;

/*
 * Standard accessor macros used to extract the name 
 * and/or invoke the functions of a module.
 */

#define ModName(mp)       ((mp)->modname)
#define ModInit(mp)       (* (mp)->modinit)
#define ModWork(mp)	  (* (mp)->modwork)
#define ModSave(mp)       (* (mp)->modsave)
#define ModRestore(mp)    (* (mp)->modrestore)
