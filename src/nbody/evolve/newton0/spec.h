/* spec.h - Celldivisionmethod, Diagnostics, Dynparam, Integrationscheme,
            Methods, N2bcalc, Nbccalc, Nfcalc, Nstep_de,
            Number_of_dynparam, Number_of_methods, Number_of_specidiag, 
            Number_of_specrdiag, Softfocus, SPECmin_pair, 
            Softparam, Specidiag, Specrdiag, Stepparam, Timestepmethod, 
            Tolsqparam, spec, specptr */

/*
 *  spec.h: for "spec", a substructure of a "state".
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
/*-----------------------------------------------------------------------------
 *  spec, specptr  --  contains parameters which control the details of the
 *                     orbit integration and the computation of those
 *                     diagnostics which are intimately interwoven with the
 *                     force calculation process; 
 *                     for global control parameters see  diag.h ;
 *                     for global diagnostics see  diag.h .
 *                     note: recipe for adding an extra element:
 *                           1) increase the appropriate array length by
 *                              adding 1 to Number_of_... ;
 *                           2) introduce a macro below to provide a handle
 *                              to access the new element;
 *                           3) provide a short description at the end of this
 *                              file to document the meaning of the new element
 *                           in addition, if the new element should be read
 *                           from the command line (which is often the case)
 *                           then:
 *                           4) add an entry to the  defv[]  string at the top
 *                              of the file  newton0.c;
 *                           5) add a line to the procedure  set_specs()
 *                              in the file  newton0.c , using one of the 
 *                              getparam tools.
 *                           no change has to be made to the procedures in
 *                           the file  save.c , unless a new type of array is
 *                           added to  spec .
 *-----------------------------------------------------------------------------
 */
#ifndef TREE
#  define  Number_of_methods	  4
#  define  Number_of_dynparam	  2
#  define  Number_of_specidiag    1
#endif
#ifdef TREE
#  define  Number_of_methods	  5
#  define  Number_of_dynparam	  3
#  define  Number_of_specidiag    4
#endif
#  define  Number_of_specrdiag    1

typedef struct
    {
    string  methods[  Number_of_methods   ];
    real   dynparam[  Number_of_dynparam  ];
    int    specidiag[ Number_of_specidiag ];
    real   specrdiag[ Number_of_specrdiag ];
    } spec, *specptr;

/*-----------------------------------------------------------------------------
 *  macros to extract arrays of components from a spec structure:
 *  NOTE:  ptr should be a pointer of type "specptr":
 *-----------------------------------------------------------------------------
 */
#define  Methods(ptr)        ((ptr)->methods)            /* type: stringptr  */
#define  Dynparam(ptr)       ((ptr)->dynparam)           /* type: realptr    */
#define  Specidiag(ptr)      ((ptr)->specidiag)          /* type: intptr     */
#define  Specrdiag(ptr)      ((ptr)->specrdiag)          /* type: realptr    */

/*-----------------------------------------------------------------------------
 *  macros to extract individual components from a spec structure (see next
 *  page for a description of their meaning):
 *  NOTE:  ptr should be a pointer of type "specptr"
 *-----------------------------------------------------------------------------
 */
#define  Softfocus(ptr)             ((ptr)->methods[0])    /* type: string   */
#define  Timestepmethod(ptr)        ((ptr)->methods[1])    /* type: string   */
#define  Integrationscheme(ptr)     ((ptr)->methods[2])    /* type: string   */
#define  Diagnostics(ptr)           ((ptr)->methods[3])    /* type: string   */
#ifdef TREE
#  define  Celldivisionmethod(ptr)  ((ptr)->methods[4])    /* type: string   */
#endif

#define  Stepparam(ptr)             ((ptr)->dynparam[0])   /* type: real     */
#define  Softparam(ptr)             ((ptr)->dynparam[1])   /* type: real     */
#ifdef TREE
#  define  Tolsqparam(ptr)          ((ptr)->dynparam[2])   /* type: real     */

#  define  Nfcalc(ptr)              ((ptr)->specidiag[1])  /* type: int      */
#  define  N2bcalc(ptr)             ((ptr)->specidiag[2])  /* type: int      */
#  define  Nbccalc(ptr)             ((ptr)->specidiag[3])  /* type: int      */
#endif
#define  Nstep_de(ptr)              ((ptr)->specidiag[0])  /* type: int      */

#define  SPECmin_pair(ptr)          ((ptr)->specrdiag[0])  /* type: real     */

/*-----------------------------------------------------------------------------
 *  short descriptions of the various components of a "spec" structure:
 *
 *      Softfocus()          choice of type of softening.
 *      Timestepmethod()     choice of timestep criterion for integration.
 *      Integrationscheme()  choice of integration scheme.
 *      Diagnostics()        choice of set of diagnostics.
 *      Celldivisionmethod() choice of cell division criterion.
 *      Stepparam()          dimensionless integration accuracy parameter,
 *			     used in determining the size of the time steps
 *      Softparam()          potential softening length
 *      Tolsqparam()         accuracy parameter for opening cells: 0.0 => exact
 *                           Tolsq is the square of the opening angle theta.
 *      Nfcalc()             number of n-on-1 force calculations  
 *      N2bcalc()	     number of 2-body force calculations  
 *      Nbccalc()	     number of body-cell force calculations
 *      Nstep_de()           number of steps after which the maximum energy
 *                           error is checked (cf.  DIAGemax_err  in  diag.h;
 *                           this DIAGemax_err to be implemented one day)
 *      SPECmin_pair()       diagnostic variable keeping registering the
 *                           minimum pair distance occurrence. NOTE: if
 *                           softening is used, the softened distance is given
 *                           as returned by the softening routines in  soften.c
 *                           as an effective distance >= the softening distance
 *-----------------------------------------------------------------------------
 */

/* endof: spec.h */
