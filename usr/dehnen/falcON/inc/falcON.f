C -*- FORTRAN -*-                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C  falcON.f                                                                    |
C                                                                              |
C==============================================================================+
C                                                                              |
C  faclON = Force ALgorithm with Complexity O(N)                               |
C                                                                              |
C==============================================================================+
C                                                                              |
C  header file for FORTRAN users                                               |
C  (C++ and C users, see files falcON.h and falcON_C.h, respectively)          |
C                                                                              |
C  Copyright Walter Dehnen, 2000-2003                                          |
C  e-mail:   wdehnen@aip.de                                                    |
C  address:  Astrophysikalisches Institut Potsdam,                             |
C            An der Sternwarte 16, D-14482 Potsdam, Germany                    |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C  FORTRAN routines implementing the code described by Dehnen (2000,2002).     |
C  These FORTRAN routines call C routines declared in file falcON_C.h, which   |
C  in turn call the original C++ functions declared in file falcON.h.          |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C  See file readme.falcON.html for some guidelines on how to install, #include |
C  and link this file and the code declared in it. You have the choice between |
C  various options determining some properties of the code.                    |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C  CONTENTS                                                                    |
C  ========                                                                    |
C                                                                              |
C  0  Preliminaries                                                            |
C     0.1  Dimensionality                                                      |
C     0.2  Data type for I/O with the tree code                                |
C     0.3  FORTRAN specific problems                                           |
C  1  Initialisation and clearing up                                           |
C     1.1  Before using the code                                               |
C     1.2  After using the code                                                |
C  2  Meaning of the body flags                                                |
C  3  Force approximation and related                                          |
C     3.1  Generating a tree structure                                         |
C     3.2  Approximating the forces                                            |
C     3.3  Estimating mass- and number-density                                 |
C  4  Sticky-particle and SPH support                                          |
C     4.1  Sticky-particle support                                             |
C     4.2  SPH support                                                         |
C  5  Other features                                                           |
C  6  Known bugs and problems                                                  |
C     6.1  Test-bodies are not possible                                        |
C     6.2  Bodies at identical positions                                       |
C  References                                                                  |
C  A  Routines used for test purposes only                                     |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  0 PRELIMINARIES                                                             |
C  ===============                                                             |
C                                                                              |
C  0.1 Dimensionaliy                                                           |
C  -----------------                                                           |
C  The code allows for 2 or 3 (default) dimensions, controlled by the macro    |
C  "TWODIMENSIONAL", which is set as compiler option, see file readme.falcON.  |
C  Below, I'll explain the 3D versions. The only difference to the 2D versions |
C  is that arguments specifying the 3rd (=z) dimension of vectors are simply   |
C  omitted.                                                                    |
C                                                                              |
C                                                                              |
C  0.2 Data type for I/O with the falcON code                                  |
C  ------------------------------------------                                  |
C  The code will take arrays of positions and accelerations. It then copies    |
C  their contents into its internal data types. The numerical representation   |
C  of the latter can be chosen independently from that of the former. Both are |
C  controlled by the macro "PRECISION" in file make.falcON:                    |
C                                                                              |
C                                       internal           arrays              |
C                                     -------------------------------          |
C  PRECISION	:= SINGLE_SINGLE        32 bit             32 bit              |
C  PRECISION	:= SINGLE_DOUBLE        32 bit             64 bit              |
C  PRECISION	:= DOUBLE_SINGLE        64 bit             32 bit              |
C  PRECISION	:= DOUBLE_DOUBLE        64 bit             64 bit              |
C                                                                              |
C  Below, we will use "INPUT_TYPE" for either "REAL" of "DOUBLE PRECISION",    |
C  depending on the value of "PRECISION". One way to ensure that any           |
C  application program has the correct input argument type is to add           |
C                                                                              |
C  #if defined(SINGLE_DOUBLE) || defined(DOUBLE_DOUBLE)                        |
C  #  define INPUT_TYPE DOUBLE PRECISION                                       |
C  #else                                                                       |
C  #  define INPUT_TYPE REAL                                                   |
C  #endif                                                                      |
C                                                                              |
C  at the top of the FORTRAN source code and run it through the C-preprocessor |
C  before compilation, see files src/FORTRAN/TestGravF.F and make.falcON for   |
C  an example.                                                                 |
C                                                                              |
C                                                                              |
C  0.3 FORTRAN specific problems                                               |
C  -----------------------------                                               |
C  Note that with FORTRAN you won't recognize a wrong argument type at all:    |
C  the code will compile, run, and give wrong result, which are not            |
C  necessarily recognisable as such. As a consequence,                         |
C                                                                              |
C      +----------------------------------------------------------------+      |
C      |  YOU MUST BE EXTERMELY CAREFUL TO GET THE ARGUMENT TYPES RIGHT |      |
C      +----------------------------------------------------------------+      |
C                                                                              |
C  This is a generic FORTRAN problem and not a problem of the code. The        |
C  best and only secure way to avoid it, is to abolish the use of FORTRAN.     |
C                                                                              |
C  You have been warned!                                                       |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  1 INITIALISATION AND CLEARING UP                                            |
C  ================================                                            |
C                                                                              |
C  1.1 BEFORE USING THE CODE                                                   |
C  -------------------------                                                   |
C  In order to use the code, one has first to initialize it. This is done      |
C  using the following routine, which actually does no computation at all.     |
C                                                                              |
      EXTERNAL FALCON_INITIALIZE         !  initialize falcON code             |
C                                                                              |
C  The syntax is a follows:                                                    |
C                                                                              |
C     INTEGER    N                                                             |
C     PARAMETER (N=10000)          number of bodies                            |
C     INTEGER    FL(N)             array with flags, see below                 |
C     INPUT_TYPE M(N)              array with masses                           |
C     INPUT_TYPE X(N),Y(N),Z(N)    arrays with position components             |
C     INPUT_TYPE E(N)              array with indiv softening lengths          |
C     INPUT_TYPE AX(N),AY(N),AZ(N) arrays for acceleration components          |
C     INPUT_TYPE P(N)              array for potentials                        |
C     INPUT_TYPE RH(N)             array for mass densities                    |
C     INPUT_TYPE EPS               global softening length, see below          |
C     INPUT_TYPE TH                opening angle, see below                    |
C     INTEGER    K                 type of softening kernel, see below         |
C                                                                              |
C     CALL FALCON_INITIALIZE(FL,M,X,Y,Z,AX,AY,AZ,P,RH, N,EPS,TH,K)             |
C                                                                              |
C  or, for the proprietary code, if the C-macro "falcON_INDI" is defined, see  |
C  line 18 of the "makefile":                                                  |
C                                                                              |
C     CALL FALCON_INITIALIZE(FL,M,X,Y,Z,E,AX,AY,AZ,P,RH, N,EPS,TH,K)           |
C                                                                              |
C  The first 11 arguments specify the sink and source properties of the        |
C  bodies. Each body has the sink properties: position (x,y,z), mass,          |
C  softening length (for the case of individual softening lengths), and flag   |
C  (see section 2 below). The source properties are: acceleration(ax,ay,az),   |
C  potential and mass-density. If fixed softening is used, a NULL pointer may  |
C  be used instead of an array holding eps_i, i.e.                             |
C                                                                              |
C     CALL FALCON_INITIALIZE(FL,M,X,Y,Z,%VAL(0),AX,AY,AZ,P,RH, N,EPS,TH,K,D)   |
C                                                                              |
C  The same applies to the arrays for density and potential: if a NULL         |
C  pointer is given, they will never be updated (but possibly computed).       |
C                                                                              |
C  The last 4 arguments are:                                                   |
C                                                                              |
C  N:                  size of arrays == number of bodies.                     |
C  EPS:   if(EPS>=0):  globally fixed softening length                         |
C         if(EPS< 0):  use individual softening lengths provided in array      |
C  TH:    if(TH>0):    theta = theta(M) with theta_min = TH (Dehnen 2002)      |
C         if(TH<0):    theta = |TH| = const                                    |
C         RECOMMENDED: TH    = 0.55                                            |
C  K:     0,1,2,3:     Ferrers sphere of this index (Dehnen 2001)              |
C         10,11,12,13  compensating kernel (Dehnen 2001) of index (KERN-10)    |
C         20,21,22,23  Plummer (20) and related (Dehnen & Teuben, 2002)        |
C         RECOMMENDED: K = 21 or 22                                            |
C                                                                              |
C  If FALCON_INIT() is called more than once, only the last initialisation     |
C  applies, the older ones will be deleted. That is, you can have only one     |
C  tree at a time -- use the C++ version if you need more.                     |
C                                                                              |
C  With the functions                                                          |
C                                                                              |
      EXTERNAL FALCON_RESETSOFTENING     !  reset softening parameters         |
C     syntax:                                                                  |
C     INTEGER    K                         softening kernel (as above)         |
C     INPUT_TYPE EPS                       global softening length             |
C     CALL FALCON_RESETSOFTENING(EPS,K)    set softening properties            |
C                                                                              |
C  and                                                                         |
C                                                                              |
      EXTERNAL FALCON_RESETOPENING       !  reset opening criterion            |
C     syntax:                                                                  |
C     INPUT_TYPE TH                        opening angle                       |
C     CALL FALCON_RESETOPENING(TH)         set opening properties              |
C                                                                              |
C  it is possible to change, after initialisation, the softening length and    |
C  kernel as well as the opening criterion.                                    |
C                                                                              |
C                                                                              |
C                                                                              |
C  1.2 AFTER USING THE CODE                                                    |
C  ------------------------                                                    |
C  After using the code, one should delete allocated memory by calling         |
C                                                                              |
      EXTERNAL FALCON_CLEARUP            !  de-allocated memory                |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  2 MEANING OF THE BODY FLAGS                                                 |
C  ===========================                                                 |
C                                                                              |
C  The body flags are integers with the following meaning:                     |
C                                                                              |
C   bit  value     meaning                                                     |
C   --------------------------------------------------------------             |
C     1      1     this body is active, i.e. wants update                      |
C     2      2     don't load this body into the tree, ignore it               |
C     3      4     this body is a SPH particle                                 |
C     4      8     this body is a sticky particle                              |
C   i>4   2^(i-1)  not used                                                    |
C                                                                              |
C                                                                              |
C  The default, flag=0, represents a plain body that is inactive, but still    |
C  source of gravity.                                                          |
C  The flag is obtained by setting the bits or, equivalently, adding the       |
C  values. The flags are not used by FALCON_INIT() and will unfold their       |
C  effect only when The functions in section 3 and 4 below are called.         |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  3 FORCE APPROXIMATION AND RELATED                                           |
C  =================================                                           |
C                                                                              |
C  3.1 Generating a tree structure                                             |
C  -------------------------------                                             |
C  In order to establish a hierarchical tree structure from all bodies that    |
C  are not flagged to be ignored, and subsequently pre-compute the centers     |
C  of mass and sizes of the cells, use                                         |
C                                                                              |
      EXTERNAL FALCON_GROW               !  grow new tree structure            |
C     syntax:                                                                  |
C     INTEGER NCRIT                                                            |
C     CALL FALCON_GROW(NCRIT)                                                  |
C                                                                              |
C  which creates a tree from scratch. Cells containing NCRIT or less bodies    |
C  will NOT be splitted. Experiments showed that for a full force calculation  |
C  NCRIT = 6-8 is about optimal and results in a few % decrease in CPU time    |
C  consumption (compared to NCRIT=1) and about 20% reduction in memory.        |
C                                                                              |
C                                                                              |
C  You may instead also re-use an old tree structure and only re-compute the   |
C  centers of mass and sizes of the cells (since the bodies have moved), by    |
C  using                                                                       |
C                                                                              |
      EXTERNAL FALCON_REUSE              !  re-use old tree structure          |
C     syntax:                                                                  |
C     CALL FALCON_REUSE()                                                      |
C                                                                              |
C  In this case, the logical linkage between cells and bodies is preserved     |
C  and no new tree structure is established. However, the code accounts for    |
C  the fact that the bodies might have moved. If the bodies have moved a lot,  |
C  the sizes of the cells will be much larger than the physical size of the    |
C  associated boxes, and the tree traversal will be very inefficient.          |
C  However, if the bodies have moved only little, the force compuation is      |
C  hardly slowed down and re-using an old tree saves the costs for             |
C  establishing a new one. Note that FALCON_REUSE() does not allow for any     |
C  change in the flags indicating whether or not a body shall be ignored:      |
C  changing this flag will have an effect only at the next call of             |
C  FALCON_GROW() or FALCON_GROW_CENTERED().                                    |
C                                                                              |
C                      +---------------------------------+                     |
C                      | USE THIS OPTION VERY CAREFULLY! |                     |
C                      +---------------------------------+                     |
C                                                                              |
C  You have been warned!                                                       |
C                                                                              |
C                                                                              |
C  3.2 Approximating the forces                                                |
C  ----------------------------                                                |
C  Once a tree structure is established, you can compute the forces of all     |
C  bodies flagged being active and due to all bodies not flagged to be ignored |
C  by                                                                          |
C                                                                              |
      EXTERNAL FALCON_APPROX_GRAV        !  approximate, e.g., forces          |
C     syntax:                                                                  |
C     FALCON_APPROX_GRAV()                                                     |
C                                                                              |
C  after a call to FALCON_GROW, FALCON_GROW_CENTERED, or FALCON_REUSE.         |
C  See file src/FORTRAN/TestGravF.f for an example application.                |
C                                                                              |
C                                                                              |
C  3.3 Estimating mass- and number-density                                     |
C  ---------------------------------------                                     |
C  There is also the possibility to obtain a rough estimate of the mass- or    |
C  number density of bodies in the neighbourhood of every body flagged being   |
C  active, via                                                                 |
C                                                                              |
      EXTERNAL FALCON_ESTIMATE_RHO       !  estimate the mass density          |
C     syntax:                                                                  |
C     INTEGER NX                                                               |
C     FALCON_ESTIMATE_RHO(NX)                                                  |
C                                                                              |
C  and                                                                         |
C                                                                              |
      EXTERNAL FALCON_ESTIMATE_N         !  estimate the number density        |
C     syntax:                                                                  |
C     INTEGER NX                                                               |
C     FALCON_ESTIMATE_N(NX)                                                    |
C                                                                              |
C  These estimates are simply the mean density within the smallest cells with  |
C  more then NX (1st arg) bodies. Note that when using test bodies (bodies     |
C  with zero or tiny mass), this guess for the mass density can have terrible  |
C  errors. Moreover, when the tree has not been grown but simply reused,       |
C  these estimate will not change. Be careful using these functions.           |
C  YOU HAVE BEEN WARNED.                                                       |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  4 STICKY-PARTICLE AND SPH SUPPORT                                           |
C  =================================                                           |
C                                                                              |
C  Once a tree structure is established (i.e. after a call to FALCON_GROW,     |
C  TREE_GROW_CENTERED, or FALCON_REUSE), you can also use it to create         |
C  interaction lists for SPH and sticky particles.                             |
C                                                                              |
C                                                                              |
C  4.1. Sticky-particle support                                                |
C  ----------------------------                                                |
C  This is provided by the routine                                             |
C                                                                              |
      EXTERNAL FALCON_STICKY              !  sticky particle support           |
C     syntax:                                                                  |
C     INTEGER    N                                                             |
C     PARAMETER (N=10000)             number of bodies                         |
C     INPUT_TYPE SIZE(N)              array with body sizes                    |
C     INPUT_TYPE VX(N),VY(N),VZ(N)    array for velocity components            |
C     INTEGER    NI                                                            |
C     PARAMETER (NI=1000)             physical size of interaction list        |
C     INTEGER    I1(NI)               array for indices: 1st of pair           |
C     INTEGER    I2(NI)               array for indices: 2nd of pair           |
C     INTEGER    NA                   actual size of interaction list          |
C     INPUT_TYPE TAU                  time step tau                            |
C     CALL FALCON_STICKY(TAU,I1,I2,NI,NA,VX,VY,VZ,SIZE)                        |
C                                                                              |
C  which fills the interaction list (2nd/3rd arg) with all pairs {i,j} of      |
C  indices which satisfy the following three conditions.                       |
C                                                                              |
C     (1) both flags indicate sticky particles (see section 2),                |
C     (2) at least one is flagged being active,                                |
C     (3) | (x_i+t*v_i)-(x_j+t*v_j) | < size_i+size_j  for some t in [0,tau]   |
C                                                                              |
C  In case of overflow, i.e. if the number of such pairs exceeds the size      |
C  (4th arg) of the list (2nd/3rd arg), the routine issues a warning to        |
C  stderr, but does not truncate the search, rather it no longer copies pairs  |
C  into the interaction list. In this case, the value returned for the actual  |
C  size (5th arg) exceeds the maximum size (4th arg), but reflects the actual  |
C  number of interactions found.                                               |
C                                                                              |
C  The arrays for velocities and sizes are assumed to be of the same size as   |
C  those used for positions (given to FALCON_INITIALIZE). However, the code    |
C  accesses velocities only for those bodies which are flagged as sticky, ie.  |
C  if only the first, say, 1000 of 10000 bodies are sticky, it is sufficient   |
C  if the size of the arrays with velocities and body sizes are 1000.          |
C                                                                              |
C  See file src/FORTRAN/TestPairF.F for an example application.                |
C                                                                              |
C                                                                              |
C  4.2 SPH particle support                                                    |
C  ------------------------                                                    |
C  This is provided by the routine                                             |
C                                                                              |
      EXTERNAL FALCON_SPH                !  SPH particle support               |
C     syntax:                                                                  |
C     INTEGER    N                                                             |
C     PARAMETER (N=10000)             number of bodies                         |
C     INPUT_TYPE SIZE(N)              array with body sizes                    |
C     INTEGER    NI                                                            |
C     PARAMETER (NI=1000)             physical size or interaction list        |
C     INTEGER    I1(NI)               array for indices: 1st of pair           |
C     INTEGER    I2(NI)               array for indices: 2nd of pair           |
C     INTEGER    NA                   actual size or interaction list          |
C     CALL FALCON_SPH(I1,I2,NI,NA,SIZE)                                        |
C                                                                              |
C  which fills the interaction list (1st/2nd arg) with all pairs {i,j} of      |
C  indices which satisfy the following three conditions.                       |
C                                                                              |
C     (1) both flags indicate SPH particles (see section 2),                   |
C     (2) at least one is flagged being active,                                |
C     (3) | x_i - x_j | < max(size_i,size_j)                                   |
C                                                                              |
C  In case of overflow, i.e. if the number of such pairs found exceeds the     |
C  size (3rd arg) of the list (1st/2nd arg), the routine issues a warning to   |
C  stderr, but does not truncate the search, rather it ceases copying pairs    |
C  into the interaction list. In this case, the value returned for the actual  |
C  size (4th arg) exceeds the maximum size (3rd arg), but reflects the actual  |
C  number of interactions found.                                               |
C                                                                              |
C  See file src/FORTRAN/TestPairF.F for an example application.                |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  5 OTHER FEATURES                                                            |
C  ================                                                            |
C                                                                              |
C  There are many further routines by which the user can obtain information    |
C  about the behaviour of the code.                                            |
C                                                                              |
      REAL     FALCON_ROOT_CENTER        !  gives Ith comp. of root center     |
      EXTERNAL FALCON_ROOT_CENTER        !                                     |
C     syntax:                                                                  |
C     INTEGER I                                                                |
C     INPUT_TYPE  C = FALCON_ROOT_CENTER(I)                                    |
C                                                                              |
      REAL     FALCON_ROOT_RADIUS        !  gives radius of root cell          |
      EXTERNAL FALCON_ROOT_RADIUS        !                                     |
C     syntax:                                                                  |
C     INPUT_TYPE    R = FALCON_ROOT_RADIUS()                                   |
C                                                                              |
      REAL     FALCON_CURRENT_EPS        !  gives global eps                   |
      EXTERNAL FALCON_CURRENT_EPS        !                                     |
C     syntax:                                                                  |
C     INPUT_TYPE    EP=FALCON_CURRENT_EPS()                                    |
C                                                                              |
      INTEGER  FALCON_CURRENT_KERNEL     !  gives kernel info                  |
      EXTERNAL FALCON_CURRENT_KERNEL     !                                     |
C     syntax:                                                                  |
C     INTEGER K = FALCON_CURRENT_KERNEL()                                      |
C                                                                              |
      INTEGER  FALCON_SOFTENING          !  gives 0/1 : fixed/individual       |
      EXTERNAL FALCON_SOFTENING          !                                     |
C     syntax:                                                                  |
C     INTEGER S = FALCON_SOFTENING()                                           |
C                                                                              |
      INTEGER  FALCON_NO_CELLS           !  gives number of cells              |
      EXTERNAL FALCON_NO_CELLS           !                                     |
C     syntax:                                                                  |
C     INTEGER NCELLS = FALCON_NO_CELLS()                                       |
C                                                                              |
      INTEGER  FALCON_DEPTH              !  gives tree depth                   |
      EXTERNAL FALCON_DEPTH              !                                     |
C     syntax:                                                                  |
C     INTEGER DEPTH = FALCON_DEPTH()                                           |
C                                                                              |
      EXTERNAL FALCON_STATS              !  writes some stats to stdout        |
C     syntax:                                                                  |
C     CALL FALCON_STATS()                                                      |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  6 KNOWN BUGS AND PROBLEMS                                                   |
C  =========================                                                   |
C                                                                              |
C  6.1 Test-bodies are not possible                                            |
C  --------------------------------                                            |
C  A body that is loaded into the tree but has zero mass, will not acquire     |
C  any acceleration. This is because the code computes first the force = mass  |
C  times acceleration (it is symmetric and hence better suited for             |
C  computation of mutual interactions) and then divides by the mass to obtain  |
C  the acceleration. The only possible work-around this problem is to set      |
C  the mass of potential test bodies much smaller than the masses of source    |
C  bodies. However, this must be done such that the gravity of test bodies     |
C  is everywhere neglible compared to that of the source bodies.               |
C  Note, however, that this work-around is wasteful: it computes the forces    |
C  generated by the test bodies. (You have been warned!)                       |
C                                                                              |
C                                                                              |
C  6.2 Bodies at identical positions                                           |
C  ---------------------------------                                           |
C  The code cannot cope with more than Ncrit bodies at the same position       |
C  (within the floating point accuracy). This situation will lead to an        |
C  infinitely deep tree, i.e. the maximum allowed tree depth will be exceeded  |
C  and the code aborts with an error message.                                  |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  REFERENCES                                                                  |
C  ==========                                                                  |
C                                                                              |
C  Dehnen W., 2000, ApJ, 536, L39                                              |
C  Dehnen W., 2001, MNRAS, 324, 273                                            |
C  Dehnen W., 2002, J. Comp. Phys., in press                                   |
C  Dehnen W. & Teuben P.J., 2002, in preparation                               |
C                                                                              |
C------------------------------------------------------------------------------+
C                                                                              |
C                                                                              |
C  A ROUTINES USED FOR TEST PURPOSES ONLY                                      |
C  ======================================                                      |
C                                                                              |
C  Routines not used by the tree, but by the test program TestGravF.f          |
C                                                                              |
C  YOU MAY COMMENT THEM OUT IF YOU NEVER USE THEM.                             |
C                                                                              |
      INTEGER           CLOCK
      DOUBLE PRECISION  DRAND48,CLOCKS_PER_SECOND
      EXTERNAL          CLOCK,DRAND48,SRAND48,CLOCKS_PER_SECOND
C------------------------------------------------------------------------------+


