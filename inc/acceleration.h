/*
 ******************************************************************************
 *
 * acceleration.h
 *
 * definition of pointer-to-functions for computing body accelerations
 *
 *  Notes:
 *  1 the 4th argument to get_acceleration() is boolean and returns whether
 *    masses are required as input for acceleration(). If they are not 
 *    required a NULL pointer may be given.
 *  2 the 5th argument to get_acceleration() is boolean and returns whether
 *    velocities are required as input for acceleration(). If they are not 
 *    required a NULL pointer may be given.
 *    Velocities may be used to compute friction forces, such as the drag
 *    a gaseous disk is generating on stars crossing it.
 *  3 arrays are passed as pointer to void. They must be either all of type
 *    float or all of type double as indicated by the last argument being 'f'
 *    or 'd', respectively.
 *  4 arrays of vector quantities are in the order x0,y0,z0, x1,y1,z1, ...
 *  5 if the pointer to flags is NULL, all bodies are supposed to be active,
 *    otherwise only those for which (f[i] & 1) is true.
 *  6 the argument "indicator" of acceleration() indicates whether the
 *    accelerations and potential shall be assigned or added.
 *    If bit 0 is set, the potential    is added, otherwise assigned,
 *    If bit 1 is set, the acceleration is added, otherwise assigned.
 *    So, 0 means both are assigned.
 *
 ******************************************************************************
 * version 0.0  17/06/2004  WD
 * version 0.1  22/06/2004  WD   allow for up to 10 fallbacks, no last_...
 *
 */

#ifndef _acceleration_h
#define _acceleration_h

#ifdef __cplusplus
extern "C" {
#endif

/*
 * declaration of pointer to acceleration()
 */

typedef void(*acc_pter)         /* return: void                            */
     (int,                      /* input:  number of dimensions            */
      double,                   /* input:  simulation time                 */
      int,                      /* input:  number bodies = size of arrays  */
      const void*,              /* input:  masses:         m[i]            */
      const void*,              /* input:  positions       (x,y,z)[i]      */
      const void*,              /* input:  velocities      (u,v,w)[i]      */
      const int *,              /* input:  flags           f[i]            */
      void*,                    /* output: potentials      p[i]            */
      void*,                    /* output: accelerations   (ax,ay,az)[i]   */
      int,                      /* input:  indicator                       */
      char);                    /* input:  type: 'f' or 'd'                */

/*
 * declaration of routine for obtaining a pointer defined above
 */

extern
acc_pter get_acceleration(               /* return: acceleration()         */
			  const char*,   /* input:  acc_name               */
			  const char*,   /* input:  acc_pars               */
			  const char*,   /* input:  acc_file               */
			  bool      *,   /* output: need masses?           */
			  bool      *);  /* output: need velocities?       */
  
#ifdef __cplusplus
}
#endif
#endif // _acceleration_h
