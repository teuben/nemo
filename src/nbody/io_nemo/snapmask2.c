/* -------------------------------------------------------------- *\
|* snapmask2.c :	JC LAMBERT	le 14-Mar-95	V1.0
|*					le 15-Nov-95    V1.1
|*                      PJT                14-Jul-01    V2.0 (Bastille day!)
|*
|* Ce programme permet de copier les particules selectionnees a
|* partir d'un snapshot en entree, vers un snapshot en sortie.
|* Ce programme est 15 a 20 fois plus rapide que le programme
|* snapmask du package NEMO lorsque le nombre de particules est tres
|* grand ( > 100000 ).
|*
|* Materiel     : Sparc STATION
|* OS           : Solaris 2.x
|* Langage      : C
|* Version Nemo : 2.0 Juillet 1994 
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Fichiers d'inclusion
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#define REALLOC

#if 0
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#endif

/* -------------------------------------------------------------- *\
|* variables globales
\* -------------------------------------------------------------- */

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			  ascii input file name ",
    "out=???\n			  snapshot output file name ",
    "select=all\n		  select string ",
    "times=all\n		  times select string ",
    "step=1\n			  copy only step by step snapshots",
    "VERSION=2.0\n		  14-Jul-01 PJT",
    NULL,
};

string usage = "mask out particles while copying snapshot data (w/fast io_nemo)";


#define TIMEFUZZ	0.00001


/* -------------------------------------------------------------- *\ 
|* save_parameter :
|* Sauvegarde le nombre de particules et la position dans un fichier
|* NEMO
\* -------------------------------------------------------------- */
int save_parameter(outstr,time_d,nbody)
stream outstr;
double time_d;
int  nbody;
{
 real time = time_d;
 put_set(outstr, SnapShotTag);		/*   start snapshot output  */
 put_set(outstr, ParametersTag);
 put_data(outstr, TimeTag, RealType, &time, 0);
 put_data(outstr, NobjTag, IntType, &nbody, 0);
 put_tes(outstr, ParametersTag);  
}
/* -------------------------------------------------------------- *\ 
|* copy_onedim :
|* Copie les particules selectionnees d'un tableau a une dimension
|* Utilise pour copier la masse et le potentiel
\* -------------------------------------------------------------- */
int copy_onedim(Qsel,oneptr,nbody)
bool * Qsel[];
real * oneptr;
int nbody;
{
  real * op, * dp;
  int  i;

  for (i=0, op=oneptr, dp = oneptr; i<nbody; i++, op++)
    { if (!Qsel[0][i])
         continue;  /* pas de copies */
      if (op == dp)
        { dp++;
          continue; /* particules sont identiques */
        }
      *dp = *op; /* Destination = Origine */
      dp++;
    }
}

/* -------------------------------------------------------------- *\ 
|* copy_twodim :
|* Copie les particules selectionnees d'un tableau a deux dimensions
|* Utilise pour copier les coordonnees spatiales
\* -------------------------------------------------------------- */
int copy_twodim(Qsel,twoptr,nbody,step)
bool * Qsel[];
real * twoptr;
int nbody,step; /* phasespace : step = 6 //  acceleration : step = 3 */
{
  real * op, * dp;
  int  i,j;

  for (i=0, op=twoptr, dp = twoptr; i<nbody; i++, op+=step)
    { if (!Qsel[0][i])
         continue;  /* pas de copies */
      for (j=0;j<step;j++)
        if (op[j]!=dp[j])
            dp[j] = op[j]; /* Destination = Origine pour x,y,z,vx,vy,vz */
      dp+=step;
    }
}

/* -------------------------------------------------------------- *\ 
|*
\* -------------------------------------------------------------- */


/* -------------------------------------------------------------- *\ 
|* programme principal
\* -------------------------------------------------------------- */ 
nemo_main()         
{
    stream instr,      /* fichier en entree */
           outstr;     /* fichier en sortie */

    string timu,      /* temps en entree */          
           infile,
           headline,
           select_pts[1];

    int coordsys = CSCode(Cartesian, NDIM, 2);

    int    nret[1];

    double * timeptr  = NULL;
    real * phaseptr = NULL;
    real * massptr  = NULL;
    real * potptr   = NULL;
    real * accptr   = NULL;

    int    i, j,nbody, bits,nbody_out=0,step,n_step;

    bool   first = TRUE, out=FALSE,
         * Qsel[1]; /* pointeur pour selectionner */
    
    int * select_i[1];     /* tableau de points selectionner select */
    bool within();
    /*bool within(double, string, double);*/

    /* recuperation des parametres en entree */
    infile=getparam("in");
    instr = stropen(getparam("in"), "r");
    select_pts[0] = getparam("select");
   
    timu = getparam("times");
    step = getiparam("step");
    n_step = step;
    /* creation du fichier de dortie */
    outstr = stropen(getparam("out"),"w");
        
    Qsel[0] = NULL;


    get_history(instr);
    put_history(outstr);

    for (;;) 
      {    /* infinite loop, broken only when ran out of snapshots */
    	get_history(instr);                    /* read history */
        while (get_tag_ok(instr,HeadlineTag))
                headline = get_string(instr,HeadlineTag);
       
        if (!get_tag_ok(instr, SnapShotTag)) 
           break; /* check if done */

        get_set(instr, SnapShotTag);

           get_set(instr, ParametersTag);
              if (get_tag_ok(instr,TimeTag))
                { 
		  if (timeptr == NULL)
                      timeptr = (double *) allocate(sizeof(double));
                  if (timeptr == NULL)
		    error("%s: pas assez de memoire\n", getargv0());
                  get_data_coerced(instr, TimeTag, DoubleType, timeptr, 0);
		  dprintf(1,"Time read : %3.10f\n",*timeptr);
                }

              if (get_tag_ok(instr, NobjTag))
		{
		  get_data(instr, NobjTag, IntType, &nbody, 0);    
		  fprintf(stderr,"Nbody = [%d]\n",nbody);
		}
           get_tes(instr, ParametersTag);

           if (get_tag_ok(instr,ParticlesTag) &&
               (timeptr == NULL || streq(timu, "all") ||
                within((real) *timeptr, timu, (real) TIMEFUZZ)))

             if (n_step == step)
              { n_step= 1;
                 
                fprintf(stderr,"Copy time steps [%f]...\n",*timeptr);

                /* insertion du select_i pour les champs de saisie */
                for (j=0; j<1; j++)
                     if (Qsel[j]==NULL)
                       { 
                          Qsel[j]   = (bool *) allocate (nbody*sizeof(bool));
                          select_i[j] = (int *) allocate (nbody*sizeof(int));
                          if (Qsel[j] == NULL || select_i[j] == NULL)
                            error("%s: not enuf memory\n", getargv0());
                        
                          for (i=0 ; i < nbody ; i++)
                             Qsel[j][i] = FALSE;

                          if (!streq("all",select_pts[j]))
                            { for (i=0; i<nbody; i++)
                              { Qsel[j][i] = FALSE;
                                select_i[j][i] =-1;
                              }
                             nret[j] = nemoinpi(select_pts[j],
                                       select_i[j],nbody);
                             for (i=0; i < nret[j]; i++)
                                Qsel[j][select_i[j][i]]=TRUE;
                            }
                          else
                            { for (i=0; i<nbody; i++)
                               Qsel[j][i] = TRUE;   
                              nret[j] = nbody;
                            }

                          free(select_i[j]);
                      } 
                
                   /* Sauvegarde des Parametres */
                   nbody_out = nret[0];
                   out = TRUE;
                   save_parameter(outstr,*timeptr,nbody_out);
                   get_set(instr, ParticlesTag);
                   put_set(outstr, ParticlesTag);

                   /* sauvegarde du systeme de coordonnees */
                   put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
                   /* recuperation de la masse */
	           if (get_tag_ok(instr, MassTag))
                     { if (massptr == NULL)
                           massptr = (real *) malloc(sizeof(real) * nbody);
                        
		       if (massptr == NULL)
		         error("%s: pas assez de memoire\n", getargv0());
		       get_data_coerced(instr, MassTag, RealType, massptr,
				 nbody, 0);

                       /* copie les masses selectionnees */
                       copy_onedim(Qsel,massptr,nbody);

                       /* sauvegarde les masses */
                       put_data(outstr,MassTag,RealType,massptr,nbody_out,0);
                       
                       /* libere la memoire allouees aux masses */
                       free(massptr);
                       massptr = NULL;

	             }
                   /* recuperation des corrdonnes spatiales */
                   if (get_tag_ok(instr, PhaseSpaceTag))
                     { if (phaseptr == NULL)
                            phaseptr = (real *) allocate(sizeof(real)*2*NDIM*
                                                        nbody);
                         
                       if (phaseptr == NULL)
                          error("%s: pas assez de memoire\n", getargv0());
                       get_data_coerced(instr, PhaseSpaceTag,
                                     RealType,phaseptr,nbody, 2, NDIM, 0);
                       /* copie les coordonnees spatiales selectionnees */
                       copy_twodim(Qsel,phaseptr,nbody,6);
                       /* sauvegardes les coordonnees spatiales */
                       put_data(outstr, PhaseSpaceTag,
	                        RealType,phaseptr,nbody_out, 2, NDIM, 0);
                       /* libere la memoire allouees aux phasespace */
                       free(phaseptr);
                       phaseptr = NULL;
                     }

                   /* recuperation du potentiel */
	           if (get_tag_ok(instr, PotentialTag))
                     { if (potptr == NULL)
		            potptr = (real *) malloc(sizeof(real) * nbody);
                          
		       if (potptr == NULL)
		         error("%s: pas assez de memoire\n", getargv0());
		       get_data_coerced(instr, PotentialTag, RealType, potptr,
				 nbody, 0);

                       /* copie les masses selectionnees */
                       copy_onedim(Qsel,potptr,nbody);

                       /* sauvegarde les masses */
                       put_data(outstr,PotentialTag,
                                RealType,potptr,nbody_out,0);
                       
                       /* libere la memoire allouees aux masses */
                       free(potptr);
                       potptr = NULL;

	             }

		   /* recuperation de l'acceleration */
	           if (get_tag_ok(instr, AccelerationTag))
                     { if (accptr == NULL)
		            accptr = (real *) malloc(sizeof(real) * nbody * 3);
                          
		       if (accptr == NULL)
		         error("%s: pas assez de memoire\n", getargv0());
		       get_data_coerced(instr, AccelerationTag, RealType, accptr,
				 nbody, 3,0);

                       /* copie les masses selectionnees */
                       copy_twodim(Qsel,accptr,nbody,3);

                       /* sauvegarde les masses */
                       put_data(outstr,AccelerationTag,
                                RealType,accptr,nbody_out, 3, 0);
                       
                       /* libere la memoire allouees aux masses */
                       free(accptr);
                       accptr = NULL;

	             }
                     /* flag du premier snapshot */
                     first = FALSE;
                     
      
                get_tes(instr, ParticlesTag);
                put_tes(outstr, ParticlesTag);

              }
            else
               n_step++;
        get_tes(instr, SnapShotTag);   
        if (out)
           put_tes(outstr,SnapShotTag);  
        out = FALSE;  
      }
    strclose(instr);    /* fermeture du fichier d'entree  */
    fclose(outstr);     /* fermeture du fichier de sortie */
    if (first) 
      {
    	warning("No snapshots processed");
        
      } 
    else  
      { /* liberation de la memoire */
        
      }
}
 
/* -------------------------------------------------------------- *\ 
|* Fin de snapmask2.c 
\* -------------------------------------------------------------- */ 
