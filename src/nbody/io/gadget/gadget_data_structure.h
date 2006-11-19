// ============================================================================
// Copyright Jean-Charles LAMBERT - 2006
// e-mail:   Jean-Charles.Lambert@oamp.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           2, place Le Verrier
//           13248 Marseille Cedex 4, France
//           CNRS U.M.R 6110
// ============================================================================
// gadget_data_structure.h                                                            
// ============================================================================
#ifndef GADGET_DATA_STRUCTURE_H
#define GADGET_DATA_STRUCTURE_H

typedef struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} t_io_header_1;

typedef struct s_data_type_header {
  int size_bytes;
  int items;
} t_data_type_header;

enum DATA_TYPE {
  INT   =4,
  FLOAT =4,
  DOUBLE=8,
  CHAR  =1
};

const int NB_DATA_HEADER = 14;
const t_data_type_header dth[NB_DATA_HEADER] = {
  {INT      ,6}, // int      npart[6];
  {DOUBLE   ,6}, // double   mass[6];
  {DOUBLE   ,1}, // double   time;
  {DOUBLE   ,1}, // double   redshift;
  {INT      ,1}, // int      flag_sfr;
  {INT      ,1}, // int      flag_feedback;
  {INT      ,6}, // int      npartTotal[6];
  {INT      ,1}, // int      flag_cooling;
  {INT      ,1}, // int      num_files;
  {DOUBLE   ,1}, // double   BoxSize;
  {DOUBLE   ,1}, // double   Omega0;
  {DOUBLE   ,1}, // double   OmegaLambda;
  {DOUBLE   ,1}, // double   HubbleParam; 
  {CHAR     ,96} // char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
};

typedef struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  int    Id;     // added by JCL for fast qsort particles reordering
#if 0
  float  Rho, U, Temp, Ne;
#endif
} t_particle_data;
#endif
// -----------------------------------------------------------------------------
