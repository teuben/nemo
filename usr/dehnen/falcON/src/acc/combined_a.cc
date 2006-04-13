//-----------------------------------------------------------------------------+
//                                                                             |
// combined.cc                                                                 |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.0    13/02/2004  WD created                                               |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
//==============================================================================
void dummy() { std::cerr<<" this is a dummy\n"; }

//=============================================================================#
// declare externally linkable C routines                                      |
//=============================================================================#
extern "C" {
  //----------------------------------------------------------------------------
  // typedefs for routines computing the potential                              
  //----------------------------------------------------------------------------
  typedef void (*pp_d)(const int*,const double*,double*,double*,const double*);
  typedef void (*pp_f)(const int*,const float *,float *,float *,const float *);
  //----------------------------------------------------------------------------
  // these are resolved in the NEMO library                                     
  //----------------------------------------------------------------------------
  void warning(char *, ...);
  void error(char *, ...);
  void nemo_dprintf(int, char *, ...);
  pp_d get_potential_double(char*, char*, char*);
  pp_f get_potential_float (char*, char*, char*);
  //----------------------------------------------------------------------------
  // these are to be defined below                                              
  //----------------------------------------------------------------------------
  void inipotential    (int*, double*, char *);
  void potential_float (int*, float *, float *, float *, float *);
  void potential_double(int*, double*, double*, double*, double*);
}
//=============================================================================#
// define C++ implementation of Galaxy potential                               |
//=============================================================================#
namespace {
  const int NPmax=10;
  int       NP;
  pp_f      Pf[NPmax];
  pp_d      Pd[NPmax];
  //----------------------------------------------------------------------------
  inline void SwallowRestofLine(std::istream& from)
  {
    char c;
    do from.get(c); while(from.good() && c !='\n');
  }
  // read line until character `#'; swallow rest of line                        
  inline void __read_line(std::istream& from, char*line, int const&n)
  {
    // 1. find first non-space character whereby:
    //    - return empty line if EOF, EOL, '#' are encountered
    do {
      from.get(*line);
      if(from.eof() || (*line)=='\n') {
	*line=0;
	return;
      }
      if((*line)=='#') {
	SwallowRestofLine(from);
	*line=0;
	return;
      }
    } while(isspace(*line));
    // 2. read line until first white-pace character
    //    swallow rest til EOL, if necessary 
    register char*l=line+1;
    const    char*L=line+n;
    while(l != L) {
      from.get(*l);
      if(from.eof() || (*l)=='\n') { *l=0; return; }
      if(isspace(*l)) {
	*l=0;
	SwallowRestofLine(from);
	return;
      }
      ++l;
    }
    error("combined: line longer than expected\n");
  }
  //----------------------------------------------------------------------------
  inline void read_line(std::istream& from, char*line, int const&n)
  {
    do {
      __read_line(from,line,n);
    } while(from.good() && *line == 0);
  }
}
//------------------------------------------------------------------------------
void inipotential(int*npar, double*par, char*file)
{
  const int size=200;
  std::ifstream in(file);
  if(!in.good())
    error("combined: couldn't open file \"%s\" for input\n",file);
  char Line[size], PotName[size], PotPars[size], PotFile[size];
  char*potname=0, *potpars=0, *potfile=0, unknown[9];
  read_line(in,Line,size);
  if(*Line == 0)
    error("combined: couldn't read data from file \"%s\"\n",file);
  if(strncmp(Line,"potname=",8))
    error("combined: first entry in file \"%s\" isn't \"potname=...\"\n",
	  file);
  NP=0;
  do {
    if(potname==0) {
      strcpy(PotName,Line+8);
      potname=PotName;
      potpars=0;
      potfile=0;
    }
    read_line(in,Line,size);
    if        (*Line == 0 || 0==strncmp(Line,"potname=",8)) {
      if(NP == NPmax)
	error("file \"%s\": more potnames than anticipated\n",file);
//       // TEST
//       std::cerr<<"loading potential "<<NP<<": \n"
// 	       <<"  potname="<<potname<<'\n';
//       if(potpars) std::cerr<<"  potpars="<<potpars<<'\n';
//       if(potfile) std::cerr<<"  potfile="<<potfile<<'\n';
//       // tensor_set
      Pd[NP] = get_potential_double(potname,potpars,potfile);
      Pf[NP] = get_potential_float (potname,potpars,potfile);
      NP++;
      potname=0;
    } else if(0==strncmp(Line,"potpars=",8)) {
      if(potpars)
	warning("additional \"potpars=\" in file \"%s\" ignored\n",file);
      else {
	strcpy(PotPars,Line+8);
	potpars = PotPars;
      }
    } else if(0==strncmp(Line,"potfile=",8)) {
      if(potfile)
	warning("additional \"potfile=\" in file \"%s\" ignored\n",file);
      else {
	strcpy(PotFile,Line+8);
	potfile = PotFile;
	if(0==strcmp(file,potfile))
	  error("recursion potfile=file not allowed\n");
      }
    } else {
      strncpy(unknown,Line,8); unknown[9]=0;
      warning("entry \"%s\" in file \"%s\" ignored\n",
	      unknown,file);
    }
  } while(*Line != 0);
}
//------------------------------------------------------------------------------
void potential_float(int*NDIM, float*X, float*F, float*P, float*T) {
  if(*NDIM == 2) {
    F[0] = F[1] = P[0] = 0.;
    float Fi[2], Pi;
    for(int IP=0; IP!=NP; ++IP) {
      (Pf[IP])(NDIM,X,Fi,&Pi,T);
      P[0] += Pi;
      F[0] += Fi[0];
      F[1] += Fi[1];
    }
  } else {
    F[0] = F[1] = F[2] = P[0] = 0.;
    float Fi[3], Pi;
    for(int IP=0; IP!=NP; ++IP) {
      (Pf[IP])(NDIM,X,Fi,&Pi,T);
      P[0] += Pi;
      F[0] += Fi[0];
      F[1] += Fi[1];
      F[2] += Fi[2];
    }
  }
}
//------------------------------------------------------------------------------
void potential_double(int*NDIM, double*X, double*F, double*P, double*T) {
  if(*NDIM == 2) {
    F[0] = F[1] = P[0] = 0.;
    double Fi[2], Pi;
    for(int IP=0; IP!=NP; ++IP) {
      (Pd[IP])(NDIM,X,Fi,&Pi,T);
      P[0] += Pi;
      F[0] += Fi[0];
      F[1] += Fi[1];
    }
  } else {
    F[0] = F[1] = F[2] = P[0] = 0.;
    double Fi[3], Pi;
    for(int IP=0; IP!=NP; ++IP) {
      (Pd[IP])(NDIM,X,Fi,&Pi,T);
      P[0] += Pi;
      F[0] += Fi[0];
      F[1] += Fi[1];
      F[2] += Fi[2];
    }
  }
}
//------------------------------------------------------------------------------
