//-*- C -*- 
// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2014
//           Centre de donneeS Astrophysiques de Marseille (CeSAM)              
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Aix Marseille Universite, CNRS, LAM 
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS UMR 7326 
// ============================================================================
// Swig python interface for UNSIO
// ============================================================================

// python module name
%module py_unsio

%{
#define SWIG_FILE_WITH_INIT
#include "uns.h"
%}

%include "numpy.i"
%include "std_string.i"

%init %{
 import_array();
%}

// getArrayX return 1D numpy array
%apply ( int* DIM1, float ** ARGOUTVIEW_ARRAY1  )
      {( int* size, float ** farray             )};

%apply ( int* DIM1, int   ** ARGOUTVIEW_ARRAY1  )
      {( int* size, int   ** iarray             )};



// getValueX return float/int value
%apply float *OUTPUT { float * fvalue };
%apply int   *OUTPUT { int   * ivalue };

// rename methods because of overloading limitations with swig c++
%rename(getValueF) getData(const std::string,float *);
%rename(getArrayF) getData(const std::string,const std::string,int *,float **);
%rename(getArrayF) getData(const std::string,int *,float **);

%rename(getValueI) getData(const std::string,int *);
%rename(getArrayI) getData(const std::string,const std::string,int *,int **);
%rename(getArrayI) getData(const std::string,int *,int **);

// OUTPUT operations
//%apply float *INPUT { float * fvalue };
//%apply int   *INPUT { int   * ivalue };

%apply (int DIM1  , float * INPLACE_ARRAY1) {( int size, float * farray)};
%apply (int DIM1  , int   * INPLACE_ARRAY1) {( int size, int   * iarray)};

// rename methods because of overloading limitations with swig c++
%rename(setValueF) setData(const std::string,float);
%rename(setArrayF_do_not_used) setData(const std::string,const std::string,int ,float *,const bool _addr=false);
%rename(setArrayF_do_not_used) setData(const std::string,int ,float *,const bool _addr=false);

%rename(setValueI) setData(const std::string,int,const bool _addr=false);
%rename(setArrayI) setData(const std::string,const std::string,int ,int *,const bool _addr=false);
%rename(setArrayI) setData(const std::string,int ,int *);

// Parse the original header file
%include "uns.h"

%extend uns::CunsIn {
  // getFileName
  std::string getFileName() {
    std::string ret="";
    if ($self->isValid() && $self->snapshot!=NULL ) {
      ret=$self->snapshot->getFileName();
    } 
    return ret;
  }
  // getFileStructure
  std::string getFileStructure() {
    std::string ret="";
    if ($self->isValid() && $self->snapshot!=NULL ) {
      ret=$self->snapshot->getFileStructure();
    } 
    return ret;
  }
  // getInterfaceType
  std::string getInterfaceType() {
    std::string ret="";
    if ($self->isValid() && $self->snapshot!=NULL ) {
      ret=$self->snapshot->getInterfaceType();
    } 
    return ret;
  }
 };

 // Extend class
%extend uns::CunsOut {
  // we rewrite setArrayF because numpy array size is different from nbody for 3D arrays
  int setArrayF(const std::string  comp,const std::string  prop,
		int  size,float * farray, const bool _addr=false) {
    if (prop=="pos" || prop=="vel" || prop=="acc") size /= 3;
    int status = $self->snapshot->setData(comp,prop,size,farray,_addr);
    return status;
  }
  // we rewrite setArrayF because numpy array size is different from nbody for 3D arrays
  int setArrayF(const std::string  prop,
		int  size,float * farray, const bool _addr=false) {
    if (prop=="pos" || prop=="vel" || prop=="acc") size /= 3;
    int status = $self->snapshot->setData(prop,size,farray,_addr);
    return status;
  }
 };

