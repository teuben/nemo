/* 
 * HDF4 interface for NEMO
 *
 * Since including <hdf.h> within NEMO programs can cause name
 * collisions, this file needs to be compiled without -I$NEMOINC,
 * or at least not include NEMO include files, just a few pure
 * ANSI headers
 *
 * Note this file compiles fine even if you don't have seem to have 
 * hdf installed!
 *
 *
 * 12-mar-97   hardcoded the relevant things from hdf.h
 * 24-nov-03   prototypes gcc3
 */

#include <stdio.h>
#include <string.h>

#if defined(HAVE_HDF)

#include <hdf.h>

#else

#define DFNT_FLOAT32	5
#define DFNT_FLOAT64	6
#define DFNT_FLOAT128	7

#define DFNT_INT8	20
#define DFNT_UINT8	21
#define DFNT_INT16	22
#define DFNT_UINT16	23
#define DFNT_INT32	24
#define DFNT_UINT32	25
#define DFNT_INT64	26
#define DFNT_UINT64	27
#define DFNT_INT128	28
#define DFNT_UINT128	29

#define DFNT_UCHAR8	3
#define DFNT_CHAR16	42
#define DFNT_UCHAR16	43

#endif

void hdf_type_info(int type, char *msg)
{
    switch (type) {
        case DFNT_FLOAT32:  strcpy(msg,"FLOAT32");  break;
        case DFNT_FLOAT64:  strcpy(msg,"FLOAT64");  break;
        case DFNT_FLOAT128: strcpy(msg,"FLOAT128"); break;
        case DFNT_INT8:     strcpy(msg,"INT8");     break;
        case DFNT_UINT8:    strcpy(msg,"UINT8");    break;
        case DFNT_INT16:    strcpy(msg,"INT16");    break;
        case DFNT_UINT16:   strcpy(msg,"UINT16");   break;
        case DFNT_INT32:    strcpy(msg,"INT32");    break;
        case DFNT_UINT32:   strcpy(msg,"UINT32");   break;
        case DFNT_INT64:    strcpy(msg,"INT64");    break;
        case DFNT_UINT64:   strcpy(msg,"UINT64");   break;
        case DFNT_INT128:   strcpy(msg,"INT128");   break;
        case DFNT_UINT128:  strcpy(msg,"UINT128");  break;
        case DFNT_UCHAR8:   strcpy(msg,"UCHAR8");   break;
        case DFNT_CHAR16:   strcpy(msg,"CHAR16");   break;
        case DFNT_UCHAR16:  strcpy(msg,"UCHAR16");  break;
        default:            strcpy(msg,"unknown");  break;
    }
}
