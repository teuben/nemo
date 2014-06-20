#ifndef XDRFUNCS_INCLUDED
#define XDRFUNCS_INCLUDED

#include <rpc/types.h>
#include <rpc/xdr.h>
#include "tipsydefs.std.h"

#ifdef __cplusplus
extern "C" {
#endif

int xdr_header(XDR *xdrs, struct dump *header);
int xdr_gas(XDR *xdrs,struct gas_particle *p);
int xdr_dark(XDR *xdrs,struct dark_particle *p);
int xdr_star(XDR *xdrs,struct star_particle *p);

#ifdef __cplusplus
}
#endif

#endif

