#include "xdrfuncs.h"

int xdr_header(XDR *xdrs, struct dump *header)
{
	int pad=0;
  
	if (xdr_double(xdrs,&header->time) != TRUE) return 0;
	if (xdr_int(xdrs,&header->nbodies) != TRUE) return 0;
	if (xdr_int(xdrs,&header->ndim) != TRUE) return 0;
	if (xdr_int(xdrs,&header->nsph) != TRUE) return 0;
	if (xdr_int(xdrs,&header->ndark) != TRUE) return 0;
	if (xdr_int(xdrs,&header->nstar) != TRUE) return 0;
	if (xdr_int(xdrs,&pad) != TRUE) return 0;
	return 1;
	}


int xdr_gas(XDR *xdrs,struct gas_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->rho) != TRUE) return 0;
	if (xdr_float(xdrs,&p->temp) != TRUE) return 0;
	if (xdr_float(xdrs,&p->hsmooth) != TRUE) return 0;
	if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


int xdr_dark(XDR *xdrs,struct dark_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


int xdr_star(XDR *xdrs,struct star_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
	if (xdr_float(xdrs,&p->tform) != TRUE) return 0;
	if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  



