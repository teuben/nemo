/*
 * yapp_long.h
 *
 * This include file will be included by yapp.h if -DYAPP_LONG is defined
 * somewhere (normally options.h or -DYAPP_LONG in the $NEMOBIN/cc script
 * In the old ``ARGS'' days It was created using the following command:
 *
 * grep ARGS $NEMOINC/yapp.h | awk '{printf("#define %s \tyapp_%s\n",$2,$2)}'
 *
 */ 
#define plinit  	yapp_plinit
#define plltype         yapp_plltype
#define plline  	yapp_plline
#define plmove  	yapp_plmove
#define plpoint         yapp_plpoint
#define plcircle        yapp_plcircle
#define plcross         yapp_plcross
#define plbox   	yapp_plbox
#define pltext  	yapp_pltext
#define pljust  	yapp_pljust
#define plflush         yapp_plflush
#define plframe         yapp_plframe
#define plstop  	yapp_plstop
#define plswap  	yapp_plswap
#define plxscale        yapp_plxscale
#define plyscale        yapp_plyscale
#define plcolor         yapp_plcolor
#define plncolor        yapp_plncolor
#define plpalette       yapp_plpalette
#define pl_interp       yapp_pl_interp
#define pl_screendump   yapp_pl_screendump
#define pl_frame        yapp_pl_frame
#define pl_getpoly      yapp_pl_getpoly
#define pl_matrix       yapp_pl_matrix
