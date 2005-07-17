/*
 *   FIE.C:  function parser and evalutor - uses a interpreter-code
 *
 *
 *   This code is taken from the Gipsy package, and as such has some
 *   copyleft's attached, which are shamefully not all reprinted here.
 *
 *   mods' for unix/NEMO:
 *      modified all short (INTEGER*2) to int (INTEGER*4)
 *
 *	lots of variables and functions are now static here
 *
 *      % can also be used as parameter reference,
 *	  in addition to the exising $ which was used in VMS
 *	  there is a clear advantage of using a % to a $ on
 *	  Unix, despite the fact that MSDOS uses them...
 *
 *	<stdinc.h> and changed all float's to real's
 *
 *      moved 'npar' up to a global (like codeptr, opcodeptr) because
 *      it defines a 'fie' (see savefie, loadfie)
 *
 *	renamed 'const' to 'constnt' (const is a reserved ANSI word)
 * 
 *      added an understanding of null's
 *
 *   updates:
 *	       15-dec-88 older version			  KGB
 *		7-jun-89 merged GIPSY and NEMO versions   PJT
 *		5-jul-89 working on f2c interface 	  PJT
 *	       15-mar-90 inline -> innline for GCC        PJT
 *             10-sep-90 allow saving/loading fie_code's  PJT
 *             18-feb-92 replaced all float's by real     PJT
 *	       25-feb-92 happy gcc2.0			  PJT
 *             18-may-92 set_xrandom for NEMO             PJT
 *             16-dec-96 moved xrandom decl. up	          PJT
 *              7-apr-01 gcc warnings                     pjt
 *             20-jun-01 gcc3                             pjt
 *             26-aug-01 added sind/cosd/tand             pjt
 *             27-nov-01 fixed cosd(), it was sind()	  pjt
 *              4-dec-02 use MAXLINE for linelength       pjt
 *             13-nov-03 make it understand NULL          pjt
 *
 */
#include <stdinc.h>   /* stdinc is NEMO's stdio =- uses real{float/double} */
#include <ctype.h>
#include <math.h>

extern double xrandom(double,double);

#define byte char
#define bool int
#define false 0
#define true  1

/*  definitions/declarations for the function code  */

#define maxfiecode 100
#define bid     8       /* the number of bytes in a double real */

#define err    -1       /* an error occurred, too bad! */
#define hlt     0
#define add     1
#define sub     2
#define mul     3
#define div     4
#define neg     5
#define pwr     6
#define ldp     7
#define ldc     8

/* function call (i) has opcode  fie + i */

#define fie     9

static char *mnem[] = { "HLT","ADD","SUB","MUL","DIV","NEG","PWR","LDP",
                        "LDC","FIE" };

static union { byte opcode[bid];
        double c;        } fiecode[maxfiecode];

static int codeptr = 0;
static int opcodeptr = 0;
static int npar = 0;
static int have_null = 0;

static void fie_gencode(int opc);
static void fie_genconst(double cst);
static void fie_nextch(void);
static void fie_nextsym(void);
static void fie_expression(void);
static void fie_term(void);
static void fie_factor(void);
static void fie_function(void);
static void fie_error(void);
static void fie_null(void);
static void fie_push(double r);
static double fie_pi(void);
static double fie_rad(double arg1);
static double fie_deg(double arg1);
static double fie_sinc(double arg1);
static double fie_max(double arg1, double arg2);
static double fie_min(double arg1, double arg2);
static double fie_erf(double arg1);
static double fie_erfc(double arg1);
static double fie_mod(double arg1, double arg2);
static double fie_int(double arg1);
static double fie_sign(double arg1);
static double fie_ran(void);
static double fie_ranu(double arg1, double arg2);
static double fie_rang(double arg1, double arg2);
static double fie_ranp(double arg1);
static double fie_pop(void);

static void fie_gencode(int opc)
{
   fiecode[codeptr].opcode[opcodeptr++] = opc;
   if (opcodeptr == bid) {
      codeptr++; opcodeptr = 0;
   }
}

static void fie_genconst(double cst)
{
	fie_gencode(ldc);
	if (opcodeptr != 0) codeptr++;
	fiecode[codeptr++].c = cst;
	opcodeptr = 0;
}


/*  functions we know  */

#define maxfuncts 48
#define maxfunlen 10
#define maxarg    4

static char *functs[] = { "SIN"  , "ASIN" , "SINH" , "COS"  , "ACOS" , "COSH" ,
                   "TAN"  , "ATAN" , "TANH" , "ATAN2", "RAD"  , "DEG"  ,
                   "PI"   , "EXP"  , "LN"   , "LOG"  , "SQRT" , "ABS"  ,
                   "SINC" , "C"    , "G"    , "M"    , "ERF"  , "ERFC" ,
                   "K"    , "H"    , "P"    , "S"    , "MAX"  , "MIN"  ,
                   "MOD"  , "INT"  , "NINT" , "SIGN" , "UNDEF", "IFGT" ,
                   "IFLT" , "IFGE" , "IFLE" , "IFEQ" , "IFNE" , "RANU" ,
		   "RANG" , "RANP" , "SIND" , "COSD" , "TAND" , "NULL"};

static int nargs[] = {    1   ,    1   ,    1   ,    1   ,    1   ,    1   ,
                      1   ,    1   ,    1   ,    2   ,    1   ,    1   ,
                      0   ,    1   ,    1   ,    1   ,    1   ,    1   ,
                      1   ,    0   ,    0   ,    0   ,    1   ,    1   ,
                      0   ,    0   ,    0   ,    0   ,    2   ,    2   ,
                      2   ,    1   ,    1   ,    1   ,    0   ,    4   ,
                      4   ,    4   ,    4   ,    4   ,    4   ,    2   ,
		      2   ,    1   ,    1   ,    1   ,    1   ,    0 };


/*  definitions/declarations for the scanner and parser  */

#define end	0
#define plus	1
#define minus	2
#define times	3
#define divide	4
#define param	5
#define constnt	6
#define funct	7
#define lpar	8
#define rpar	9
#define comma	10
#define power   11

#define MAXPAR  32

#define MAXLINE 1024

static bool   parused[MAXPAR];
static char   innline[MAXLINE];
static int    pos, errorpos;
static char   ch;
static int    sym;
static int    curpar, curfun;
static double curconst;
static int    oddran, errorlev;

static void fie_nextch()
{
	if (ch != 0) ch = innline[pos++];
}
static void fie_nextsym()
{
	int    tenp;
	double f, frac;
	int    i;
	char   fun[maxfunlen];

	while (ch == ' ') fie_nextch();
	if (isdigit(ch)||(ch == '.')){
		curconst = 0.0;
		while (isdigit(ch)) {
			curconst = 10.0 * curconst + (ch - '0');
			fie_nextch();
		}
		if (ch == '.') {
			fie_nextch();
			f = 1;
			frac = 0;
			while (isdigit(ch)){
				frac = 10 * frac + (ch - '0');
				f = 10 * f;
				fie_nextch();
			}
			curconst = curconst + ((real) frac)/((real) f);	
		}
		if ((ch == 'E')||(ch == 'D')||(ch =='e')||(ch == 'd')) {
			fie_nextch();
			tenp = 1;
			frac = 0;
			if (ch == '+') {
				fie_nextch();
			} else if (ch == '-') {
				tenp = -tenp;
				fie_nextch();
			}
			while (isdigit(ch)) {
				frac = 10 * frac + (ch - '0');
				fie_nextch();
			}
			frac = frac * tenp;
			curconst = curconst * pow(10.0,(double) frac);
		}
		sym = constnt;
	} else if (isalpha(ch)) {
		i = 0;
		while ((isalpha(ch)||isdigit(ch)) && (i < maxfunlen)){
			fun[i++] = toupper(ch);
			fie_nextch();
		}
		fun[i] = 0;
		curfun = 0;
		while ((curfun < maxfuncts) && strcmp(fun,functs[curfun])){
			curfun++;
		}
		if (curfun == maxfuncts) fie_error();
		sym = funct;
	} else switch(ch){
		case '\0' : sym = end; fie_nextch(); break;
		case '+'  : sym = plus; fie_nextch(); break;
		case '-'  : sym = minus; fie_nextch(); break;
		case '*'  : sym = times; fie_nextch();
			    if (ch == '*') {
			    	sym = power; fie_nextch();
			     }
			    break;
		case '/'  : sym = divide; fie_nextch(); break;
		case '('  : sym = lpar; fie_nextch(); break;
		case ')'  : sym = rpar; fie_nextch(); break;
		case ','  : sym = comma; fie_nextch(); break;
		case '%'  :				/* added for UNIX */
		case '$'  : fie_nextch();
			    if (!isdigit(ch)) fie_error();
			    curpar = 0;
			    while (isdigit(ch)) {
				curpar = curpar * 10 + ch - '0';
				fie_nextch();
			    }
			    if (curpar >= MAXPAR) fie_error();
			    if (curpar > npar) npar = curpar;
			    parused[curpar] = true;
			    sym = param;
			    break;
		default: fie_error(); fie_nextch(); break;
	  }
}

int inifie(char *expr)			/* PJT: now int instead of short */
{
        int i, n = strlen(expr);
	
	for ( i=0 ; i<MAXPAR ; i++)
	    parused[i] = false;
	if (n > MAXLINE) error("inifie: No room (%d) to copy %s",n,expr);
	strcpy(innline,expr);
	oddran = 0;
	errorlev = 0;
	pos = 0;
	errorpos = 0;
	codeptr = 0;
	npar = 0;
	opcodeptr = 0;
	ch = ' ';
	fie_nextsym();
	fie_expression();
	if (sym != end) fie_error();
	if (errorpos == 0){
		errorpos = npar;
		/* for (i=1; i<=npar; i++)
		        if (!parused[i]) errorpos = 0; */
		}
	fie_gencode(hlt);
	return(errorpos);
}


static void fie_expression()
{
	int s;

	fie_term();
	while ((sym == plus) || (sym == minus)){
		s = sym;
		fie_nextsym();
		fie_term();
		if (s == plus) fie_gencode(add); else fie_gencode(sub);
	}
}

static void fie_term()
{
	int s;

	fie_factor();
	while ((sym == times) || (sym == divide)){
		s = sym;
		fie_nextsym();
		fie_factor();
		if (s == times) fie_gencode(mul); else fie_gencode(div);
	}
}


static void fie_factor()
{
	switch (sym) {
	case lpar   : fie_nextsym(); fie_expression();
		      if (sym == rpar) fie_nextsym(); else fie_error();
		      break;
	case minus  : fie_nextsym(); fie_factor(); fie_gencode(neg); break;
	case param  : fie_gencode(ldp); fie_gencode(curpar); fie_nextsym(); break;
	case constnt: fie_genconst(curconst); fie_nextsym(); break;
	case funct  : fie_function(); break;
	default     : fie_error(); break;
	}
	if (sym == power){
		fie_nextsym();
		fie_factor();
		fie_gencode(pwr);
	}
}



static void fie_function()
{
	int f = curfun;
	int n = nargs[curfun];

	fie_nextsym();
	if (n > 0) {
		if (sym == lpar) fie_nextsym(); else fie_error();
		while (n > 0) {
			fie_expression();
			if (--n > 0) {
				if (sym == comma) fie_nextsym();
				else fie_error();
			}
		}
		if (sym == rpar) fie_nextsym(); else fie_error();
	}
	fie_gencode( fie + f );
}


static void fie_error()
{
	if (errorpos == 0) errorpos = -pos;
	sym = end;
}

void dmpfie(void)
{
	int c,o,opc,op;

	c = 0;
	o = 0;

	do {
		op = fiecode[c].opcode[o++];
		if (o == bid) { c++ ; o = 0; }
		opc = op>fie ? fie : op;
		printf("     %s",mnem[opc]);
		if (opc == ldp) {
			opc = fiecode[c].opcode[o++];
			if (o == bid) { c++ ; o = 0; }
			printf("   %d",opc++);
		} else if (opc == fie) {
			printf("   %s",functs[op-opc]);
		} else if (opc == ldc) {
			if (o != 0) c++;
			printf("   %f",fiecode[c++].c);
			o = 0;
		}
		printf("\n");
	} while (opc != hlt);
}
		


#define stackmax 20

static double stack[stackmax];

static int sp;


static void fie_null(void)
{
  have_null = 1;
  warning("fie_null: i've seen null");
}

static void fie_push(double r)
{
	stack[++sp] = r;
}

static double fie_pi()
{
	double val;
	val = (double) 3.141592653589793238462643;
	return(val);
}

static double fie_rad(double arg1)
{
	double val;
	val = arg1 * 0.017453292519943295769237;
	return(val);
}

static double fie_deg(double arg1)
{
	double val;
	val = arg1 * 57.295779513082320876798155;
	return(val);
}

static double fie_sinc(double arg1)
{
	double val;
	if (fabs(arg1) < 1.0e-30) val = 1.0;
	else val = sin(arg1) / arg1;
	return(val);
}

static double fie_max(double arg1,double arg2)
{
	double val;
	if (arg1 > arg2) val = arg1; else val = arg2;
	return(val);
}

static double fie_min(double arg1,double arg2)
{
	double val;
	if (arg1 < arg2) val = arg1; else val = arg2;
	return(val);
}

static double fie_erf(double arg1)
{
   double p  =  0.327591100;
   double a1 =  0.254829592;
   double a2 = -0.284496736;
   double a3 =  1.421413741;
   double a4 = -1.453152027;
   double a5 =  1.061405429;
   double t1 = 1.0 / ( 1.0 + p * fabs(arg1));
   double t2 = t1*t1, t3 = t1*t2, t4 = t1*t3, t5 = t4*t1;
   if (arg1 > 0.0) {
      return(1.0 - (a1*t1+a2*t2+a3*t3+a4*t4+a5*t5)*exp(-arg1*arg1));
   } else {
      return((a1*t1+a2*t2+a3*t3+a4*t4+a5*t5)*exp(-arg1*arg1) - 1.0);
   }
}

static double fie_erfc(double arg1)
{
   return(1.0 - fie_erf(arg1));
}

static double fie_mod(double arg1,double arg2)
{
	double val;
	int   xxx = (int)(arg1/arg2);
	val = arg1 - xxx * arg2;
	return(val);
}

static double fie_int(double arg1)
{
	double val;
	int xxx = (int)arg1;
	val = xxx;
	return(val);
}

static double fie_sign(double arg1)
{
	if (arg1 == 0.0) return(0.0);
	else if (arg1 > 0.0) return(1.0);
	     else return(-1.0);
}

static double fie_ran(){
#if defined(NEMO_INC)
	return xrandom(0.0,1.0);
#else
	int xxx = rand();
	double val;
	val = (xxx+1) / 2147483648.0;
	return(val);
#endif
}

static double fie_ranu(double arg1,double arg2)
{
	double val;
        val = arg1 + fie_ran() * (arg2 - arg1);
        return(val);
}

static double fie_rang(double arg1,double arg2)
{
	double val, r1, r2;
	r1 = fie_ran();
	r2 = fie_ran();
	if (oddran == 0) {
		val = sqrt(-2*log(r1))*cos(6.283185307179586476925286*r2);
		oddran = 1;
	} else {
		val = sqrt(-2*log(r1))*cos(6.283185307179586476925286*r2);
		oddran = 0;
	}
	val = arg1 + fabs(arg2) * val;
	return(val);
}

static double fie_ranp(double arg1)
{
	double val, cum, p, f;
	if (arg1 < 40) {
		int xxx = (int)(fie_rang(arg1,sqrt(arg1))+0.5);
		val = xxx;
	} else {
		cum = exp(-arg1);
		p = cum;
		val = 0.0;
		f = fie_ran();
		while ( f >= cum) {
			val = val + 1.0;
			p = p * arg1 / val;
			cum = cum + p;
		}
	}
	return(val);
}
	
static double fie_pop()
{
	return( stack[sp--] );
}

void dofie(real *data, int *nop, real *results, real *errorval)
{
	int iop, i, c, o, opc;
	double undef = *errorval;
	double s,r;
	double params[MAXPAR], arg[maxarg];

	for ( iop = 1 ; iop <= *nop ; iop++ ){
		
	c = 0;
	o = 0;
	sp = 0;
	for (i=1; i<=npar; i++) {
		params[i] = data[iop + *nop * (i-1) - 1];
	}
	
	do {
		opc = fiecode[c].opcode[o++];
		if (o == bid) { c++ ; o = 0; }
		if (opc >= fie){
			int narg = nargs[opc-fie], n;
			for ( n=1 ; n<=narg ; n++ ) arg[narg-n] = fie_pop();
		}
		switch (opc){
		case hlt: break;
		case add: fie_push(fie_pop()+fie_pop()); break;
		case sub: r = fie_pop();
			  fie_push(fie_pop()-r);
			  break;
		case mul: fie_push(fie_pop()*fie_pop()); break;
		case div: r = fie_pop();
			  if (r == 0.0) opc = err;
			  else fie_push(fie_pop()/r);
			  break;
		case neg: fie_push(-fie_pop()); break;
		case pwr: r = fie_pop();
		          s = fie_pop();
		          if (s >= 0) fie_push(pow(s,r));
		          else {
				int t;
		          	int p = (int) r;
		          	double epsilon = 0.000001;
		          	if (fabs(r - p) <= epsilon){
		          		t = (p % 2 == 0) ? 1 : -1;
		          		fie_push( t * pow(fabs(s),r) );
		          	} else opc = err;
		          }
			  break;
		case ldp: opc = fiecode[c].opcode[o++];
			  if (o == bid) { c++ ; o = 0; }
			  fie_push( params[opc] );
			  break;
		case ldc: if (o != 0) c++;	
			  fie_push(fiecode[c++].c);
			  o = 0;
			  break;
		default:  switch(opc-fie){
			  case  0: fie_push(sin(arg[0])); break;
			  case  1: if (fabs(arg[0]) > 1) opc = err;
			  	   else fie_push(asin(arg[0]));    
			  	   break;
			  case  2: if (fabs(arg[0]) > 70) opc = err;
			  	   else fie_push(sinh(arg[0])); break;
			  case  3: fie_push(cos(arg[0])); break;
			  case  4: if (fabs(arg[0]) > 1) opc = err;
			  	   else fie_push(acos(arg[0]));
			  	   break;
			  case  5: if (fabs(arg[0]) > 70) opc = err;
			  	   else fie_push(cosh(arg[0])); break;
			  case  6: fie_push(tan(arg[0])); break;
			  case  7: fie_push(atan(arg[0])); break;
			  case  8: if (fabs(arg[0]) > 70) opc = err;
			  	   else fie_push(tanh(arg[0])); break;
			  case  9: fie_push(atan2(arg[0],arg[1])); break;
			  case 10: fie_push(fie_rad(arg[0])); break;
			  case 11: fie_push(fie_deg(arg[0])); break;
			  case 12: fie_push(fie_pi()); break;
			  case 13: if (fabs(arg[0]) > 70) opc =err;
			  	   else fie_push(exp(arg[0])); break;
			  case 14: if (arg[0] > 0) fie_push(log(arg[0]));
			  	   else opc = err;
			  	   break;
			  case 15: if (arg[0] > 0) fie_push(log10(arg[0]));
			           else opc = err;
			           break;
			  case 16: if (arg[0] < 0) opc = err;
			  	   else fie_push(sqrt(arg[0]));
			  	   break;
			  case 17: fie_push(fabs(arg[0])); break;
			  case 18: fie_push(fie_sinc(arg[0])); break;
			  case 19: fie_push( 2.997925e+8 ); break;
			  case 20: fie_push( 6.6732e-11 ); break;
			  case 21: fie_push( 1.99e30 ); break;
			  case 22: fie_push(fie_erf(arg[0])); break;
			  case 23: fie_push(fie_erfc(arg[0])); break;
			  case 24: fie_push( 1.380622e-23 ); break;
			  case 25: fie_push( 6.6256196e-34 ); break;
			  case 26: fie_push( 3.086e16 ); break;
			  case 27: fie_push( 5.66961e-8 ); break;
			  case 28: fie_push(fie_max(arg[0],arg[1])); break;
			  case 29: fie_push(fie_min(arg[0],arg[1])); break;
			  case 30: if (arg[1] == 0.0) opc = err; else
			           fie_push(fie_mod(arg[0],arg[1])); break;
			  case 31: fie_push(fie_int(arg[0])); break;
			  case 32: fie_push(fie_int(arg[0]+0.5)); break;
			  case 33: fie_push(fie_sign(arg[0])); break;
			  case 34: fie_push(undef); break;
			  case 35: if (arg[0] > arg[1]) fie_push(arg[2]);
			           else fie_push(arg[3]); break;
			  case 36: if (arg[0] < arg[1]) fie_push(arg[2]);
			           else fie_push(arg[3]); break;
			  case 37: if (arg[0] >= arg[1]) fie_push(arg[2]);
			           else fie_push(arg[3]); break;
			  case 38: if (arg[0] <= arg[1]) fie_push(arg[2]);
			           else fie_push(arg[3]); break;
			  case 39: if (arg[0] == arg[1]) fie_push(arg[2]);
			           else fie_push(arg[3]); break;
			  case 40: if (arg[0] != arg[1]) fie_push(arg[2]);
			           else fie_push(arg[3]); break;
			  case 41: fie_push(fie_ranu(arg[0],arg[1])); break;
			  case 42: fie_push(fie_rang(arg[0],arg[1])); break;
			  case 43: if (arg[0] < 0) opc = err;
			  	   else fie_push(fie_ranp(arg[0])); 
			  	   break;
			  case 44: fie_push(sin(PI*arg[0]/180.0)); break;
			  case 45: fie_push(cos(PI*arg[0]/180.0)); break;
			  case 46: fie_push(tan(PI*arg[0]/180.0)); break;
		          case 47: fie_null(); break;
			  default: opc = err; break;
			  }
			  break;
		}
	} while ((opc != hlt) && (opc != err));
	if (opc == err)
		results[iop-1] = *errorval;
	else
		results[iop-1] = fie_pop();
	}
}

/* 
 * SAVEFIE, LOADFIE:  Quickly save and load fie's when multiple fie's
 *                    have to be 'online'
 *
 *  PJT - sep 90    written to replace default loadobj in NEMO
 */
static struct fie_slot {
    char            fiecode[bid*maxfiecode];
    int             npar;
    int             codeptr;
    int             opcodeptr;
    int             slot;
    struct fie_slot *fwd;
} *save_fie = NULL;
    
/* savefie:
 *      input: slot  -- requested slot (>0) to save in 
 *      return:         actual slot saved in (-1 means error)
 */

int savefie(int slot)
{
    struct fie_slot *psfie, *newf;
    int found=0, count=0;
    
    if (slot<=0) return(-1);    /* not supported yet */

    if (save_fie==NULL) {
        save_fie = (struct fie_slot *) allocate(sizeof(struct fie_slot));
        /* if (save_fie==NULL) return -1;  */   /* couln't allocate */
        save_fie->slot = 1;
        save_fie->fwd = NULL;
    }

    for(psfie = save_fie;;) {            /* loop to find the right one */
        if (psfie->slot == slot) {
            found = 1;
            break;
        }
        if (psfie->fwd == NULL)
            break;
        else
            psfie = psfie->fwd;
    }
    if (found==0 && psfie->fwd==NULL) {   /* allocate new one */
        newf = (struct fie_slot *) malloc(sizeof(struct fie_slot));
        if (newf==NULL) return(-1);     /* couln't allocate */
        newf->slot = psfie->slot + 1;
        newf->fwd = NULL;
        psfie->fwd = newf;
        psfie = newf;        /* work on this new one */
    }
    
    bcopy(fiecode,psfie->fiecode,bid*maxfiecode);
    psfie->npar = npar;
    psfie->codeptr = codeptr;
    psfie->opcodeptr = opcodeptr;
    if (slot==0)
        psfie->slot = count+1;

    dprintf(1,"SAVEFIE: slot %d, npar = %d codeptr=%d opcodeptr=%d ***\n",
         psfie->slot, npar, codeptr, opcodeptr);


    return(slot);
}

/* loadfie:
 *      input: slot  -- requested slot to reload from
 *      return:         1 = OK.
 *                      0 = didn't find it, use old default 
 */

int loadfie(int slot)
{
    struct fie_slot *psfie;

    if (save_fie==NULL) return(0);  /* didn't save, take default */

    psfie = save_fie;       /* start at root of linear list */
    while (psfie) {
        if (psfie->slot == slot) {
            bcopy(psfie->fiecode,fiecode,bid*maxfiecode);
            npar = psfie->npar;
            codeptr = psfie->codeptr;
            opcodeptr = psfie->opcodeptr;
            dprintf(1,"LOADFIE: slot %d, npar=%d codeptr=%d opcodeptr=%d ###\n",
                    slot, npar, codeptr, opcodeptr);
            return(1);      /* OK, found slot */
        } else
            psfie = psfie->fwd;
    }
    return(0);      /* no slot found !! */
}


#if defined(TESTBED)

main (argc, argv)
int argc;
char *argv[];
{
   real  x[256], y[256];
   char   expr[256];
   int    i, j;
   int    nout = 256, nop, ierd;
   char   type = 'f';
   real  errval = 0.0;

   for (i=1; i<argc; i++) {
      nop = inifie(argv[i]);
      if (nop<0) {
         printf("Error code %d with inifie(%s)\n",nop,argv[i]);
         exit(1);
      }
      printf ("*** argv[%d] : inifie returns nop=%d *** \n",i,nop);
      printf ("  DMPFIE displays operations to be done:\n");
      dmpfie();
      savefie(i);
   }


   printf ("\n Array x was initialized as 0..255\n");
   for (i=0; i<256; i++)
        x[i] = i;


   for (i=1; i<argc; i++) {
      j=loadfie(i);
      printf("loading slot %d loadfie -> %d\n",i,j);
      if (j<0) continue;
      
      j = 1;
      dofie (x,&j,y,&errval);
   
      printf ("y=%f (%f %f %f)\n",y[0],x[0],x[1],x[2]);
   }
}
#endif
