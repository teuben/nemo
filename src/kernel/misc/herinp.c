/*
 *   name:      HERINP
 *
 *   author:    K.G. Begeman
 *
 *   routines:  herinp   (C-callable - expert mode) 
 *              dcdchar  (fortran-callable) - to be generated with f2cvv
 *              dcddble  (fortran-callable) - to be generated with f2cvv
 *              dcdreal  (fortran-callable) - to be generated with f2cvv
 *              dcdint   (fortran-callable) - to be generated with f2cvv
 *              dcdlog   (fortran-callable) - to be generated with f2cvv
 *
 *   f2c interface definitions:
 *
 * @ integer function dcdchar( character, character, integer, integer )
 * @ integer function dcdint ( character, integer, integer, integer )
 * @ integer function dcddble( character, double precision, integer, integer )
 * @ integer function dcdreal( character, real, integer, integer )
 * @ integer function dcdlog ( character, logical, integer, integer )
 *
 *   errors:     #    explanation
 *              -11   bad call
 *              -12   unknown function
 *              -13   syntax error
 *              -14   illegal character
 *              -15   wrong repeat argument
 *              -16   wrong number of arguments
 *              -17   arithmetic error
 *              -18   not enough internal memory
 *              -19   conversion error
 *              -20   unequal list length
 *              -21   empty list
 *              -22   nested lists
 *              -23   output buffer overflow
 *              -24   floating overflow/underflow in conversion
 *
 *   updates:   16/dec/87 : Type C implemented (KGB)
 *              23/mar/88 : Argument NCHR added (KGB)
 *              14/apr/88 : Error if after , no values (KGB)
 *              13/oct/88 : Implemented ANSI F77 calls (KGB)
 *               1/mar/89 : Bug in decoding logicals removed (KGB)
 *               7/jun/89 : merged GIPSY and NEMO version   (PJT)
 *		19-jun-89 : F2C interface added (PJT)
 *		25-jan-90 : sun4 bug about fie.arg.execution order fixed (PJT)
 *              15-mar-90 : made GCC happy about fie return codes (PJT)
 *		11-nov-90 : made TAB equivalent to SPACE (PJT)
 *                          alternatively, at each "if (ch==' ')" instance one
 *                          could use the macro 'ispace(ch)' from <ctype.h>
 *              10-aug-93: allow ``1:1'' also to generate a list (1)       PJT
 *               5-dec-93: -DBIGLOOP allows do-loops with no upper limit   PJT
 *		 7-jan-94: atand() bug resolved				   pjt
 *		           and atand2() also with dcd_ calls instead	   pjt
 *		16-feb-97: removed some nexted external decl's             pjt
 *               7-apr-01: gcc warning                                     pjt
 *		20-jun-01: gcc3, const->hconst				   pjt
 *               3-nov-01: use NEMO's random numbers, not rand()           pjt
 *              10-nov-03: allow and handle "null" more gracefully         pjt
 *              24-jan-04: STACKMAX check
 *              26-jul-2016:    double precision log(min/max) = 308        pjt
 */

#define BIGLOOP /* comment this our if you want MAXSHORT as largest count */

/* if we use NEMO.H here, string.h seems to cause a linux parse problem */
#if 1
# include <stdinc.h>
#else
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
extern void warning(char *,...), error(char *,...);
#endif

#include <ctype.h>
#include "gipsyc.h"


#if defined(toupper)
#undef  toupper
#define toupper(ch) (((ch >= 'a') && (ch <= 'z')) ?  (ch - 'a' + 'A') : (ch))
#endif

extern logical fblank_(int *);	
extern void    setfblank_(void *);

#define byte    char
#define bool    int
#define DEFAULT 1

#define NATIVE_RAND  0

static void dcd_inifblank(void);
static int  dcd_round(double arg);
static void dcd_error(int errnum);
static void dcd_gencode(int opc);
static void dcd_genconst(double cst);
static void dcd_inilist(void);
static void dcd_beginlist(void);
static void dcd_endlist(void);
static void dcd_putlist(void);
static void dcd_nextch(void);
static void dcd_nextsym(void);
static void dcd_nextwr(void);
static void dcd_movenum(void);
static void dcd_loop(void);
static void dcd_expression(void);
static void dcd_term(void);
static void dcd_factor(void);
static void dcd_function(void);
static void dcd_list(void);
static void dcd_dump(void);
static void dcd_push(double r);
static double dcd_pop(void);
static double dcd_add(double arg1, double arg2);
static double dcd_sub(double arg1, double arg2);
static double dcd_mul(double arg1, double arg2);
static double dcd_div(double arg1, double arg2);
static double dcd_neg(double arg1);
static double dcd_pwr(double arg1, double arg2);
static double dcd_sin(double arg1);
static double dcd_asin(double arg1);
static double dcd_sinh(double arg1);
static double dcd_cos(double arg1);
static double dcd_acos(double arg1);
static double dcd_cosh(double arg1);
static double dcd_tan(double arg1);
static double dcd_atan(double arg1);
static double dcd_tanh(double arg1);
static double dcd_atan2(double arg1, double arg2);
static double dcd_rad(double arg1);
static double dcd_deg(double arg1);
static double dcd_pi(void);
static double dcd_exp(double arg1);
static double dcd_ln(double arg1);
static double dcd_log(double arg1);
static double dcd_sqrt(double arg1);
static double dcd_abs(double arg1);
static double dcd_sinc(double arg1);
static double dcd_max(double arg1, double arg2);
static double dcd_min(double arg1, double arg2);
static double dcd_erf(double arg1);
static double dcd_erfc(double arg1);
static double dcd_mod(double arg1, double arg2);
static double dcd_int(double arg1);
static double dcd_nint(double arg1);
static double dcd_sign(double arg1);
static double dcd_ifgt(double arg1, double arg2, double arg3, double arg4);
static double dcd_iflt(double arg1, double arg2, double arg3, double arg4);
static double dcd_ifge(double arg1, double arg2, double arg3, double arg4);
static double dcd_ifle(double arg1, double arg2, double arg3, double arg4);
static double dcd_ifeq(double arg1, double arg2, double arg3, double arg4);
static double dcd_ifne(double arg1, double arg2, double arg3, double arg4);
static double dcd_ran(void);
static double dcd_ranu(double arg1, double arg2);
static double dcd_rang(double arg1, double arg2);
static double dcd_ranp(double arg1);
static void   dcd_null(void);
static void   dcd_evaluate(int q);

/*  definitions/declarations for the function code  */

#define maxfiecode 1024

#define bif     2       /* number of floats in a double real                 */
#define bii     2       /* number of integers in a double real               */
#define bis     4       /* number of shorts in a double real                 */
#define bid     8       /* number of bytes in a double real                  */

#define err    -1       /* an error occurred, too bad!                       */
#define hlt     0
#define add     1
#define sub     2
#define mul     3
#define div     4
#define neg     5
#define pwr     6
#define ldc     7
#define lst     8
#define fie     9       /* function call (i) has opcode  fie + i             */

static char *mnem[] = {
   "HLT","ADD","SUB","MUL","DIV","NEG","PWR","LDC","LST","FIE"
};

static union {
   byte opcode[bid];
   double c;
} fiecode[maxfiecode], lstfiecode[maxfiecode];

static int codeptr;
static int opcodeptr;
static int lstcodeptr;
static int lstopcodeptr;

/*  functions we know  */

#define maxfuncts  52
#define maxfunlen  10
#define maxarg      4
#define maxbools    8
#define maxboollen  6

static char *functs[] = {
   "SIN"   , "ASIN"  , "SINH"  , "COS"   , "ACOS"  , "COSH"  ,
   "TAN"   , "ATAN"  , "TANH"  , "ATAN2" , "RAD"   , "DEG"   ,
   "PI"    , "EXP"   , "LN"    , "LOG"   , "SQRT"  , "ABS"   ,
   "SINC"  , "C"     , "G"     , "M"     , "ERF"   , "ERFC"  ,
   "K"     , "H"     , "P"     , "S"     , "MAX"   , "MIN"   ,
   "MOD"   , "INT"   , "NINT"  , "SIGN"  , "BLANK" , "IFGT"  ,
   "IFLT"  , "IFGE"  , "IFLE"  , "IFEQ"  , "IFNE"  , "RANU"  ,
   "RANG"  , "RANP"  , "SIND"  , "ASIND" , "COSD"  , "ACOSD" ,
   "TAND"  , "ATAND" , "ATAND2", "NULL"
};

static int nargs[]    = {
   1   ,    1   ,    1   ,    1   ,    1   ,    1   ,
   1   ,    1   ,    1   ,    2   ,    1   ,    1   ,
   0   ,    1   ,    1   ,    1   ,    1   ,    1   ,
   1   ,    0   ,    0   ,    0   ,    1   ,    1   ,
   0   ,    0   ,    0   ,    0   ,    2   ,    2   ,
   2   ,    1   ,    1   ,    1   ,    0   ,    4   ,
   4   ,    4   ,    4   ,    4   ,    4   ,    2   ,
   2   ,    1   ,    1   ,    1   ,    1   ,    1   ,
   1   ,    1   ,    2   ,    0
};

static char *bools[] = {
   "TRUE", "FALSE", "YES", "NO", "JA", "NEE", "1", "0"
};

static bool boolv[] = {
   true, false, true, false, true, false, true, false
};

/*  definitions/declarations for the scanner and parser  */

#define end      0
#define plus     1
#define minus    2
#define times    3
#define divide   4
#define hconst   5
#define funct    6
#define lpar     7
#define rpar     8
#define comma    9
#define power   10
#define blank   11
#define colon   12
#define lbrac   13
#define rbrac   14

static int    listlen[2];
static int    lists;
static bool   list;
static int    pos, mpos, errorpos, errornum;
static char   ch;
static int    sym;
static int    curfun, npar;
static double curconst;
static int    oddran;
static char   *cptr, *dptr;
static char   ctyp;
static int    ilen, nmax;
static double dval;
static int    nval;
static int    have_null;

static union {
   byte    b[bid];
   short   s;
   int     i;
   float   f;
   double  d;
} unie;

/* define internal blank value */
static double DCDBLANK;

static union {
   byte   bb[sizeof(double)];
   double dd;
} dcd_blank;

/* some systems have these already defined....but we need them here*/
/* as long as these are not out of bounds with your local machine  */
#define MAXFLOAT        1.2e+37  /* cray */
#define MINFLOAT        0.8e-37  /* cray */
#define MAXSHORT        32767.5  /* cray */
#define MINSHORT       -32768.5
#define MAXINT     2147483647.5  /* cray */
#define MININT    -2147483648.5
/* these are usually not defined */

#if 0
/* single precision */
#define MAXLOG             38.0
#define MINLOG            -38.0
#else
#define MAXLOG             308.0
#define MINLOG            -308.0
#endif

static void dcd_inifblank()
{
   int i = 0;
#if defined(VMS)
   for (i = 0; i < sizeof(double); i++) dcd_blank.bb[i] = 0xff;
#else
   for (i = 0; i < sizeof(double); i++) dcd_blank.bb[i] = 0x77;
#endif
   DCDBLANK = dcd_blank.dd;
}

static int dcd_round(double arg)
{
   int val;
   if (arg > 0.0) {
      val = (int) (arg + 0.5);
   } else {
      val = (int) (arg - 0.5);
   }
   return(val);
}

static void dcd_error(int errnum)
{
   if (errornum == 0) {
      if (errorpos == 0) errorpos = pos;
      sym = err;
      errornum = errnum;
   }
}

static void dcd_gencode(int opc)
{
   if (errornum != 0) return;
   if (list) {
      lstfiecode[lstcodeptr].opcode[lstopcodeptr++] = opc;
      if (lstopcodeptr == bid) {
         lstcodeptr++;
         lstopcodeptr = 0;
      }
      if (lstcodeptr == maxfiecode) dcd_error(-18);  /* buffer overflow */
   } else {
      fiecode[codeptr].opcode[opcodeptr++] = opc;
      if (opcodeptr == bid) {
         codeptr++;
         opcodeptr = 0;
      }
      if (codeptr == maxfiecode) dcd_error(-18);     /* buffer overflow */
   }
}

static void dcd_genconst(double cst)
{
   dcd_gencode(ldc);
   if (errornum != 0) return;
   if (list) {
      if (lstopcodeptr != 0) lstcodeptr++;
      if (lstcodeptr == maxfiecode) {
         dcd_error(-18);
         return;
      }
      lstfiecode[lstcodeptr++].c = cst;
      lstopcodeptr = 0;
   } else {
      if (opcodeptr != 0) codeptr++;
      if (codeptr == maxfiecode) {
         dcd_error(-18);
         return;
      }
      fiecode[codeptr++].c = cst;
      opcodeptr = 0;
   }
}

static void dcd_inilist()
{
   lists = 0;
   list = false;
}

static void dcd_beginlist()
{
   int l;
   if (errornum != 0) return;
   dcd_gencode(lst);
   lstcodeptr = 0;
   lstopcodeptr = 0;
   list = true;
   if (lists++ > 0) {
      l = 1;
   } else {
      l = 0;
   }
   listlen[l] = 0;
}

static void dcd_endlist()
{
   int l;
   list = false;
   if (lists > 1) l = 1; else l = 0;
   if (listlen[l] == 0) dcd_error(-21);          /* empty list */
   if (listlen[0] != listlen[l]) dcd_error(-20); /* unequal list length */
}

static void dcd_putlist()
{
   int l;
   if (errornum != 0) return;
   if (opcodeptr != 0) {
      codeptr++;
      opcodeptr = 0; 
   }
   if (codeptr == maxfiecode) {
      dcd_error(-18);
      return;
   }
   fiecode[codeptr++].c = dval;
   if (codeptr == maxfiecode) {
      dcd_error(-18);
      return;
   }
   if (lists > 1) l = 1; else l = 0;
   listlen[l]++;
}

static void dcd_nextch()
{
   if (pos++ < mpos) ch = *cptr++; else ch = '\0';
}

static void dcd_nextsym()
{
   double f, frac;
   int    tenp;
   int    i;
   char   fun[maxfunlen];
   if (sym == err) return;
   if (isdigit(ch)||(ch == '.')) {
      curconst = 0.0;
      while (isdigit(ch)) {
         if (errornum == 0) curconst = 10.0 * curconst + ch - '0';
         if (curconst > MAXFLOAT) dcd_error(-24);        /* floating overflow */
         dcd_nextch();
      };
      if (ch == '.') {
         dcd_nextch();
         f = 1;
         frac = 0;
         while (isdigit(ch)) {
            if (errornum == 0) {
               frac = 10 * frac + (ch - '0');
               f = 10 * f;
            }
            if ((frac > MAXFLOAT)||(f > MAXFLOAT)) dcd_error(-24);
            dcd_nextch();
         };
         if (errornum == 0) curconst = curconst + frac/f;
      };
      if ((ch == 'E')||(ch == 'e')||(ch == 'D')||(ch == 'd')) {
         dcd_nextch();
         tenp = 1;
         frac = 0;
         if (ch == '+') {
            dcd_nextch();
         }
         else if (ch == '-') {
            tenp = -tenp;
            dcd_nextch();
         }
         while (isdigit(ch)) {
            if (errornum == 0) frac = 10 * frac + (ch - '0');
            if (frac > MAXLOG) dcd_error(-24);
            dcd_nextch();
         }
         if (errornum == 0) {
            double sumlog = 0.0;
            if (curconst != 0.0) sumlog = log10(fabs(curconst));   /* PJT */
            frac = frac * tenp;
            if ((MINLOG < frac) && (frac < MAXLOG)) {
               sumlog += frac;
               if ((MINLOG < sumlog)&&(sumlog < MAXLOG)) {
                  curconst = curconst * pow(10.0,(double) frac);
               } else dcd_error(-24);
            } else dcd_error(-24);
         }
      }
      sym = hconst;
   } else if (isalpha(ch)) {
      i = 0;
      while ((isalpha(ch)||isdigit(ch)) && (i < maxfunlen)) {
         fun[i++] = toupper(ch);
         dcd_nextch();
      }
      fun[i] = 0;
      curfun = 0;
      while ((curfun < maxfuncts) && strcmp(fun,functs[curfun])) {
         curfun++;
      }
      sym = funct;
      if (curfun == maxfuncts) dcd_error(-12);  /* unknown function */
   } else switch(ch) {
      case '\0' : sym = end; dcd_nextch(); break;
      case '+'  : sym = plus; dcd_nextch(); break;
      case '-'  : sym = minus; dcd_nextch(); break;
      case '*'  : {
         sym = times; dcd_nextch();
         if (ch == '*') {
            sym = power; dcd_nextch();
         }
         break;
      }
      case '/'  : sym = divide; dcd_nextch(); break;
      case '('  : {
         sym = lpar;
         do {
            dcd_nextch();
         } while (ch == ' ');
         break;
      }
      case ')'  : sym = rpar; dcd_nextch(); break;
      case ','  : {
         sym = comma;
         do {
            dcd_nextch();
         } while (ch == ' ');
         if (ch == '\0') dcd_error(-13);     /* syntax error */
         break;
      }
      case '\t' :                        /* make TAB equivalent to SPACE */
      case ' '  : {
         sym = blank;
         while (ch == ' ') dcd_nextch(); 
         if (ch == ')') {
            sym = rpar;
            break;
         }
         if (ch == ',') {
            sym = comma;
            do {
               dcd_nextch();
            } while (ch == ' ');
            if (ch == '\0') dcd_error(-13);     /* syntax error */
            break;
         }
         if (list&&(ch == ']')) {
            sym = rbrac;
            dcd_nextch();
            break;
         }
         if (ch != ':') break;
      }
      case ':'  : {
         sym = colon; 
         do dcd_nextch(); while (ch == ' '); 
         if (ch == '[') dcd_error(-22);       /* nested lists */
         break;
      }
      case '['  : {
         sym = lbrac; 
         do dcd_nextch(); while (ch == ' '); 
         if (list) dcd_error(-22);            /* nested lists */            
         break;
      }
      case ']'  : {
         sym = rbrac;
         dcd_nextch();
         if (!list) dcd_error(-13);           /* syntax error */
         break;
      }
      default   : {
         dcd_error(-14);                      /* illegal character */
         dcd_nextch();
         break;
      }
   }
}

static void dcd_nextwr()                   /* put ilen bytes in output buffer */
{
   int i;
   if (nval++ < nmax) {
      for ( i=0 ; i<ilen ; *dptr++ = unie.b[i++]);
   }
}

static void dcd_movenum()                  /* send next item to output buffer */
{
   if (sym == err) return;
   if (list) dcd_putlist(); else {
      if ((ctyp == 'I')&&(ilen == 2)) {
         if (dval == DCDBLANK) dval = 0.0;
         if ((dval > MINSHORT)&&(dval < MAXSHORT)) unie.s = dcd_round(dval);
         else dcd_error(-19);                             /* conversion error */
      } else if ((ctyp == 'I')&&(ilen == 4)) {
         if (dval == DCDBLANK) dval = 0.0;
         if ((dval > MININT)&&(dval < MAXINT)) unie.i = dcd_round(dval);
         else dcd_error(-19);                             /* conversion error */
      } else if ((ctyp == 'F')&&(ilen == 4)) {
         if (dval == DCDBLANK) setfblank_(&unie.f); else unie.f = (float) dval; 
      }
      else if ((ctyp == 'F')&&(ilen == 8)) {
         if (dval == DCDBLANK) setfblank_(&unie.f); else unie.d = dval;
      }
      if (errorpos == 0) dcd_nextwr();
   }
}

static void dcd_loop()
{
   double a, b, c;
   int    p, q = 0, r;

   if (sym == err) return;
   dcd_expression();
   if (sym == colon) {
      dcd_gencode(hlt);
      dcd_evaluate(q);
      a = dval;
      dcd_nextsym();
      if (sym == colon) {
         dcd_nextsym();
         dcd_expression();
         dcd_gencode(hlt);
         dcd_evaluate(q);
#ifdef BIGLOOP
         if (dval > 0.5) {
#else
         if ((dval > 0.5)&&(dval < MAXSHORT)) {
#endif
            r = dcd_round(dval);
            dval = a;
            for ( p=0; p<r; p++ ) { dcd_movenum();}
         } else dcd_error(-15);               /* wrong repeat/loop argument */
      } else {
         dcd_expression();
         dcd_gencode(hlt);
         dcd_evaluate(q);
         b = dval;
         if (sym == colon) {
            dcd_nextsym();
            dcd_expression();
            dcd_gencode(hlt);
            dcd_evaluate(q);
            c = dval;
         } else { c = 1.0; }
         if (c != 0.0) {
            int    n,i;
            double d;
            d = (b-a)/c;
#ifdef BIGLOOP
            if (d >= 0.0) {	/* 1:1  hence use '>= 0.0' */
#else
            if ((d >= 0.0)&&(d < MAXSHORT)) {	/* 1:1  hence use '>= 0.0' */
#endif
               n = (int) (d+0.00001);
               for ( i=0; i<=n; i++) {
                  dval = a+i*c;
                  dcd_movenum();
               }
            } else dcd_error(-15);
/*
         if (c > 0) {
            for ( dval = a; dval <= b; dval = dval + c) dcd_movenum();
         } else if ( c < 0) {
            for ( dval = a; dval >= b; dval = dval + c) dcd_movenum();
*/
         } else dcd_error(-15);               /* repeat/loop argument 0 */
      }
   } else {
      int i;
      dcd_gencode(hlt);
      if ((lists > 0) && (!list)) {
         for ( i=0; i<listlen[0]; i++ ) {
            dcd_evaluate(i);
            dcd_movenum();
         }
      } else {
         dcd_evaluate(q);
         dcd_movenum();
      }
   }
}

static void dcd_expression()
{
   int s;

   if (sym == err) return;
   dcd_term();
   while ((sym == plus) || (sym == minus)) {
      s = sym;
      dcd_nextsym();
      dcd_term();
      if (s == plus) dcd_gencode(add); else dcd_gencode(sub);
   }
}

static void dcd_term()
{
   int s;

   if (sym == err) return;
   dcd_factor();
   while ((sym == times) || (sym == divide)) {
      s = sym;
      dcd_nextsym();
      dcd_factor();
      if (s == times) dcd_gencode(mul); else dcd_gencode(div);
   }
}

static void dcd_factor()
{
   if (sym == err) return;
   switch (sym) {
      case lpar  : {
         dcd_nextsym(); dcd_expression();
         if (sym == rpar) dcd_nextsym(); else dcd_error(-13);
         break;
      }
      case minus : dcd_nextsym(); dcd_factor(); dcd_gencode(neg); break;
      case plus  : dcd_nextsym(); dcd_factor(); break;
      case hconst: dcd_genconst(curconst); dcd_nextsym(); break;
      case funct : dcd_function(); break;
      case lbrac : {
         dcd_beginlist();
         dcd_list();
         if (sym == rbrac) {
            dcd_endlist();
            dcd_nextsym();
         } else dcd_error(-13);
         break;
      }
      default    : dcd_error(-13); break;
   }
   if (sym == power) {
      dcd_nextsym();
      dcd_factor();
      dcd_gencode(pwr);
   }
}

static void dcd_function()
{
   int f = curfun;
   int n = nargs[curfun];
   
   if (sym == err) return;
   dcd_nextsym();
   if (n > 0) {
      if (sym == lpar) dcd_nextsym(); else dcd_error(-16);
      while (n > 0) {
         dcd_expression();
         if (--n > 0) {
            if (sym == comma) dcd_nextsym(); else dcd_error(-16);
         }
      }
      if (sym == rpar) dcd_nextsym(); else dcd_error(-16);
   }
   dcd_gencode( fie + f );
}

static void dcd_list()
{
   if (sym == err) return;
   do {
      while (ch == ' ') dcd_nextch();
      dcd_nextsym();
      dcd_loop();
   } while ((sym == blank)||(sym == comma));
}

static void dcd_dump()  /* never used */
{
   int    c, o, opc, op;
   if (sym == err) return;
   c = o = 0;
   do {
      if (list) op = lstfiecode[c].opcode[o++]; 
      else      op = fiecode[c].opcode[o++];
      if (o == bid) { c++ ; o = 0; }
      opc = op>fie ? fie : op;
      printf("     %s",mnem[opc]);
      if (opc == fie) {
         printf("   %s",functs[op-opc]);
      } else if (opc == ldc) {
         if (o != 0) c++;
         if (list) printf("   %f",lstfiecode[c++].c);
         else      printf("   %f",fiecode[c++].c);
         o = 0;
      } else if (opc == lst) {
         int i;
         if (o != 0) c++;
         printf("   %d",listlen[0]);
         for ( i=0; i<listlen[0]; i++) {
            printf("\n");
            printf("           %f",fiecode[c++].c);
         };
         o = 0;
      };
      printf("\n");
   } while ((opc != hlt)&&(c < maxfiecode));
}

#define STACKMAX 32

static double stack[STACKMAX];

static int sp;

static void dcd_push(double r)
{
  if (sp==STACKMAX) error("dcd_push: stack exceeded %d",STACKMAX);
  stack[++sp] = r;
}

static double dcd_pop()
{
  if (sp < 0) error("dcd_pop: empty stack");
  return(stack[sp--]);
}

static double dcd_add(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else {
      return(arg1+arg2);
   }
}

static double dcd_sub(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else {
      return(arg1-arg2);
   }
}

static double dcd_mul(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if ((arg1 == 0.0) || (arg2 == 0.0)) {
      return(0.0);
   } else {
      double intsum = (log10(fabs(arg1))+log10(fabs(arg2)));
      if ((MINLOG < intsum) && (intsum < MAXLOG)) {
         return(arg1*arg2);
      } else {
         dcd_error(-17);
         return(DCDBLANK);
      }
   }
}

static double dcd_div(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg2 == 0.0) {
      dcd_error(-17);
      return(DCDBLANK);
   } else if (arg1 == 0.0) {
      return(0.0);
   } else {
      double intsum = (log10(fabs(arg2))-log10(fabs(arg1)));
      if ((MINLOG < intsum)&&(intsum < MAXLOG)) {
         return(arg1/arg2);
      } else {
         dcd_error(-17);
         return(DCDBLANK);
      }
   }
}

static double dcd_neg(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(-arg1);
   }
}

static double dcd_pwr(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 >= 0.0) {
      return(pow(arg1,arg2));
   } else {
      int p = (int) arg2, t;
      double epsilon = 0.000001;
      if (fabs(arg2 - p) <= epsilon) {
         t = (p % 2 == 0) ? 1 : -1;
         return(t * pow(fabs(arg1),arg2));
      } else {
         dcd_error(-17);
         return(DCDBLANK);
      }
   }
}

static double dcd_sin(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(sin(arg1));
   }
}

static double dcd_asin(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) > 1) {
      dcd_error(-17);
   } else {
      return(asin(arg1));
   }
   return(0.0);     /* make stringent compilers happy */
}

static double dcd_sinh(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) > 70) {
      dcd_error(-17);
   } else {
      return(sinh(arg1));
   }
   return(0.0);     /* make stringent compilers happy */
}

static double dcd_cos(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(cos(arg1));
   }
}

static double dcd_acos(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) > 1) {
      dcd_error(-17);
   } else {
      return(acos(arg1));
   }
   return(0.0);     /* make stringent compilers happy */
}

static double dcd_cosh(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) > 70) {
      dcd_error(-17);
   } else {
      return(cosh(arg1));
   }
   return(0.0);     /* make stringent compilers happy */
}

static double dcd_tan(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(tan(arg1));
   }
}

static double dcd_atan(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(atan(arg1));
   }
}
static double dcd_tanh(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) > 70) {
      dcd_error(-17);
   } else {
      return(tanh(arg1));
   }
   return(0.0);     /* make stringent compilers happy */
}

static double dcd_atan2(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else {
      return(atan2(arg1,arg2));
   }
}

static double dcd_rad(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(arg1 * 0.017453292519943295769237);
   }
}

static double dcd_deg(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(arg1 * 57.295779513082320876798155);
   }
}

static double dcd_pi()
{
   double val;
   val = (double) 3.141592653589793238462643;
   return(val);
}

static double dcd_exp(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) > 70) {
      dcd_error(-17);
      return(DCDBLANK);
   } else {
      return(exp(arg1));
   }
}

static double dcd_ln(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (arg1 > 0) {
      return(log(arg1));
   } else {
      dcd_error(-17);
      return(DCDBLANK);
   }
}

static double dcd_log(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (arg1 > 0) {
      return(log10(arg1));
   } else {
      dcd_error(-17);
      return(DCDBLANK);
   }
}

static double dcd_sqrt(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (arg1 < 0) {
      dcd_error(-17);
      return(DCDBLANK);
   } else {
      return(sqrt(arg1));
   }
}

static double dcd_abs(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(fabs(arg1));
   }
}

static double dcd_sinc(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (fabs(arg1) < 1.0e-30) {
      return(1.0);
   } else {
      return(sin(arg1)/arg1);
   }
}

static double dcd_max(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 > arg2) {
      return(arg1);
   } else {
      return(arg2);
   }
}

static double dcd_min(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 < arg2) {
      return(arg1);
   } else {
      return(arg2);
   }
}

static double dcd_erf(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
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
}

static double dcd_erfc(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
      return(1.0 - dcd_erf(arg1));
   }
}

static double dcd_mod(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 == 0.0) {
      dcd_error(-17);
      return(DCDBLANK);
   } else {
     int   xxx = (int) (arg1/arg2);
     return(arg1 - xxx * arg2);
   }
}

static double dcd_int(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
     int xxx = (int) arg1;  /* this could be dangerous */
     return((double) xxx);
   }
}

static double dcd_nint(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else {
     int xxx = (int)(arg1 + 0.5);  /* this could be dangerous */
     return((double) xxx);
   }
}

static double dcd_sign(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (arg1 == 0.0) {
      return(0.0);
   } else if (arg1 > 0.0) {
      return (1.0);
   } else {
      return(-1.0);
   }
}

static double dcd_ifgt(double arg1, double arg2, double arg3, double arg4)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 > arg2) {
      return(arg3);
   } else {
      return(arg4);
   }
}

static double dcd_iflt(double arg1, double arg2, double arg3, double arg4)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 < arg2) {
      return(arg3);
   } else {
      return(arg4);
   }
}

static double dcd_ifge(double arg1, double arg2, double arg3, double arg4)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 >= arg2) {
      return(arg3);
   } else {
      return(arg4);
   }
}

static double dcd_ifle(double arg1, double arg2, double arg3, double arg4)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else if (arg1 <= arg2) {
      return(arg3);
   } else {
      return(arg4);
   }
}

static double dcd_ifeq(double arg1, double arg2, double arg3, double arg4)
{
   if (arg1 == arg2) {
      return(arg3);
   } else {
      return(arg4);
   }
}

static double dcd_ifne(double arg1, double arg2, double arg3, double arg4)
{
   if (arg1 != arg2) {
      return(arg3);
   } else {
      return(arg4);
   }
}

static double dcd_ran()
{
#if NATIVE_RAND
   int xxx = rand();
   return((xxx+1) / 2147483648.0);
#else
   extern double xrandom(double, double);
   return xrandom(0.0,1.0);
#endif
}

static double dcd_ranu(double arg1, double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else {
      return(arg1 + dcd_ran() * (arg2 - arg1));
   }
}

static double dcd_rang(double arg1,double arg2)
{
   if ((arg1 == DCDBLANK) || (arg2 == DCDBLANK)) {
      return(DCDBLANK);
   } else {
     double val, r1, r2;
     r1 = dcd_ran();
     r2 = dcd_ran();
     if (oddran == 0) {
        val = sqrt(-2*log(r1))*cos(6.283185307179586476925286*r2); oddran = 1;
     } else {
        val = sqrt(-2*log(r1))*cos(6.283185307179586476925286*r2); oddran = 0;
     }
     val = arg1 + fabs(arg2) * val;
     return(val);
   }
}

static double dcd_ranp(double arg1)
{
   if (arg1 == DCDBLANK) {
      return(DCDBLANK);
   } else if (arg1 < 0) {
      dcd_error(-17);
      return(DCDBLANK);
   } else {
      double val, cum, p, f;
      if (arg1 < 40) {
         int xxx = (int)(dcd_rang(arg1,sqrt(arg1))+0.5);
         val = xxx;
      } else {
         cum = exp(-arg1);
         p = cum;
         val = 0.0;
         f = dcd_ran();
         while ( f >= cum) {
            val = val + 1.0;
            p = p * arg1 / val;
            cum = cum + p;
         }
      }
     return(val);
   }
}
static void dcd_null(void)
{
  have_null = 1;
  warning("dcd_null: have a null");
}
     
static void dcd_evaluate(int q)
{
   int c, o, opc;
   double arg[maxarg];
   double arg1,arg2;		/* fix order evalution bug PJT */
      
   if (sym == err) return;
   c = o = sp = 0;
   do {
      if (list) {
         opc = lstfiecode[c].opcode[o++];
      } else {
         opc = fiecode[c].opcode[o++];
      }
      if (o == bid) { c++ ; o = 0; };
      if (opc >= fie) {
         int narg = nargs[opc-fie], n;
         for (n = 1 ; n<=narg ; n++ ) arg[narg-n] = dcd_pop();
      };
      switch (opc) {
         case hlt: break;
         case add: dcd_push(dcd_add(dcd_pop(),dcd_pop())); break;
         case sub: 
                arg2 = dcd_pop(); arg1 = dcd_pop();
                dcd_push(dcd_sub(arg1,arg2)); break;
         case mul: dcd_push(dcd_mul(dcd_pop(),dcd_pop())); break;
         case div: 
                arg2 = dcd_pop(); arg1 = dcd_pop();
                dcd_push(dcd_div(arg1,arg2)); break;
         case neg: dcd_push(dcd_neg(dcd_pop())); break;
         case pwr: 
                arg2 = dcd_pop(); arg1 = dcd_pop();
                dcd_push(dcd_pwr(arg1,arg2)); break;
         case ldc: {
            if (o != 0) c++; if (list) {
               dcd_push(lstfiecode[c++].c);
            } else {
               dcd_push(fiecode[c++].c);
            }
            o = 0; break;
         }
         case lst: {
            if (o != 0) c++; o = 0; if (list) {
               dcd_push(lstfiecode[q+c].c);
            } else {
               dcd_push(fiecode[q+c].c);
            }
            c = c + listlen[0]; break;
         }
         default:  switch(opc-fie) {
/* sin   */ case  0: dcd_push(dcd_sin(arg[0])); break;
/* asin  */ case  1: dcd_push(dcd_asin(arg[0])); break;
/* sinh  */ case  2: dcd_push(dcd_sinh(arg[0])); break;
/* cos   */ case  3: dcd_push(dcd_cos(arg[0])); break;
/* acos  */ case  4: dcd_push(dcd_acos(arg[0])); break;
/* cosh  */ case  5: dcd_push(dcd_cosh(arg[0])); break;
/* tan   */ case  6: dcd_push(dcd_tan(arg[0])); break;
/* atan  */ case  7: dcd_push(dcd_atan(arg[0])); break;
/* tanh  */ case  8: dcd_push(dcd_tanh(arg[0])); break;
/* atan2 */ case  9: dcd_push(dcd_atan2(arg[0],arg[1])); break;
/* rad   */ case 10: dcd_push(dcd_rad(arg[0])); break;
/* deg   */ case 11: dcd_push(dcd_deg(arg[0])); break;
/* pi    */ case 12: dcd_push(dcd_pi()); break;
/* exp   */ case 13: dcd_push(dcd_exp(arg[0])); break;
/* ln    */ case 14: dcd_push(dcd_ln(arg[0])); break;
/* log   */ case 15: dcd_push(dcd_log(arg[0])); break;
/* sqrt  */ case 16: dcd_push(dcd_sqrt(arg[0])); break;
/* abs   */ case 17: dcd_push(dcd_abs(arg[0])); break;
/* sinc  */ case 18: dcd_push(dcd_sinc(arg[0])); break;
/* c     */ case 19: dcd_push( 2.997925e+8 ); break;
/* g     */ case 20: dcd_push( 6.6732e-11 ); break;
/* m     */ case 21: dcd_push( 1.99e30 ); break;
/* erf   */ case 22: dcd_push(dcd_erf(arg[0])); break;
/* erfc  */ case 23: dcd_push(dcd_erfc(arg[0])); break;
/* k     */ case 24: dcd_push( 1.380622e-23 ); break;
/* h     */ case 25: dcd_push( 6.6256196e-34 ); break;
/* p     */ case 26: dcd_push( 3.086e16 ); break;
/* s     */ case 27: dcd_push( 5.66961e-8 ); break;
/* max   */ case 28: dcd_push(dcd_max(arg[0],arg[1])); break;
/* min   */ case 29: dcd_push(dcd_min(arg[0],arg[1])); break;
/* mod   */ case 30: dcd_push(dcd_mod(arg[0],arg[1])); break;
/* int   */ case 31: dcd_push(dcd_int(arg[0])); break;
/* nint  */ case 32: dcd_push(dcd_nint(arg[0])); break;
/* sign  */ case 33: dcd_push(dcd_sign(arg[0])); break;
/* undef */ case 34: dcd_push(DCDBLANK); break;
/* ifgt  */ case 35: dcd_push(dcd_ifgt(arg[0],arg[1],arg[2],arg[3])); break;
/* iflt  */ case 36: dcd_push(dcd_iflt(arg[0],arg[1],arg[2],arg[3])); break;
/* ifge  */ case 37: dcd_push(dcd_ifge(arg[0],arg[1],arg[2],arg[3])); break;
/* ifle  */ case 38: dcd_push(dcd_ifle(arg[0],arg[1],arg[2],arg[3])); break;
/* ifeq  */ case 39: dcd_push(dcd_ifeq(arg[0],arg[1],arg[2],arg[3])); break;
/* ifne  */ case 40: dcd_push(dcd_ifne(arg[0],arg[1],arg[2],arg[3])); break;
/* ranu  */ case 41: dcd_push(dcd_ranu(arg[0],arg[1])); break;
/* rang  */ case 42: dcd_push(dcd_rang(arg[0],arg[1])); break;
/* ranp  */ case 43: dcd_push(dcd_ranp(arg[0])); break;
/* sind  */ case 44: dcd_push(sin(dcd_rad(arg[0]))); break;
/* asind */ case 45: dcd_push(dcd_deg(dcd_asin(arg[0]))); break;
/* cosd  */ case 46: dcd_push(cos(dcd_rad(arg[0]))); break;
/* acosd */ case 47: dcd_push(dcd_deg(dcd_acos(arg[0]))); break;
/* tand  */ case 48: dcd_push(tan(dcd_rad(arg[0]))); break;
/* atand */ case 49: dcd_push(dcd_deg(dcd_atan(arg[0]))); break;
/* atand2*/ case 50: dcd_push(dcd_deg(dcd_atan2(arg[0],arg[1]))); break;
/* null  */ case 51: dcd_null(); break;
            default: opc = err; break;
         }; break;
      };
   } while ((opc != hlt) && (opc != err) && (errornum == 0));
   if (opc == err) {
      dcd_error(-17);
   }
   if (errornum != 0) {
      dval = DCDBLANK;
   } else {
      dval = dcd_pop();
   }
   if (list) {
      lstcodeptr = 0;
      lstopcodeptr = 0;
   } else {
      codeptr = 0;
      opcodeptr = 0;
   }
}

/* this is the routine */
void herinp(
    char   *expr,
    int    *nchr,
    char   *type,
    int    *length,
    char   *outv,
    int    *nout,
    int    *nret,
    int    *ierd)
{
   int    i;
   
   cptr = expr;
   mpos = *nchr;
   dptr = outv;
   ilen = *length;
   nmax = *nout;
   ctyp = toupper(*type);
   nval = 0;
   oddran = 0;
   errornum = 0;
   errorpos = 0;
   pos = 0;
   npar = 0;
   codeptr = 0;
   opcodeptr = 0;
   have_null = 0;
   ch = ' ';
   dcd_inifblank();
   switch(ctyp) {
      case 'A': {
         dcd_nextch();
         if (ch == '\0') break;
         do {
            for ( i=0 ; i<ilen; i++ ) {
               *dptr++ = ch;
               dcd_nextch();
            };
            nval++;
         } while ((ch != '\0') && (nval < nmax));
/*       if ((nval == nmax) && (ch != '\0')) dcd_error(-23);  check disabled */
         /* fill with end-of-string */
         if (errornum == 0) {
            int n;
            for ( n = nval; n < nmax; n++) {
               for ( i = 0; i < ilen; i++) *dptr++ = '\0';
            }
         }
         break;
      }
      case 'C': {
         dcd_nextch();
         if (ch == '\0') break;
         do {
            while (ch == ' ') dcd_nextch(); 
            if (ch == ',') dcd_error(-13);
            else {
               for ( i=0 ; i<ilen; i++ ) {
                  if ((ch == '\0') || (ch == ',') || (ch == ' ')) {
                     *dptr++ = ' '; 
                  } else {
                     *dptr++ = ch;
                     dcd_nextch();
                  }
               }
               if ((ch != '\0') && (ch != ',') && (ch != ' ')) dcd_error(-13);
               else {
                  nval++;
                  while (ch == ' ') dcd_nextch();
                  if (ch == ',') {
                     do {
                        dcd_nextch();
                     } while (ch == ' ');
                     if (ch == '\0') dcd_error(-13);
                  }
               }
           };
         } while ((ch != '\0') && (errornum == 0) && (nval < nmax));
         if ((nval == nmax) && ( ch != '\0')) dcd_error(-23);
         /* fill with end-of-string */
         if (errornum == 0) {
            int n;
            for ( n = nval; n < nmax; n++) {
               for ( i = 0; i < ilen; i++) {*dptr++ = ' ';}
            }
         }
         break;
      }
      case 'L': {
         char logc[maxboollen];
         int  i, curlog;
         dcd_nextch();
         if (ch == '\0') break;
         do {
            while (ch == ' ') dcd_nextch();
            if (isalpha(ch)) {
               i = 0;
               while ((isalpha(ch)) && (i < maxboollen)) {
                  logc[i++] = toupper(ch);
                  dcd_nextch();
               };
               curlog = 0;
               while (( curlog < maxbools) && strncmp(logc,bools[curlog],i)) {
                  curlog++;
               }
               if (curlog == maxbools) dcd_error(-13); /* syntax error */
               else {
                  switch(ilen) {
                     case 1 : unie.b[0] = boolv[curlog]; break;
                     case 2 : unie.s = boolv[curlog]; break;
                     case 4 : unie.i = boolv[curlog]; break;
                     default: break;
                  };
                  dcd_nextwr();
                  while (ch == ' ') dcd_nextch();
                  if (ch == ',') {
                     do {
                        dcd_nextch();
                     } while (ch == ' ');
                     if (ch == '\0') dcd_error(-13);
                  }
               }
            } else dcd_error(-13);
         } while ((ch != '\0') && (errornum == 0) && (nval < nmax));
         if ((nval == nmax) && (ch != '\0')) dcd_error(-23);
         break;
      }
      case 'I':
      case 'F': {
         dcd_nextch();
         if (ch == '\0') break;
         sym = blank;
         do {
            dcd_inilist();
            while (ch == ' ') dcd_nextch(); 
            if (ch != '\0') {
               dcd_nextsym();
               dcd_loop();
            } else sym = end;
         } while ((errornum == 0)&&((sym == blank)||(sym == comma)));
         if ((errornum == 0) && (ch != '\0')) dcd_error(-13);
         if ((errornum == 0) && (nval > nmax)) dcd_error(-23);
         break;
      }
      default : errornum = -11;                  /* bad call */
   }
   *nret = nval;
   *ierd = errornum;
}

#if NEED_DCD

/* for integers */
int dcdint_(
    char    *expr,
    int *outv, int *nout, int *ierd, int nchr)
{
   int nret, length = sizeof(int);
   char    type = 'I';
   herinp(expr,&nchr,&type,&length,outv,nout,&nret,ierd);   /* BAD_PARM */
   return(nret);
}

/* for logicals */
int dcdlog_( 
    char    *expr,
    int     *outv, int *nout, int *ierd, int nchr)
{
   int     nret, length = sizeof(logical);
   char    type = 'L';
   herinp(expr,&nchr,&type,&length,outv,nout,&nret,ierd);   /* BAD_PARM */
   return(nret);
}
 
/* for reals */
int dcdreal_( 
    char    *expr,
    int *outv, int *nout, int *ierd, int nchr)
{
   int     nret, length = sizeof(real);
   char    type = 'F';
   herinp(expr,&nchr,&type,&length,outv,nout,&nret,ierd);   /* BAD_PARM */
   return(nret);
}
 

/* for double precision */
int dcddble_( 
    char    *expr,
    int *outv, int *nout, int *ierd, int nchr)
{
   int     nret, length = sizeof(double);
   char    type = 'F';
   herinp(expr,&nchr,&type,&length,outv,nout,&nret,ierd);   /* BAD_PARM */
   return(nret);
}
 
/* for characters */
int dcdchar_( 
    char    *expr, char *outv,
    int     *nout, int *ierd, int nchr, int nlen)
{
   integer nret, length = nlen;
   char    type = 'C';
   herinp(expr,&nchr,&type,&length,outv,nout,&nret,ierd);	/* OK_PARM */
   return(nret);
}

#endif


#if defined(TESTBED)

char *defv[] = {
    "expr=\n      Expression to parse",
    "VERSION=1\n  11-nov-2003 PJT",
    NULL,
};

nemo_main()
{
   double outv[256];
   char *expr;
   int    carg;
   int    nout = 256, nret, ierd;
   char   type = 'f';

   expr = getparam("expr");
   nret = dcddble_(expr,outv,&nout,&ierd,strlen(expr));
   switch(ierd) {
      case -11: printf("bad call\n"); break;
      case -12: printf("unknown function\n"); break;
      case -13: printf("syntax error\n"); break;
      case -14: printf("illegal character\n"); break;
      case -15: printf("wrong repeat argument\n"); break;
      case -16: printf("wrong number of arguments\n"); break;
      case -17: printf("arithmetic error\n"); break;
      case -18: printf("not enough internal memory\n"); break;
      case -19: printf("conversion error\n"); break;
      case -20: printf("unequal list length\n"); break;
      case -21: printf("empty list\n"); break;
      case -22: printf("nested lists\n"); break;
      case -23: printf("output buffer overflow\n"); break;
      case -24: printf("floating overflow/underflow in conversion\n"); break;
      default: {
         int n;
         int p;
         for (n = 0; n < nret; n++) {
#if 0         	
            if (fblank_(&outv[n])) printf("BLANK\n"); else {
#else
            if (have_null) printf("NULL\n"); else {
#endif            	
               if (outv[n] != 0.0) {
                  p = log10(fabs(outv[n]));
                  outv[n] = outv[n]*pow(10.0,(double) -p);
                  if (p > 9) printf("  %.10fE+%d\n",outv[n],p);
                  else if (p >=0) printf("  %.10fE+0%d\n",outv[n],p);
                  else if (p > -9) printf("  %.10fE-0%d\n",outv[n],-p);
                  else printf("  %.10fE-%d\n",outv[n],-p);
               } else printf("  0.0000000000E+00\n");
            }
         }
      }
      break;
   }
}
#endif
