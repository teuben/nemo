#include "catutil.h"
#define TRUE  1
#define FALSE 0
int readqn(qnstr, iqn, n)
char *qnstr;
short *iqn;
int n;
{ /* read n quanta from string to iqn*/
  char *pqn;
  int i, ich, ic, ival;

  pqn = qnstr;
  for (i = 0; i < n; ++i) {
    if (*pqn == 0) break;
    ich = (*pqn & 0xff); ++pqn;
    ic  = (*pqn & 0xff); ++pqn;
    if (ic == ' '){
      iqn[i] = 0; continue;
    }
    ival = ic - '0';
    if (ival < 0 || ival > 9) break;
    if (ich != ' ') {
      if (ich == '-') {
        ival = -ival;
      } else if (ich >= '0' && ich <= '9') {
        ival += (ich - '0') * 10;
      } else if (ich >= 'a' && ich <= 'z') {
        ival = -ival - 10 * (ich - ('a' - 1));
      } else if (ich >= 'A' && ich <= 'Z') {
        ival += 10 * (ich - ('A' - 10));
      } else {
        break;
      }
    }
    iqn[i] = (short)ival;
  }
  return (n - i);
} /* readqn */

int pcard(card, val, nval, fmtlen)
char *card;
double *val;
int nval;
int *fmtlen;
{ /*  parses a card for NVAL numbers */
  /*  ndec = -2 for comma detected, ndec = -1 indicates looking for  */
  /*         beginning of field, ndec >= 0 indicates in field */
  /*  pwrflg = 0 indicates mantissa decimal point not detected */
  /*  pwrflg > 0 indicates mantissa decimal point detected */
  /*  pwrflg = -1 indicates e, E, d, or D found */
  /*  pwrflg = -2 indicates in exponent field */

  static double pten[]={1.,10.,100.,1000.,1e4,1e5,1e6,1e7,1e8,1e9};
  static int mxdec = 9;
  char *newfield;
  double tmp, fac;
  int ndec, kval, itmp, ipwr, pwrflg, ich, npwr;
  int newnum, neg, negpwr; /* boolean */

  if (sizeof(int) == 2) mxdec = 4; /* allow for 16-bit integer */
  newfield = (char *)fmtlen;
  if (newfield) newfield = card + fmtlen[0];
  newnum = TRUE; neg = negpwr = FALSE;
  itmp= ipwr= pwrflg= kval= npwr= 0; tmp = fac = 0.;
  ndec= -1;
  while (kval < nval) {
    ich = (*card & 0xff);
    if (ich != 0) {
      if (card == newfield){
        ich = ','; 
      } else { 
        card++;
      }
    }
    if (ich >= '0' && ich <= '9') { /*  character is a number */
      if (pwrflg == 0){ 
        ++npwr; /* count integer part of mantissa */
      } else if (pwrflg > 0){ 
        --ipwr; /* count fraction part of mantissa */
      }else {/* pwrflg < 0 */
        pwrflg = -2; /* flag indicates digit found in exponent */
      }
      if (ndec <= 0) {
        itmp = ich - '0'; ndec = 1;  /* first digit */
      } else if (ndec == mxdec) {
        if (pwrflg < 0) return kval; /* exponent field is too big */
        /* now is a good time to convert integer to real */
        if (newnum) {
          tmp = (double) itmp; newnum = FALSE;
        } else {
          tmp = tmp * pten[ndec] + itmp;
        }
        itmp = ich - '0'; ndec = 1;
      } else { /* accumulate integer */
        ++ndec;
        itmp = itmp * 10 + ich - '0';
      }
    } else if (ich == '.' && pwrflg == 0) { /* first decimal point */
      if (ndec < 0) ndec = 0;
      pwrflg = 1;
    } else if (ich == '-' && pwrflg == -1) { /* leading - in exponent field */
      negpwr = TRUE; ndec = 0; pwrflg = -2;
    } else if (ich == '+' && pwrflg == -1) { /* leading + in exponent field */
      ndec = 0; pwrflg = -2;
    } else {/* character is not a number or decimal point, E+, E- */
      if (ndec >= 0) { /* save results from number decoding */
        if (pwrflg < 0) {/* integer follows 'E' */
          pwrflg = 0;
          if (negpwr) itmp = -itmp;
          npwr += itmp; ipwr += itmp;
        } else {/* finish up mantissa */
          pwrflg = 0;
          if (newnum) {
            tmp = (double) itmp;
          } else {
            tmp = itmp + tmp * pten[ndec];
          }
          if (ich == 'E' || ich == 'e' || ich == 'D' || ich == 'd') {
            /*  look for exponent */
            pwrflg = -1; ndec = 0; negpwr = FALSE;
          }
        }
        if (pwrflg == 0) { /* number finished */
          if (npwr < -38) {
            ipwr = 0; tmp = 0.;
          } else if (npwr > 37) {
            ipwr = 0; tmp = 1.e+37;
          }
          if (ipwr != 0) { /* scale by powers of 10 */
            if (ipwr < 0) {
              ipwr = -ipwr;
              if (ipwr > 7) fac = 1. / pten[8];
              itmp = ipwr & 7;
              if (itmp != 0) tmp /= pten[itmp];
            } else {
              fac = pten[8];
              itmp = ipwr & 7;
              if (itmp != 0) tmp *= pten[itmp];
            }
            if (ipwr > 7){
              ipwr = ipwr >> 3;
              if ((ipwr & 1) != 0) tmp *= fac;
              while (ipwr > 1) {
                ipwr = ipwr >> 1; fac = fac * fac;
                if ((ipwr & 1) != 0) tmp *= fac; 
              }
            }
            ipwr = 0;
          } /* end scale */
          if (neg) tmp = -tmp;
          val[kval] = tmp;
          ++kval;
          if (newfield && kval < nval) newfield += fmtlen[kval];
          ndec= -1; npwr = 0; neg = FALSE;
        }
        itmp = 0; newnum = TRUE;
      } /* finished save results */ 
      if (ich == 0) return kval;
      /* check for delimiters in new field */
      if (ich == '.') { /* decimal afer exponent field, assume new */
        ndec = 0; pwrflg = 1;
      } else if (ich == '-') { /* leading - in mantissa field */
        ndec = 0; neg = TRUE; 
      } else if (ich == '+') { /* leading + in mantissa field */
        ndec = 0; neg = FALSE;
      } else if (ich == ',') {
        if (ndec == -2) { /* second comma */
          ++kval;
          if (newfield && kval < nval) newfield += fmtlen[kval];
        }
        ndec = -2;
      } else if (ich == '/') { /* end of line character */
        kval = nval;
      }
    }
  }
  return kval;
} /* pcard */
