/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1987 Wayne A. Christopher, U. C. Berkeley CAD Group
**********/

/*
 * All the functions used in the B-source parse tree.  These functions return HUGE
 * if their argument is out of range.
 */

#include "ngspice/ngspice.h"
#include <stdio.h>
#include "ngspice/fteext.h"
#include "ngspice/ifsim.h"
#include "ngspice/inpptree.h"
#include "ngspice/cktdefs.h"
#include "inpxx.h"
#include "ngspice/compatmode.h"
#include "ngspice/cpextern.h"

double PTfudge_factor;

#define MODULUS(NUM,LIMIT) ((NUM) - ((int) ((NUM) / (LIMIT))) * (LIMIT))

double
PTabs(double arg)
{
    return fabs(arg);
}

double
PTsgn(double arg)
{
    return arg > 0.0 ? 1.0 : arg < 0.0 ? -1.0 : 0.0;
}

double
PTplus(double arg1, double arg2)
{
    return (arg1 + arg2);
}

double
PTminus(double arg1, double arg2)
{
    return (arg1 - arg2);
}

double
PTtimes(double arg1, double arg2)
{
    return (arg1 * arg2);
}

double
PTdivide(double arg1, double arg2)
{
    if (arg2 >= 0.0)
        arg2 += PTfudge_factor;
    else
        arg2 -= PTfudge_factor;

    if (arg2 == 0.0)
        return (HUGE);

    return (arg1 / arg2);
}

double
PTpower(double arg1, double arg2)
{
    double res;

    if (arg1 == 0)
        res = 0;
    else if(arg1 > 0)
        res = pow(arg1, arg2);
    else {
        /* If arg2 is quasi an integer, round it to have pow not fail
            when arg1 is negative. Takes into account the double
            representation which sometimes differs in the last digit(s). */
        if (AlmostEqualUlps(nearbyint(arg2), arg2, 10))
            res = pow(arg1, round(arg2));
        else
            /* As per LTSPICE specification for ** */
            res = 0;
    }
    return res;
}

double
PTpowerH(double arg1, double arg2)
{
    double res;

    if (newcompat.hs) {
        if (arg1 < 0)
            res = pow(arg1, round(arg2));
        else if (arg1 == 0){
            res = 0;
        }
        else
        {
            res = pow(arg1, arg2);
        }
    }
    else {
        if (arg1 >= 0)
            res = pow(arg1, arg2);
        else {
            /* If arg2 is quasi an integer, round it to have pow not fail
               when arg1 is negative. Takes into account the double
               representation which sometimes differs in the last digit(s). */
            if (AlmostEqualUlps(nearbyint(arg2), arg2, 10))
                res = pow(arg1, round(arg2));
            else
                /* As per LTSPICE specification for ** */
                res = 0;
        }
    }
    return res;
}


double
PTpwr(double arg1, double arg2)
{
    /* if PSPICE device is evaluated */
    if (arg1 == 0.0 && arg2 < 0.0 && newcompat.ps)
        arg1 += PTfudge_factor;

    if (arg1 < 0.0)
        return (-pow(-arg1, arg2));
    else
        return (pow(arg1, arg2));
}

double
PTmin(double arg1, double arg2)
{
    return arg1 > arg2 ? arg2 : arg1;
}

double
PTmax(double arg1, double arg2)
{
    return arg1 > arg2 ? arg1 : arg2;
}

double
PTacos(double arg)
{
    return (acos(arg));
}

double
PTacosh(double arg)
{
    return (acosh(arg));
}

double
PTasin(double arg)
{
    return (asin(arg));
}

double
PTasinh(double arg)
{
    return (asinh(arg));
}

double
PTatan(double arg)
{
    return (atan(arg));
}

double
PTatanh(double arg)
{
    return (atanh(arg));
}

double
PTustep(double arg)
{
    if (arg < 0.0)
        return 0.0;
    else if (arg > 0.0)
        return 1.0;
    else
        return 0.5; /* Ick! */
}

/* MW. PTcif is like "C" if - 0 for (arg<=0), 1 elsewhere */

double
PTustep2(double arg)
{
    if (arg <= 0.0)
        return 0.0;
    else if (arg <= 1.0)
        return arg;
    else /* if (arg > 1.0) */
        return 1.0;
}

double
PTeq0(double arg)
{
    return (arg == 0.0) ? 1.0 : 0.0;
}

double
PTne0(double arg)
{
    return (arg != 0.0) ? 1.0 : 0.0;
}

double
PTgt0(double arg)
{
    return (arg > 0.0) ? 1.0 : 0.0;
}

double
PTlt0(double arg)
{
    return (arg < 0.0) ? 1.0 : 0.0;
}

double
PTge0(double arg)
{
    return (arg >= 0.0) ? 1.0 : 0.0;
}

double
PTle0(double arg)
{
    return (arg <= 0.0) ? 1.0 : 0.0;
}

double
PTuramp(double arg)
{
    if (arg < 0.0)
        return 0.0;
    else
        return arg;
}

double
PTcos(double arg)
{
    return (cos(MODULUS(arg, 2 * M_PI)));
}

double
PTcosh(double arg)
{
    return (cosh(arg));
}

/* Limit the exp: If arg > EXPARGMAX (arbitrarily selected to 14), continue with linear output,
   if compatmode PSPICE is selected.
   If arg exceeds 227.9559242, output its exp value 1e99. */
double
PTexp(double arg)
{
    if (newcompat.ps && arg > EXPARGMAX)
        return EXPMAX * (arg - EXPARGMAX + 1.);
    else if (arg > 227.9559242)
        return 1e99;
    else
        return (exp(arg));
}

/* If arg < , returning HUGE will lead to an error message.
   If arg == 0, don't bail out, but return an arbitrarily very negative value (-1e99).
   Arg 0 may happen, when starting iteration for op or dc simulation. */
double
PTlog(double arg)
{
    if (arg < 0.0)
        return (HUGE);
    if (arg == 0)
        return -1e99;
    return (log(arg));
}

double
PTlog10(double arg)
{
    if (arg < 0.0)
        return (HUGE);
    if (arg == 0)
        return -1e99;
    return (log10(arg));
}

double
PTsin(double arg)
{
    return (sin(MODULUS(arg, 2 * M_PI)));
}

double
PTsinh(double arg)
{
    return (sinh(arg));
}

double
PTsqrt(double arg)
{
    if (arg < 0.0)
        return (HUGE);
    return (sqrt(arg));
}

double
PTtan(double arg)
{
    return (tan(MODULUS(arg, M_PI)));
}

double
PTtanh(double arg)
{
    return (tanh(arg));
}

double
PTuminus(double arg)
{
    return (- arg);
}

double
PTpwl(double arg, void *data)
{
  struct pwldata { int n; double *vals; } *thing = (struct pwldata *) data;

  double y;

  int k0 = 0;
  int k1 = thing->n/2 - 1;

  /* monotonically increasing abscissa */
  if (thing->vals[0] < thing->vals[2]) {
      while (k1 - k0 > 1) {
          int k = (k0 + k1) / 2;
          if (thing->vals[2 * k] > arg)
              k1 = k;
          else
              k0 = k;
      }
  }
  /* monotonically decreasing abscissa */
  else {
      while (k1 - k0 > 1) {
          int k = (k0 + k1) / 2;
          if (thing->vals[2 * k] < arg)
              k1 = k;
          else
              k0 = k;
      }
  }
  /* interpolate the ordinate */
  y = thing->vals[2*k0+1] +
    (thing->vals[2*k1+1] - thing->vals[2*k0+1]) *
    (arg - thing->vals[2*k0]) / (thing->vals[2*k1] - thing->vals[2*k0]);

  return y;
}

double
PTpwl_derivative(double arg, void *data)
{
  struct pwldata { int n; double *vals; } *thing = (struct pwldata *) data;

  double y;

  int k0 = 0;
  int k1 = thing->n/2 - 1;

  while(k1-k0 > 1) {
    int k = (k0+k1)/2;
    if(thing->vals[2*k] > arg)
      k1 = k;
    else
      k0 = k;
  }

  y =
    (thing->vals[2*k1+1] - thing->vals[2*k0+1]) /
    (thing->vals[2*k1]   - thing->vals[2*k0]);

  return y;
}

double
PTceil(double arg1)
{
    return (ceil(arg1));
}

double
PTfloor(double arg1)
{
    return (floor(arg1));
}

double
PTnint(double arg1)
{
    /* round to "nearest integer",
     *   round half-integers to the nearest even integer
     *   rely on default rounding mode of IEEE 754 to do so
     */
    return nearbyint(arg1);
}


/* ddt() until s158 (set ddtlegacy): evaluated once per new time point, at its first
   iteration, as (y_n - y_n-1) / (t_n-1 - t_n-2).
   If time == 0, return 0.
   If not transient sim, return 0.
   The derivative is then (y2-y1)/(t2-t1).
   */
static double
PTddt_legacy(double arg, void* data)
{
    struct ddtdata { int n; double* vals; } *thing = (struct ddtdata*)data;
    double y, time;

    CKTcircuit* ckt = ft_curckt->ci_ckt;

    time = ckt->CKTtime;

    if (time == 0) {
        thing->vals[3] = arg;
        return 0;
    }

    if (!(ckt->CKTmode & MODETRAN))
        return 0;

    if (time > thing->vals[0]) {
        thing->vals[4] = thing->vals[2];
        thing->vals[5] = thing->vals[3];
        thing->vals[2] = thing->vals[0];
        thing->vals[3] = thing->vals[1];
        thing->vals[0] = time;
        thing->vals[1] = arg;

/*      // Some less effective smoothing option
        if (thing->vals[2] > 0) {
            thing->vals[6] = 0.5 * ((arg - thing->vals[3]) / (time - thing->vals[2]) + thing->vals[6]);
        }
*/
        if (thing->n > 1) {
            thing->vals[6] = (thing->vals[1] - thing->vals[3]) / (thing->vals[2] - thing->vals[4]);
        }
        else {
            thing->vals[6] = 0;
            thing->vals[3] = arg;
        }
        thing->n += 1;
    }

    y = thing->vals[6];

    return y;
}


/* ddt(): the derivative of the converged waveform, by the circuit's own
   integration formula (the one capacitors use, NIintegrate), so that
   C * ddt(v(a,b)) equals the current of a capacitor C between a and b.

   vals[0..6] is the state of ddt() as until s158 (PTddt_legacy, selected by
   vals[DDT_LEGACY]); this one uses
     vals[8]  time of the latest evaluation       vals[9]  its argument
     vals[10] last accepted time                   vals[19] latest result
     vals[11..17] arguments at the last 7 accepted points (newest first)
     vals[18] result at the last accepted point    vals[20] accepted points held
   A transient time point is t = t_acc + CKTdelta. When the time changes and
   the latest evaluation was made at t - CKTdelta, that point was accepted
   and its last (converged) argument and result enter the history; a
   rejected point leaves it alone. The result is computed at every
   iteration, with its Jacobian CKTag[0] (PTddt_derivative). */

#define DDT_T    8
#define DDT_X    9
#define DDT_TACC 10
#define DDT_XH   11   /* 11 .. 17 */
#define DDT_YACC 18
#define DDT_Y    19
#define DDT_NH   20

static void
ddt_history(double *v, CKTcircuit *ckt)
{
    double t = ckt->CKTtime;

    if (t != v[DDT_T]) {
        double tacc = t - ckt->CKTdelta;
        if (fabs(v[DDT_T] - tacc) <= 1e-6 * ckt->CKTdelta) {
            int i;
            for (i = 6; i > 0; i--)
                v[DDT_XH + i] = v[DDT_XH + i - 1];
            v[DDT_XH] = v[DDT_X];
            v[DDT_YACC] = v[DDT_Y];
            v[DDT_TACC] = v[DDT_T];
            if (v[DDT_NH] < 7)
                v[DDT_NH] += 1;
        }
    }
}

/* result and Jacobian at the present time point; returns 0 when the history
   does not reach back to this point's previous accepted one */
static int
ddt_eval(double *v, CKTcircuit *ckt, double arg, double *y, double *dy)
{
    double *ag = ckt->CKTag;
    int order = ckt->CKTorder, i;

    ddt_history(v, ckt);
    if (v[DDT_NH] < 1 || fabs(v[DDT_TACC] - (ckt->CKTtime - ckt->CKTdelta)) > 1e-6 * ckt->CKTdelta)
        return 0;

    if (ckt->CKTintegrateMethod == GEAR && order > 1 && v[DDT_NH] >= order) {
        *y = ag[0] * arg;
        for (i = 1; i <= order; i++)
            *y += ag[i] * v[DDT_XH + i - 1];
        *dy = ag[0];
    } else if (ckt->CKTintegrateMethod == TRAPEZOIDAL && order == 2) {
        *y = ag[0] * (arg - v[DDT_XH]) - ag[1] * v[DDT_YACC];
        *dy = ag[0];
    } else {        /* first order (also gear without enough history yet) */
        *y = (arg - v[DDT_XH]) / ckt->CKTdelta;
        *dy = 1 / ckt->CKTdelta;
    }
    return 1;
}

double
PTddt(double arg, void* data)
{
    struct ddtdata { int n; double* vals; } *thing = (struct ddtdata*)data;
    double *v = thing->vals;
    CKTcircuit* ckt = ft_curckt->ci_ckt;
    double t, y, dy;

    t = ckt->CKTtime;

    /* operating point, and the first point of a transient: the start of the
       history; `.option ddtlegacy' or `set ddtlegacy' selects ddt() as
       until s158 */
    if (t == 0 || !(ckt->CKTmode & MODETRAN)) {
        v[DDT_LEGACY] = cp_getvar("ddtlegacy", CP_BOOL, NULL, 0) ? 1 : 0;
        if (v[DDT_LEGACY] != 0)
            return PTddt_legacy(arg, data);
        v[DDT_T] = t;
        v[DDT_X] = arg;
        v[DDT_Y] = 0;
        v[DDT_NH] = 0;
        return 0;
    }

    if (v[DDT_LEGACY] != 0)
        return PTddt_legacy(arg, data);

    if (!ddt_eval(v, ckt, arg, &y, &dy))
        y = v[DDT_Y];   /* time went back past the history (a periodic
                           analysis restarting its period): keep the latest
                           result until a point is accepted again */

    v[DDT_T] = t;
    v[DDT_X] = arg;
    v[DDT_Y] = y;
    return y;
}

/* d ddt(arg) / d arg, for the Jacobian */
double
PTddt_derivative(double arg, void* data)
{
    struct ddtdata { int n; double* vals; } *thing = (struct ddtdata*)data;
    double *v = thing->vals;
    CKTcircuit* ckt = ft_curckt->ci_ckt;
    double y, dy;

    if (v[DDT_LEGACY] != 0)
        return 0;
    if (ckt->CKTtime == 0 || !(ckt->CKTmode & MODETRAN))
        return 0;
    return ddt_eval(v, ckt, arg, &y, &dy) ? dy : 0;
}
