#ifndef ngspice_MIF_H
#define ngspice_MIF_H

/* ===========================================================================
FILE    MIF.h

MEMBER OF process XSPICE

Copyright 1991
Georgia Tech Research Corporation
Atlanta, Georgia 30332
All Rights Reserved

PROJECT A-8503

AUTHORS

    9/12/91  Bill Kuhn

MODIFICATIONS

    <date> <person name> <nature of modifications>

SUMMARY

    This file structure definitions global data used with the MIF package.
    The global data structure is used to circumvent the need to modify
    argument lists in existing SPICE 3C1 functions.

INTERFACES

    None.

REFERENCED FILES

    None.

NON-STANDARD FEATURES

    None.

=========================================================================== */

#include  "ngspice/miftypes.h"
#include  "ngspice/mifdefs.h"
#include  "ngspice/cktdefs.h"


typedef struct {
   Mif_Boolean_t   init;        /* TRUE if first call to model */
   Mif_Boolean_t   anal_init;   /* TRUE if first call for this analysis type */
   Mif_Analysis_t  anal_type;   /* The type of analysis being performed */
   Mif_Call_Type_t call_type;   /* Type of call to code model - analog or event-driven */
   double          evt_step;    /* The current DC step or time in event analysis */
} Mif_Circuit_Info_t;


typedef struct {
   double           current;    /* The current dynamic breakpoint time */
   double           last;       /* The last used dynamic breakpoint time */
} Mif_Bkpt_Info_t;


typedef struct {
    Mif_Boolean_t   global;     /* Set by .option to force all models to use auto */
    Mif_Boolean_t   local;      /* Set by individual model to request auto partials */
} Mif_Auto_Partial_t;


typedef struct {
   Mif_Circuit_Info_t circuit;    /* Circuit data that will be needed by MIFload */
   MIFinstance        *instance;  /* Current instance struct */
   CKTcircuit         *ckt;       /* The ckt struct for the circuit */
   char               *errmsg;    /* An error msg from a cm_... function */
   Mif_Bkpt_Info_t    breakpoint; /* Data used by dynamic breakpoints */
   Mif_Auto_Partial_t auto_partial; /* Flags to enable auto partial computations */
} Mif_Info_t;


/* These are defined in mif.c */
extern int MIFiSize;
extern int MIFmSize;


extern Mif_Info_t      g_mif_info;
/* card #551: model-time offset (s) of the current run: a `.pss method=shoot` period run restarts the circuit clock at 0 while the
 * XSPICE models (event queue, stored change times, oscillator phases) must see a CONTINUING time; 0 except inside such an analysis. */
extern double          g_xs_toff;
extern double          g_xs_ag0_fix;   /* card #557: when nonzero cm_analog_integrate() uses it instead of CKTag[0] (td small-signal C extraction) */
extern Mif_Private_t   g_mif_noise_cm_data;

/* card #714: digital edge-timing tracking for the td small-signal analyses (pacshoot.c).  While g_xs_dig_on is set, every queued
 * output event carries a tag (0: its time does not depend on any analog quantity; > 0: a timing record of pacshoot.c), every event
 * node remembers the tag and time of its last change, and an instance called because inputs changed passes the tag of those inputs to
 * its outputs.  g_xs_dig_hook (set by pacshoot.c) is told about input edges (what 0), decides the tag of each changed output (what 1)
 * and makes a record for simultaneous input edges of different origin (what 2).  Defined in evtload.c. */
extern int     g_xs_dig_on;
extern int     g_xs_dig_cur_tag;
extern int     g_xs_dig_nnode;
extern int    *g_xs_dig_node_tag;
extern double *g_xs_dig_node_t;
extern int   (*g_xs_dig_hook)(void *ckt, void *inst, int what, int conn, int port, int intag, double t);
void XSdigReset(int num_nodes);
extern int    *g_xs_dig_carry_tag;   /* [num_outputs] tag of an output event put back from the XSPICE snapshot (card #714) */
extern int     g_xs_dig_ncarry;
int  XSdigPendingTag(void *ckt, int (*allowed)(int tag));

/* card #714: exact location of adc_bridge threshold crossings (evtsetup.c), called from dctran.c */
extern int XSadcLocateActive;
int  XSadcLocateWanted(void);
int  XSadcLocate(void *ckt, double *newdelta);
void XSadcAccept(void *ckt);
extern int XSperiodicRun;   /* set by the shooting / td recording runs (pssshoot.c, pacshoot.c) */

#endif
