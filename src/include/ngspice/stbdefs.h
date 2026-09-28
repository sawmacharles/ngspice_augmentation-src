/**********
ngspice_augmentation private fork -- card #456 (plan T3.2/T3.3), 2026-09-28.

Native `.stb` loop-stability analysis (Tian two-injection return ratio).
See research/cadence-stb-pstb.md secs 2.4/3b and research/ngspice-internals-pss-pac.md
sec 2.6 in the ngspice_augmentation repo for the design. Structurally modeled on
acdefs.h/ACAN (frequency sweep bookkeeping) and, for the probe/branch pattern, on
the unmerged 2020 loop_analysis branch's loopdefs.h/LOOPAN (not merged upstream;
this is an independent, smaller implementation, not a port of that file).
**********/

#ifndef ngspice_STBDEFS_H
#define ngspice_STBDEFS_H

#include "ngspice/acdefs.h"   /* reuses the DECADE/OCTAVE/LINEAR step-type enum */
#include "ngspice/jobdefs.h"

/* structure used to describe an STB analysis to be performed */
typedef struct {
    int JOBtype;
    JOB *JOBnextJob;     /* pointer to next thing to do */
    IFuid JOBname;       /* name of this job */
    double STBstartFreq;
    double STBstopFreq;
    double STBfreqDelta;  /* multiplier for decade/octave stepping, step for linear */
    double STBsaveFreq;   /* frequency at which we left off last time (pause/resume) */
    int STBstepType;      /* DECADE / OCTAVE / LINEAR, see acdefs.h */
    int STBnumberSteps;
    IFuid STBprobeName;    /* name of the probe V-source instance (probe=) */
    CKTnode *STBlocalGnd;  /* optional local-ground reference node (localgnd=) */
    unsigned int STBlocalGndGiven : 1;
} STBAN;

enum {
    STB_START = 1,
    STB_STOP,
    STB_STEPS,
    STB_DEC,
    STB_OCT,
    STB_LIN,
    STB_PROBE,
    STB_LOCALGND,
};

#endif
