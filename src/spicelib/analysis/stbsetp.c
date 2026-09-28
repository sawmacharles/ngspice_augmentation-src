/**********
ngspice_augmentation private fork -- card #456 (plan T3.2/T3.3), 2026-09-28.
Native `.stb` loop-stability analysis: parameter table. See stbdefs.h.
**********/

#include "ngspice/ngspice.h"
#include "ngspice/ifsim.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/stbdefs.h"
#include "ngspice/cktdefs.h"

#include "analysis.h"

/* ARGSUSED */
int
STBsetParm(CKTcircuit *ckt, JOB *anal, int which, IFvalue *value)
{
    STBAN *job = (STBAN *) anal;

    NG_IGNORE(ckt);

    switch (which) {

    case STB_START:
        if (value->rValue <= 0.0) {
            errMsg = copy("Frequency of <= 0 is invalid for STB start");
            job->STBstartFreq = 1.0;
            return (E_PARMVAL);
        }
        job->STBstartFreq = value->rValue;
        break;

    case STB_STOP:
        if (value->rValue <= 0.0) {
            errMsg = copy("Frequency of <= 0 is invalid for STB stop");
            job->STBstopFreq = 1.0;
            return (E_PARMVAL);
        }
        job->STBstopFreq = value->rValue;
        break;

    case STB_STEPS:
        job->STBnumberSteps = value->iValue;
        break;

    case STB_DEC:
        if (value->iValue) {
            job->STBstepType = DECADE;
        } else if (job->STBstepType == DECADE) {
            job->STBstepType = 0;
        }
        break;

    case STB_OCT:
        if (value->iValue) {
            job->STBstepType = OCTAVE;
        } else if (job->STBstepType == OCTAVE) {
            job->STBstepType = 0;
        }
        break;

    case STB_LIN:
        if (value->iValue) {
            job->STBstepType = LINEAR;
        } else if (job->STBstepType == LINEAR) {
            job->STBstepType = 0;
        }
        break;

    case STB_PROBE:
        job->STBprobeName = value->uValue;
        break;

    case STB_LOCALGND:
        job->STBlocalGnd = value->nValue;
        job->STBlocalGndGiven = 1;
        break;

    default:
        return (E_BADPARM);
    }
    return (OK);
}

static IFparm STBparms[] = {
    { "start",    STB_START,    IF_SET | IF_ASK | IF_REAL,     "starting frequency" },
    { "stop",     STB_STOP,     IF_SET | IF_ASK | IF_REAL,     "ending frequency" },
    { "numsteps", STB_STEPS,    IF_SET | IF_ASK | IF_INTEGER,  "number of frequencies" },
    { "dec",      STB_DEC,      IF_SET | IF_ASK | IF_FLAG,     "step by decades" },
    { "oct",      STB_OCT,      IF_SET | IF_ASK | IF_FLAG,     "step by octaves" },
    { "lin",      STB_LIN,      IF_SET | IF_ASK | IF_FLAG,     "step linearly" },
    { "probe",    STB_PROBE,    IF_SET | IF_ASK | IF_INSTANCE, "loop-break probe: an existing 0V V-source" },
    { "localgnd", STB_LOCALGND, IF_SET | IF_ASK | IF_NODE,     "local ground reference node (optional)" },
};

SPICEanalysis STBinfo = {
    {
        "STB",
        "Loop stability (Tian two-injection return ratio) analysis",

        NUMELEMS(STBparms),
        STBparms
    },
    sizeof(STBAN),
    FREQUENCYDOMAIN,
    1,
    STBsetParm,
    STBaskQuest,
    NULL,
    STBan
};
