/**********
ngspice_augmentation private fork -- card #456 (plan T3.2/T3.3), 2026-09-28.
Native `.stb` loop-stability analysis: parameter readback. See stbdefs.h.
**********/

#include "ngspice/ngspice.h"
#include "ngspice/ifsim.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/stbdefs.h"
#include "ngspice/cktdefs.h"

/* ARGSUSED */
int
STBaskQuest(CKTcircuit *ckt, JOB *anal, int which, IFvalue *value)
{
    STBAN *job = (STBAN *) anal;

    NG_IGNORE(ckt);

    switch (which) {

    case STB_START:
        value->rValue = job->STBstartFreq;
        break;

    case STB_STOP:
        value->rValue = job->STBstopFreq;
        break;

    case STB_STEPS:
        value->iValue = job->STBnumberSteps;
        break;

    case STB_DEC:
        value->iValue = (job->STBstepType == DECADE);
        break;

    case STB_OCT:
        value->iValue = (job->STBstepType == OCTAVE);
        break;

    case STB_LIN:
        value->iValue = (job->STBstepType == LINEAR);
        break;

    case STB_PROBE:
        value->uValue = job->STBprobeName;
        break;

    case STB_LOCALGND:
        value->nValue = job->STBlocalGnd;
        break;

    case STB_PROBE2:
        value->uValue = job->STBprobeName2;
        break;

    case STB_MODE:
        value->iValue = job->STBmode;
        break;

    case STB_SENS:
        /* card #475: readback is informational only -- the actual list
         * lives in STBsensNames/STBsensAll, not reachable through the
         * single-scalar IFvalue this callback returns; report whether ANY
         * sens= was given, and how many explicit names (0 if sens=all). */
        value->iValue = job->STBsensAll ? -1 : job->STBsensCount;
        break;

    default:
        return (E_BADPARM);
    }
    return (OK);
}
