/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1988 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

#include "ngspice/ngspice.h"
#include <stdio.h>
#include <string.h>
#include "ngspice/ifsim.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/inpdefs.h"
#include "ngspice/inpmacs.h"
#include "ngspice/fteext.h"
#include "inpxx.h"
#include "ngspice/cpdefs.h"
#include "ngspice/tskdefs.h"
#include "ngspice/stringutil.h"
/* card #457: probe=<inst>:<terminal> auto-insert needs the CKT device/node
 * manipulation API (CKTfndDev, CKTmkVolt, CKTmodCrt, CKTcrtElt, CKTbindNode,
 * CKTnum2nod, CKTtypelook, CKTpName) and the device table (DEVices[]). */
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "ngspice/gendefs.h"
#include "ngspice/stbdefs.h"   /* STB_MODE_DM / STB_MODE_CM */

static int
dot_noise(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
          TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int which;			/* which analysis we are performing */
    int error;			/* error code temporary */
    char *name;			/* the resistor's name */
    char *nname1;		/* the first node's name */
    char *nname2;		/* the second node's name */
    CKTnode *node1;		/* the first node's node pointer */
    CKTnode *node2;		/* the second node's node pointer */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *steptype;		/* ac analysis, type of stepping function */

    int found;
    char *point;

    /* .noise V(OUTPUT,REF) SRC {DEC OCT LIN} NP FSTART FSTOP <PTSPRSUM> */
    which = ft_find_analysis("NOISE");
    if (which == -1) {
        LITERR("Noise analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Noise Analysis", &foo, task));
    INPgetTok(&line, &name, 1);

    /* Make sure the ".noise" command is followed by V(xxxx).  If it
       is, extract 'xxxx'.  If not, report an error. */

    if (name != NULL) {

        if ((*name == 'V' || *name == 'v') && !name[1]) {

            INPgetNetTok(&line, &nname1, 0);
            INPtermInsert(ckt, &nname1, tab, &node1);
            ptemp.nValue = node1;
            GCA(INPapName, (ckt, which, foo, "output", &ptemp));

            if (*line != ')') {
                INPgetNetTok(&line, &nname2, 1);
                INPtermInsert(ckt, &nname2, tab, &node2);
                ptemp.nValue = node2;
            } else {
                ptemp.nValue = gnode;
            }
            GCA(INPapName, (ckt, which, foo, "outputref", &ptemp));

            tfree(name);
            INPgetTok(&line, &name, 1);
            INPinsert(&name, tab);
            ptemp.uValue = name;
            GCA(INPapName, (ckt, which, foo, "input", &ptemp));

            INPgetTok(&line, &steptype, 1);
            ptemp.iValue = 1;
            error = INPapName(ckt, which, foo, steptype, &ptemp);
            tfree(steptype);
            if (error)
                current->error = INPerrCat(current->error, INPerror(error));
            parm = INPgetValue(ckt, &line, IF_INTEGER, tab);
            error = INPapName(ckt, which, foo, "numsteps", parm);
            if (error)
                current->error = INPerrCat(current->error, INPerror(error));
            parm = INPgetValue(ckt, &line, IF_REAL, tab);
            error = INPapName(ckt, which, foo, "start", parm);
            if (error)
                current->error = INPerrCat(current->error, INPerror(error));
            parm = INPgetValue(ckt, &line, IF_REAL, tab);
            error = INPapName(ckt, which, foo, "stop", parm);
            if (error)
                current->error = INPerrCat(current->error, INPerror(error));

            /* now see if "ptspersum" has been specified by the user */

            for (found = 0, point = line; (!found) && (*point != '\0'); found = ((*point != ' ') && (*(point++) != '\t')))
                ;
            if (found) {
                parm = INPgetValue(ckt, &line, IF_INTEGER, tab);
                error = INPapName(ckt, which, foo, "ptspersum", parm);
                if (error)
                    current->error = INPerrCat(current->error, INPerror(error));
            } else {
                ptemp.iValue = 0;
                error = INPapName(ckt, which, foo, "ptspersum", &ptemp);
                if (error)
                    current->error = INPerrCat(current->error, INPerror(error));
            }
        } else
            LITERR("bad syntax "
                   "[.noise v(OUT) SRC {DEC OCT LIN} "
                   "NP FSTART FSTOP <PTSPRSUM>]\n");
    } else {
        LITERR("bad syntax "
               "[.noise v(OUT) SRC {DEC OCT LIN} "
               "NP FSTART FSTOP <PTSPRSUM>]\n");
    }
    return 0;
}


static int
dot_op(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
       TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int which;			/* which analysis we are performing */
    int error;			/* error code temporary */

    NG_IGNORE(line);
    NG_IGNORE(tab);
    NG_IGNORE(gnode);

    /* .op */
    which = ft_find_analysis("OP");
    if (which == -1) {
        LITERR("DC operating point analysis unsupported\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Operating Point", &foo, task));
    return (0);
}


static int
dot_disto(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
          TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int which;			/* which analysis we are performing */
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *steptype;		/* ac analysis, type of stepping function */

    NG_IGNORE(gnode);

    /* .disto {DEC OCT LIN} NP FSTART FSTOP <F2OVERF1> */
    which = ft_find_analysis("DISTO");
    if (which == -1) {
        LITERR("Small signal distortion analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Distortion Analysis", &foo, task));
    INPgetTok(&line, &steptype, 1);	/* get DEC, OCT, or LIN */
    ptemp.iValue = 1;
    GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);	/* number of points */
    GCA(INPapName, (ckt, which, foo, "numsteps", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* fstart */
    GCA(INPapName, (ckt, which, foo, "start", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* fstop */
    GCA(INPapName, (ckt, which, foo, "stop", parm));
    if (*line) {
        parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* f1phase */
        GCA(INPapName, (ckt, which, foo, "f2overf1", parm));
    }
    return (0);
}


static int
dot_ac(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
       TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */
    char *steptype;		/* ac analysis, type of stepping function */
    bool pdef = FALSE;  /* issue a warning if default parameters are used */
    char* mline = line; /* for debug printout */
    double startval, stopval;

    NG_IGNORE(gnode);

    /* .ac {DEC OCT LIN} NP FSTART FSTOP */
    which = ft_find_analysis("AC");
    if (which == -1) {
        LITERR("AC small signal analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "AC Analysis", &foo, task));
    INPgetTok(&line, &steptype, 1);	/* get DEC, OCT, or LIN */
    if (!*steptype || (!ciprefix("dec", steptype) && !ciprefix("oct", steptype) && !ciprefix("lin", steptype))) {
        LITERR("Missing DEC, OCT, or LIN.\n");
        return (0);
    }
    ptemp.iValue = 1;
    GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
    tfree(steptype);

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab); /* number of points */
    if (parm->iValue < 1) {
        pdef = TRUE;
        parm->iValue = 10;
    }
    GCA(INPapName, (ckt, which, foo, "numsteps", parm));

    if((*line != '.' && !isdigit(*line)) || (*line == '.' && !isdigit(*(line+1))))
        pdef = TRUE;
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* fstart */
    startval = parm->rValue;
    if (startval < 0) {
        pdef = TRUE;
        startval = parm->rValue = 1.;
    }
    GCA(INPapName, (ckt, which, foo, "start", parm));

    if((*line != '.' && !isdigit(*line)) || (*line == '.' && !isdigit(*(line+1))))
        pdef = TRUE;
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* fstop */
    stopval = parm->rValue;
    if (stopval < startval) {
        pdef = TRUE;
        parm->rValue = 1000. * startval;
    }
    GCA(INPapName, (ckt, which, foo, "stop", parm));

    if (pdef) {
        fprintf(stderr, "Warning, ngspice assumes default parameter(s) for ac simulation\n");
        fprintf(stderr, "    Check your input line '.ac %s'\n\n", mline);
    }
    return (0);
}

/* ngspice_augmentation card #457: probe=<inst>:<terminal> auto-insert.
 *
 * Splits the named terminal of an existing device instance onto a new
 * internal node, in series with a freshly-inserted 0V V-source, so a loop
 * can be probed without hand-editing the netlist to add a probe device.
 * Modeled on (not copied from, and not compiled from) the unmerged
 * loop_analysis branch's LOOPpreset (src/spicelib/analysis/loopan.c on
 * branch loop_analysis, commit fb93986e, fetched read-only as
 * origin/loop_analysis per card #457's brief -- research/
 * ngspice-internals-pss-pac.md sec 2.6). LOOPpreset does the equivalent
 * splice from an analysis's an_init hook, which in this fork's cktdojob.c
 * flow runs AFTER CKTsetup has already sized/factored the matrix once for
 * the run -- inserting a device there would need the (KLU-capable, per
 * buildrecipe/README.md) sparse matrix to grow post-setup, which was not
 * risked here. Doing the splice HERE instead, at netlist-parse time (both
 * ".stb ... probe=inst:term" on a netlist card, and the interactive "stb
 * ... probe=inst:term" command -- if_run() in spiceif.c builds the latter
 * into a ".stb ..." line and feeds it through this SAME dot_stb, so both
 * forms hit this code), is unconditionally before CKTsetup: CKTsetup runs
 * once per CKTdoJob call, itself only reachable after the deck (or the
 * interactive command's synthesized one-line deck) has been fully
 * translate()d, i.e. after dot_stb has already returned. This means: this
 * instance must already exist in the netlist by the time the .stb card
 * naming it is reached (normal SPICE deck order: analysis cards after
 * device cards); a forward reference is COULD-NOT-RUN, not silently
 * ignored (see the "instance not found" case below).
 *
 * a/b convention (documented per card #457's requirement -- matches
 * stban.c's existing probe=<Vname> doc and LOOPpreset's own node1/node2
 * assignment): the new V-source's node 1 (a, "downstream": where run 1's
 * voltage is read back) is the NEW internal node, now bound to the named
 * device's terminal; node 2 (b, "upstream": where run 2's current is
 * injected) is the ORIGINAL node the terminal used to connect straight to.
 * Unlike LOOPpreset (whose node-name suffix is the constant "loopinj",
 * regardless of which terminal), the suffix here includes the terminal
 * name, so probing two different terminals of the SAME instance in one
 * deck cannot collide on one synthesized node name.
 *
 * Idempotent like LOOPpreset: probing the same inst:terminal a second
 * time (e.g. re-running the interactive "stb" command in one session)
 * finds the already-inserted probe by its deterministic name and reuses
 * it rather than splicing again.
 *
 * Red-proof (card #457): a nonexistent instance or terminal name attaches
 * an error to the card via LITERR (current->error) -- the standard ngspice
 * parse-error path, which refuses to run the circuit -- never a silent
 * default probe. */
static int
stb_probe_splice(CKTcircuit *ckt, struct card *current, char *spec, IFuid *out_probe_uid)
{
    char *colon, *instname, *termname;
    GENinstance *inst;
    IFdevice *dev;
    int termidx, ntermnames;
    IFuid eltUid, modUid;
    GENinstance *probesrc;
    GENmodel *modvsrc;
    CKTnode *nodeinj, *extnode;
    IFvalue ptemp;
    int error, vtype;
    char sbuf[128];
    char msg[256];

    *out_probe_uid = NULL;

    colon = strchr(spec, ':');
    if (!colon)
        return 1;  /* not this syntax; caller treats spec as a plain instance name */

    instname = spec;
    *colon = '\0';
    termname = colon + 1;

    if (!*instname || !*termname) {
        current->error = INPerrCat(current->error,
            INPmkTemp("'.stb' probe=inst:terminal needs both an instance name and a terminal name\n"));
        return 2;
    }

    inst = CKTfndDev(ckt, (IFuid) instname);
    if (!inst || inst->GENmodPtr->GENmodType < 0) {
        snprintf(msg, sizeof msg, "stb probe: instance '%s' not in circuit (forward references "
                 "to a probe=inst:terminal are not supported -- the instance must appear "
                 "before the .stb card)\n", instname);
        current->error = INPerrCat(current->error, INPmkTemp(msg));
        return 2;
    }

    dev = &(DEVices[inst->GENmodPtr->GENmodType]->DEVpublic);
    ntermnames = dev->numNames ? *(dev->numNames) : 0;
    for (termidx = 0; termidx < ntermnames; termidx++)
        if (strcasecmp(termname, dev->termNames[termidx]) == 0)
            break;
    if (termidx >= ntermnames) {
        snprintf(msg, sizeof msg, "stb probe: instance '%s' has no terminal '%s'\n", instname, termname);
        current->error = INPerrCat(current->error, INPmkTemp(msg));
        return 2;
    }

    /* idempotent probe name, keyed on instance + terminal index */
    snprintf(sbuf, sizeof sbuf, "stbprobe_%s", termname);
    error = SPfrontEnd->IFnewUid(ckt, &eltUid, inst->GENname, sbuf, UID_INSTANCE, NULL);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not allocate probe instance name\n"));
        return 2;
    }
    probesrc = CKTfndDev(ckt, eltUid);
    if (probesrc && probesrc->GENmodPtr->GENmodType == CKTtypelook("Vsource")) {
        /* already spliced in (repeat interactive "stb ... probe=inst:term") */
        *out_probe_uid = probesrc->GENname;
        return 0;
    }

    extnode = CKTnum2nod(ckt, GENnode(inst)[termidx]);

    error = CKTmkVolt(ckt, &nodeinj, inst->GENname, sbuf);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not create the internal node\n"));
        return 2;
    }

    error = SPfrontEnd->IFnewUid(ckt, &modUid, inst->GENname, "stbprobemod", UID_MODEL, NULL);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not allocate probe model name\n"));
        return 2;
    }
    vtype = CKTtypelook("Vsource");
    modvsrc = NULL;
    error = CKTmodCrt(ckt, vtype, &modvsrc, modUid);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not create the probe's model\n"));
        return 2;
    }
    error = CKTcrtElt(ckt, modvsrc, &probesrc, eltUid);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not create the probe instance\n"));
        return 2;
    }
    ptemp.rValue = 0;
    error = CKTpName("dc", &ptemp, ckt, vtype, "probe", &probesrc);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not set the probe's DC value\n"));
        return 2;
    }
    error = CKTbindNode(ckt, probesrc, 1, nodeinj);
    if (!error) error = CKTbindNode(ckt, probesrc, 2, extnode);
    if (!error) error = CKTbindNode(ckt, inst, termidx + 1, nodeinj);
    if (error) {
        current->error = INPerrCat(current->error, INPmkTemp("stb probe: could not wire the probe in\n"));
        return 2;
    }

    printf("stb: probe=%s:%s -- auto-inserted 0V probe '%s' "
           "(a=new internal node bound to %s's terminal %s, b=%s, the original net)\n",
           instname, termname, (char *) probesrc->GENname, instname, termname,
           extnode->name ? (char *) extnode->name : "0");

    *out_probe_uid = probesrc->GENname;
    return 0;
}

/* ngspice_augmentation card #456/#457: native .stb loop-stability analysis.
 * Syntax: .stb dec|oct|lin <n> <fstart> <fstop> probe=<probespec> [probe2=<Vn>] [localgnd=<node>] [mode=dm|cm]
 *   <probespec> ::= <Vname>            | an existing 0V V-source already in the loop (#456)
 *                  | <inst>:<terminal> | auto-insert a probe by splitting that terminal (#457,
 *                                        stb_probe_splice() above)
 *                  | <Vp>,<Vn>         | native DM/CM (#457, requires mode=dm|cm) -- best-effort;
 *                                        may not survive the interactive command's own comma-
 *                                        splitting tokenizer, see the probe2= case below
 * Native DM/CM (#457, requires mode=dm|cm together with either probe=<Vp>,<Vn> or the
 * RECOMMENDED probe=<Vp> probe2=<Vn>): Vp and Vn are each an existing 0V V-source, driven in
 * lockstep (antiphase for dm, in phase for cm); T is read from Vp only (branch p), see stban.c
 * and stb/diffprobe.lib's "reading branch 1 alone is enough" symmetry note.
 * (spec.md stb_syntax; see src/spicelib/analysis/stban.c for the a/b terminal convention).
 * The dec/oct/lin/n/fstart/fstop prefix is parsed exactly like .ac
 * (dot_ac above); the probe=/localgnd=/mode= keyword tail is parsed here (a
 * minimal, .stb-specific key=value reader -- not the generic dot_params()
 * helper from the unmerged loop_analysis branch, which this fork does not
 * carry).
 *
 * QUOTING (found while validating card #456): on a *netlist* ".stb" card,
 * the numparam brace-expression preprocessor (src/frontend/numparam/
 * xpressn.c) runs over every dot-card line before dot_stb ever sees it,
 * and tries to evaluate the RHS of any bare (unquoted) "key=value" as a
 * number -- so "probe=Vprobe" fails with "Undefined parameter [vprobe]"
 * unless quoted: "probe="Vprobe"". This is not specific to .stb -- the
 * unmerged loop_analysis branch's README.loop documents the same
 * requirement for its own string-valued keyword params ("insrc="<src>":
 * ... must be quoted inside \"\""). INPgetTok() (used below) already
 * strips the quotes, so both quoted and unquoted forms parse the SAME way
 * once past this preprocessing stage. The interactive "stb" control-block
 * command (com_stb, runcoms.c) is not run through this preprocessor, so
 * probe=Vprobe (unquoted) works there either way. See
 * tests/test_stb_native.py and the card #456 report for both forms
 * exercised against the same build. */
static int
dot_stb(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
        TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int error;
    IFvalue ptemp;
    IFvalue *parm;
    int which;
    char *steptype;
    bool pdef = FALSE;
    char *mline = line;
    double startval, stopval;
    bool sawProbe = FALSE;
    bool sawProbe2 = FALSE, sawMode = FALSE;

    NG_IGNORE(gnode);

    which = ft_find_analysis("STB");
    if (which == -1) {
        LITERR("STB loop-stability analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "STB Loop Stability Analysis", &foo, task));

    INPgetTok(&line, &steptype, 1);
    if (!*steptype || (!ciprefix("dec", steptype) && !ciprefix("oct", steptype) && !ciprefix("lin", steptype))) {
        LITERR("Missing DEC, OCT, or LIN.\n");
        tfree(steptype);
        return (0);
    }
    ptemp.iValue = 1;
    GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
    tfree(steptype);

    /* card #466: n/fstart/fstop are read the same way .ac reads them (see
     * dot_ac above), but .ac's own "clamp any bad value to a default,
     * print only the generic pdef warning" behaviour is wrong for .stb --
     * it silently swallowed the E_PARMVAL check STBsetParm (stbsetp.c)
     * already has for fstart/fstop <= 0, so that check's own message
     * ("Frequency of <= 0 is invalid for STB start/stop") could never
     * appear: by the time GCA(INPapName, ...) reached STBsetParm, the bad
     * value had already been overwritten with a valid default (found
     * 2026-09-29 while checking the independent reader's report, card
     * #466's brief). The fix: only fill in a default when the field was
     * genuinely OMITTED -- i.e. INPgetValue/INPevaluate did not advance
     * `line` at all, so no token was there to read -- and otherwise pass
     * whatever the user actually typed straight through to the setter, so
     * an explicitly-bad value (<=0, or fstop < fstart, or n <= 0) is
     * REFUSED (E_PARMVAL, no analysis, no vectors -- see if_run()'s
     * "if (deck.error) return 2" for the interactive form and inp.c's
     * per-card error loop for the netlist-card form; both run off the
     * same current->error/deck.error this GCA() macro fills in via
     * INPerror(), which prefers the analysis's own errMsg when set) rather
     * than silently swapped for a default. STBsetParm's STB_STOP case now
     * also checks fstop < fstart (job->STBstartFreq, already set by the
     * "start" GCA() call above by the time "stop" runs) since that
     * ordering check needs both values at once and has nowhere else
     * analysis-specific to live; STB_STEPS gained the same "n <= 0"
     * E_PARMVAL check numsteps never had. */
    {
        char *before = line;
        parm = INPgetValue(ckt, &line, IF_INTEGER, tab);
        if (line == before) {
            /* omitted -- default, as before */
            pdef = TRUE;
            parm->iValue = 10;
        }
        /* else: pass the value through as-is, even if <= 0 --
         * STBsetParm's STB_STEPS case now refuses it. */
        GCA(INPapName, (ckt, which, foo, "numsteps", parm));
        if (error)
            return (0);
    }

    {
        char *before = line;
        parm = INPgetValue(ckt, &line, IF_REAL, tab);
        startval = parm->rValue;
        if (line == before) {
            pdef = TRUE;
            startval = parm->rValue = 1.;
        }
        /* else: pass startval through as-is, even if <= 0 --
         * STBsetParm's STB_START case refuses it (E_PARMVAL). */
        GCA(INPapName, (ckt, which, foo, "start", parm));
        if (error)
            return (0);
    }

    {
        char *before = line;
        parm = INPgetValue(ckt, &line, IF_REAL, tab);
        stopval = parm->rValue;
        if (line == before) {
            pdef = TRUE;
            stopval = parm->rValue = 1000. * startval;
        }
        /* else: pass stopval through as-is, even if <= 0 or < startval --
         * STBsetParm's STB_STOP case refuses both (E_PARMVAL). */
        GCA(INPapName, (ckt, which, foo, "stop", parm));
        if (error)
            return (0);
    }

    if (pdef) {
        fprintf(stderr, "Warning, ngspice assumes default parameter(s) for stb simulation\n");
        fprintf(stderr, "    Check your input line '.stb %s'\n\n", mline);
    }

    /* tail: probe=<Vname> [localgnd=<node>], any order */
    while (*line) {
        char *key;
        while (*line == ' ' || *line == '\t' || *line == '+')
            line++;
        if (!*line)
            break;
        INPgetTok(&line, &key, 0);
        if (!*key) { tfree(key); break; }
        if (*line != '=') {
            fprintf(stderr, "Error: expected key=value on .stb, got '%s'\n", key);
            LITERR("bad syntax [.stb dec|oct|lin n fstart fstop probe=<Vname> [localgnd=<node>]]\n");
            tfree(key);
            return (0);
        }
        line++; /* skip '=' */
        if (cieq(key, "probe")) {
            /* card #457: probe= now takes three forms --
             *   probe=<Vname>              plain existing 0V V-source (unchanged, card #456)
             *   probe=<inst>:<terminal>    auto-insert (stb_probe_splice() above)
             *   probe=<Vp>,<Vn> mode=dm|cm native differential/common-mode (branch p is Vp)
             * Two tokenizer quirks (inpgtok.c) both matter here, so this uses
             * INPgetNetTok(), not INPgetTok():
             *   (1) INPgetTok() splits tokens AT a comma (and gobbles it when
             *       gobble=1), so "Vp,Vn" cannot be told apart from plain "Vp"
             *       AFTER calling it once -- the comma is looked for in the
             *       raw line FIRST, below, specifically to decide this.
             *   (2) INPgetTok() also treats a bare '+'/'-' as a possible
             *       floating-point exponent sign and can end a token on one
             *       (e.g. terminal name "C+" -- a capacitor's IFdevice
             *       termNames[0], see stb_probe_splice() -- reads back as
             *       just "C", found while validating card #457's
             *       probe=Cn:C+). INPgetNetTok() (used for localgnd= below
             *       already) has no such sign-state logic. */
            char *name, *name2 = NULL;
            char *p;
            bool hascomma = FALSE;
            for (p = line; *p && *p != ' ' && *p != '\t' && *p != '\r'; p++)
                if (*p == ',') { hascomma = TRUE; break; }

            INPgetNetTok(&line, &name, 1);
            if (hascomma) {
                INPgetNetTok(&line, &name2, 1);
                if (!*name || !*name2) {
                    LITERR("'.stb' probe=<Vp>,<Vn> needs two non-empty probe names\n");
                    tfree(name); tfree(name2);
                    return (0);
                }
                INPinsert(&name, tab);
                INPinsert(&name2, tab);
                ptemp.uValue = name;
                GCA(INPapName, (ckt, which, foo, "probe", &ptemp));
                ptemp.uValue = name2;
                GCA(INPapName, (ckt, which, foo, "probe2", &ptemp));
                sawProbe = TRUE;
                sawProbe2 = TRUE;
            } else if (strchr(name, ':')) {
                IFuid spliced = NULL;
                int serr = stb_probe_splice(ckt, current, name, &spliced);
                tfree(name);
                if (serr) return (0);  /* error already attached to current->error */
                ptemp.uValue = spliced;
                GCA(INPapName, (ckt, which, foo, "probe", &ptemp));
                sawProbe = TRUE;
            } else {
                INPinsert(&name, tab);
                ptemp.uValue = name;
                GCA(INPapName, (ckt, which, foo, "probe", &ptemp));
                sawProbe = TRUE;
            }
        } else if (cieq(key, "probe2")) {
            /* card #457: probe2=<Vn>, an explicit alternative to
             * probe=<Vp>,<Vn> for native mode=dm|cm. Documented as the
             * RECOMMENDED form (spec.md stb_syntax / this card's report):
             * the interactive "stb ..." control command is tokenized into
             * a wordlist by ngspice's OWN command-line lexer (cp_event
             * loop) BEFORE if_run() ever reassembles it into a ".stb ..."
             * line for dot_stb -- and that lexer ALSO splits on comma, so
             * probe=Vp,Vn arrives here already mangled into "probe=vp" +
             * a stray "vn" token (found validating card #457; a netlist
             * ".stb" card's numparam preprocessing is a separate stage and
             * was not exhaustively checked either). probe=<Vp> probe2=
             * <Vn> mode=dm|cm, three independent key=value pairs, needs no
             * character that any layer of ngspice's tokenizing treats as a
             * separator, so it is reliable on both the interactive command
             * and the netlist card. The comma form above is left in for
             * whichever form it happens to survive on, but is not the
             * form this card's own tests (or its docs) rely on. */
            char *name2;
            INPgetNetTok(&line, &name2, 1);
            if (!*name2) {
                LITERR("'.stb' probe2= needs a probe name\n");
                tfree(name2);
                return (0);
            }
            INPinsert(&name2, tab);
            ptemp.uValue = name2;
            GCA(INPapName, (ckt, which, foo, "probe2", &ptemp));
            sawProbe2 = TRUE;
        } else if (cieq(key, "mode")) {
            char *mval;
            INPgetTok(&line, &mval, 1);
            if (cieq(mval, "dm")) {
                ptemp.iValue = STB_MODE_DM;
            } else if (cieq(mval, "cm")) {
                ptemp.iValue = STB_MODE_CM;
            } else {
                LITERR("'.stb' mode= must be dm or cm\n");
                tfree(mval);
                return (0);
            }
            GCA(INPapName, (ckt, which, foo, "mode", &ptemp));
            tfree(mval);
            sawMode = TRUE;
        } else if (cieq(key, "localgnd")) {
            char *nname;
            CKTnode *nnode;
            INPgetNetTok(&line, &nname, 1);
            INPtermInsert(ckt, &nname, tab, &nnode);
            ptemp.nValue = nnode;
            GCA(INPapName, (ckt, which, foo, "localgnd", &ptemp));
        } else {
            fprintf(stderr, "Error: unknown parameter '%s' on .stb - ignored\n", key);
            INPgetTok(&line, &mline, 1); /* consume the value so parsing can continue */
            tfree(mline);
        }
        tfree(key);
    }
    if (!sawProbe) {
        LITERR("'.stb' needs probe=<Vname> (an existing 0V V-source that breaks the loop), "
               "probe=<inst>:<terminal>, or probe=<Vp>,<Vn> mode=dm|cm\n");
    }
    if (sawProbe2 != sawMode) {
        LITERR("'.stb' probe=<Vp>,<Vn> and mode=dm|cm must be given together\n");
    }
    return (0);
}

static int
dot_pz(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
       TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */
    char *steptype;		/* ac analysis, type of stepping function */

    NG_IGNORE(gnode);

    /* .pz nodeI nodeG nodeJ nodeK {V I} {POL ZER PZ} */
    which = ft_find_analysis("PZ");
    if (which == -1) {
        LITERR("Pole-zero analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Pole-Zero Analysis", &foo, task));
    parm = INPgetValue(ckt, &line, IF_NODE, tab);
    GCA(INPapName, (ckt, which, foo, "nodei", parm));
    parm = INPgetValue(ckt, &line, IF_NODE, tab);
    GCA(INPapName, (ckt, which, foo, "nodeg", parm));
    parm = INPgetValue(ckt, &line, IF_NODE, tab);
    GCA(INPapName, (ckt, which, foo, "nodej", parm));
    parm = INPgetValue(ckt, &line, IF_NODE, tab);
    GCA(INPapName, (ckt, which, foo, "nodek", parm));
    INPgetTok(&line, &steptype, 1);	/* get V or I */
    ptemp.iValue = 1;
    GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
    INPgetTok(&line, &steptype, 1);	/* get POL, ZER, or PZ */
    ptemp.iValue = 1;
    GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
    return (0);
}


static int
dot_dc(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
       TSKtask *task, CKTnode *gnode, JOB *foo)
{
    char *name;			/* the resistor's name */
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */

    NG_IGNORE(gnode);

    /* .dc SRC1NAME Vstart1 Vstop1 Vinc1 [SRC2NAME Vstart2 Vstop2 Vinc2]
       Return 1 upon error because of bad syntax (missing tokens).*/
    which = ft_find_analysis("DC");
    if (which == -1) {
        LITERR("DC transfer curve analysis unsupported\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "DC transfer characteristic", &foo, task));
    INPgetTok(&line, &name, 1);
    if (*name == '\0')
        return 1;
    INPinsert(&name, tab);
    ptemp.uValue = name;
    GCA(INPapName, (ckt, which, foo, "name1", &ptemp));
    if (*line == '\0')
        return 1;
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* vstart1 */
    GCA(INPapName, (ckt, which, foo, "start1", parm));
    if (*line == '\0')
        return 1;
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* vstop1 */
    GCA(INPapName, (ckt, which, foo, "stop1", parm));
    if (*line == '\0')
        return 1;
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* vinc1 */
        if (parm->rValue == 0)
            return 1;
    GCA(INPapName, (ckt, which, foo, "step1", parm));
    if (*line) {
        INPgetTok(&line, &name, 1);
        if (*line == '\0')
            return 1;
        INPinsert(&name, tab);
        ptemp.uValue = name;
        GCA(INPapName, (ckt, which, foo, "name2", &ptemp));
        parm = INPgetValue(ckt, &line, IF_REAL, tab); /* vstart2 */
        if (*line == '\0')
            return 1;
        GCA(INPapName, (ckt, which, foo, "start2", parm));
        parm = INPgetValue(ckt, &line, IF_REAL, tab); /* vstop2 */
        if (*line == '\0')
            return 1;
        GCA(INPapName, (ckt, which, foo, "stop2", parm));
        parm = INPgetValue(ckt, &line, IF_REAL, tab); /* vinc2 */
        if (parm->rValue == 0)
            return 1;
        GCA(INPapName, (ckt, which, foo, "step2", parm));
    }
    return 0;
}


static int
dot_tf(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
       TSKtask *task, CKTnode *gnode, JOB *foo)
{
    char *name;			/* the resistor's name */
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    int which;			/* which analysis we are performing */
    char *nname1;		/* the first node's name */
    char *nname2;		/* the second node's name */
    CKTnode *node1;		/* the first node's node pointer */
    CKTnode *node2;		/* the second node's node pointer */

    /* .tf v( node1, node2 ) src */
    /* .tf vsrc2             src */
    which = ft_find_analysis("TF");
    if (which == -1) {
        LITERR("Transfer Function analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Transfer Function", &foo, task));
    INPgetTok(&line, &name, 0);
    /* name is now either V or I or a serious error */
    if (*name == 'v' && strlen(name) == 1) {
        if (*line != '(' ) {
            /* error, bad input format */
        }
        INPgetNetTok(&line, &nname1, 0);
        INPtermInsert(ckt, &nname1, tab, &node1);
        ptemp.nValue = node1;
        GCA(INPapName, (ckt, which, foo, "outpos", &ptemp));
        if (*line != ')') {
            INPgetNetTok(&line, &nname2, 1);
            INPtermInsert(ckt, &nname2, tab, &node2);
            ptemp.nValue = node2;
            GCA(INPapName, (ckt, which, foo, "outneg", &ptemp));
            ptemp.sValue = tprintf("V(%s,%s)", nname1, nname2);
            GCA(INPapName, (ckt, which, foo, "outname", &ptemp));
        } else {
            ptemp.nValue = gnode;
            GCA(INPapName, (ckt, which, foo, "outneg", &ptemp));
            ptemp.sValue = tprintf("V(%s)", nname1);
            GCA(INPapName, (ckt, which, foo, "outname", &ptemp));
        }
    } else if (*name == 'i' && strlen(name) == 1) {
        INPgetTok(&line, &name, 1);
        INPinsert(&name, tab);
        ptemp.uValue = name;
        GCA(INPapName, (ckt, which, foo, "outsrc", &ptemp));
    } else {
        LITERR("Syntax error: voltage or current expected.\n");
        return 0;
    }
    INPgetTok(&line, &name, 1);
    INPinsert(&name, tab);
    ptemp.uValue = name;
    GCA(INPapName, (ckt, which, foo, "insrc", &ptemp));
    return (0);
}


static int
dot_tran(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
         TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */
    double dtemp;		/* random double precision temporary */
    char *word;			/* something to stick a word of input into */

    NG_IGNORE(gnode);

    /* .tran Tstep Tstop <Tstart <Tmax> > <UIC> */
    which = ft_find_analysis("TRAN");
    if (which == -1) {
        LITERR("Transient analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Transient Analysis", &foo, task));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* Tstep */
    GCA(INPapName, (ckt, which, foo, "tstep", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* Tstop */
    GCA(INPapName, (ckt, which, foo, "tstop", parm));
    if (*line) {
        dtemp = INPevaluate(&line, &error, 1);	/* tstart? */
        if (error == 0) {
            ptemp.rValue = dtemp;
            GCA(INPapName, (ckt, which, foo, "tstart", &ptemp));
            dtemp = INPevaluate(&line, &error, 1);	/* tmax? */
            if (error == 0) {
                ptemp.rValue = dtemp;
                GCA(INPapName, (ckt, which, foo, "tmax", &ptemp));
            }
        }
    }
    if (*line) {
        INPgetTok(&line, &word, 1);	/* uic? */
        if (strcmp(word, "uic") == 0) {
            ptemp.iValue = 1;
            GCA(INPapName, (ckt, which, foo, "uic", &ptemp));
        } else {
            LITERR(" Error: unknown parameter on .tran - ignored\n");
        }
        tfree(word);
    }
    return (0);
}


static int
dot_sens(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
         TSKtask *task, CKTnode *gnode, JOB *foo)
{
    char *name;			/* the resistor's name */
    int error;			/* error code temporary */
    int filters, fidx;          /* Filter allocation and index. */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */
    char *nname1;		/* the first node's name */
    char *nname2;		/* the second node's name */
    CKTnode *node1;		/* the first node's node pointer */
    CKTnode *node2;		/* the second node's node pointer */
    char *steptype;		/* ac analysis, type of stepping function */
    char *cp;                   /* Scan for filters. */

    extern char **Sens_filter;  /* cktsens.c */

    which = ft_find_analysis("SENS");
    if (which == -1) {
        LITERR("Sensitivity unsupported.\n");
        return (0);
    }

    IFC(newAnalysis, (ckt, which, "Sensitivity Analysis", &foo, task));

    /* Format is:
     *      .sens <output> [<filter strings>]
     *      + [ac [dec|lin|oct] <pts> <low freq> <high freq> | dc ]
     */

    /* Get the output voltage or current */
    INPgetTok(&line, &name, 0);
    /* name is now either V or I or a serious error */
    if (*name == 'v' && strlen(name) == 1) {
        if (*line != '(') {
            LITERR("Syntax error: '(' expected after 'v'\n");
            return 0;
        }
        INPgetNetTok(&line, &nname1, 0);
        INPtermInsert(ckt, &nname1, tab, &node1);
        ptemp.nValue = node1;
        GCA(INPapName, (ckt, which, foo, "outpos", &ptemp));

        if (*line != ')') {
            INPgetNetTok(&line, &nname2, 1);
            INPtermInsert(ckt, &nname2, tab, &node2);
            ptemp.nValue = node2;
            GCA(INPapName, (ckt, which, foo, "outneg", &ptemp));
            ptemp.sValue = tprintf("V(%s,%s)", nname1, nname2);
            GCA(INPapName, (ckt, which, foo, "outname", &ptemp));
        } else {
            ptemp.nValue = gnode;
            GCA(INPapName, (ckt, which, foo, "outneg", &ptemp));
            ptemp.sValue = tprintf("V(%s)", nname1);
            GCA(INPapName, (ckt, which, foo, "outname", &ptemp));
        }
    } else if (*name == 'i' && strlen(name) == 1) {
        INPgetTok(&line, &name, 1);
        INPinsert(&name, tab);
        ptemp.uValue = name;
        GCA(INPapName, (ckt, which, foo, "outsrc", &ptemp));
    } else {
        LITERR("Syntax error: voltage or current expected.\n");
        return 0;
    }

    /* Scan for filters for the parameter names to be varied.
     * INPgetTok() breaks on '*' so scan by hand.
     */

    if (Sens_filter)
        FREE(Sens_filter);
    fidx = 0;
    filters = -1; // Ensure room for NULL.
    name = NULL;
    while (*line && *line > ' ')
            ++line;
    for (;;) {
        int l;

        while (*line && *line <= ' ')
            ++line;
        if (!*line)
            break;
        cp = line;
        while (*cp && *cp > ' ')
            ++cp;
        l = (int)(cp - line);
        name = TMALLOC(char, l + 1);
        strncpy(name, line, l);
        name[l] = 0;
        line = cp;
        if (!strcmp(name, "ac") || !strcmp(name, "dc"))
            break;
        if (fidx >= filters)
            Sens_filter = TREALLOC(char *, Sens_filter, filters + 8);
        filters += 8;
        Sens_filter[fidx++] = name;
        name = NULL;
    }
    if (Sens_filter) {
        Sens_filter[fidx] = NULL;
    }

    if (name && !strcmp(name, "ac")) {
        INPgetTok(&line, &steptype, 1);	/* get DEC, OCT, or LIN */
        ptemp.iValue = 1;
        GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
        parm = INPgetValue(ckt, &line, IF_INTEGER, tab); /* number of points */
        GCA(INPapName, (ckt, which, foo, "numsteps", parm));
        parm = INPgetValue(ckt, &line, IF_REAL, tab); /* fstart */
        GCA(INPapName, (ckt, which, foo, "start", parm));
        parm = INPgetValue(ckt, &line, IF_REAL, tab); /* fstop */
        GCA(INPapName, (ckt, which, foo, "stop", parm));
        return (0);
    } else if (name && *name && strcmp(name, "dc")) {
        /* Bad flag */
        LITERR("Syntax error: 'ac' or 'dc' expected.\n");
    }
    if (name)
        FREE(name);
    return 0;
}


#ifdef WANT_SENSE2
static int
dot_sens2(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
          TSKtask *task, CKTnode *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */
    char *token;		/* a token from the line */

    NG_IGNORE(gnode);

    /* .sens {AC} {DC} {TRAN} [dev=nnn parm=nnn]* */
    which = ft_find_analysis("SENS2");
    if (which == -1) {
        LITERR("Sensitivity-2 analysis unsupported\n");
        return (0);
    }

    IFC(newAnalysis, (ckt, which, "Sensitivity-2 Analysis", &foo, task));

    while (*line) {

        IFparm *if_parm;

        /* read the entire line */
        INPgetTok(&line, &token, 1);

        if_parm = ft_find_analysis_parm(which, token);

        if (!if_parm) {
            /* didn't find it! */
            LITERR(" Error: unknown parameter on .sens-ignored \n");
            continue;
        }

        /* found it, analysis which, parameter i */
        if (if_parm->dataType & IF_FLAG) {

            /* one of the keywords! */
            ptemp.iValue = 1;
            error = ft_sim->setAnalysisParm (ckt, foo,
                                             if_parm->id,
                                             &ptemp,
                                             NULL);
            if (error)
                current->error = INPerrCat(current->error, INPerror(error));

        } else {

            parm = INPgetValue(ckt, &line, if_parm->dataType, tab);
            error = ft_sim->setAnalysisParm (ckt, foo,
                                             if_parm->id,
                                             parm,
                                             NULL);
            if (error)
                current->error = INPerrCat(current->error, INPerror(error));
        }
    }

    return (0);
}
#endif

#ifdef WITH_PSS
/*SP: Steady State Analyis */
static int
dot_pss(char *line, void *ckt, INPtables *tab, struct card *current,
        void *task, void *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *nname;		/* the oscNode name */
    CKTnode *nnode;		/* the oscNode node */
    int which;			/* which analysis we are performing */
    char *word;			/* something to stick a word of input into */

    NG_IGNORE(gnode);

    /* .pss Fguess StabTime OscNode <UIC>*/
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic steady state analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Periodic Steady State Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "fguess", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "stabtime", parm));

    INPgetNetTok(&line, &nname, 0);
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));	/* OscNode given as string */

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    /* card #458: zero or more trailing tokens, `uic` and/or `method=shoot`
     * (order-independent), on top of the existing (required, positional)
     * arguments above -- the pre-existing single-optional-"uic" behaviour
     * is unchanged for decks that don't use `method=`.
     *
     * Trap (found running the first shoot-mode smoke test): INPgetTok()
     * treats '=' as a plain separator character (inpgtok.c), not part of
     * either token either side of it -- "method=shoot" and "method shoot"
     * both come back as two independent tokens "method" and "shoot", with
     * no '=' anywhere. A single-token strncmp(word,"method=",7) check
     * therefore never matches either half; read "method" as its own
     * keyword and pull the value as a second, separate token. */
    while (*line) {
        INPgetTok(&line, &word, 1);
        if (!*word)
            break;
        if (strcmp(word, "uic") == 0) {
            ptemp.iValue = 1;
            GCA(INPapName, (ckt, which, foo, "uic", &ptemp));
        } else if (strcmp(word, "method") == 0) {
            char *mword;
            INPgetTok(&line, &mword, 1);
            ptemp.sValue = mword;
            GCA(INPapName, (ckt, which, foo, "method", &ptemp));
        } else {
            fprintf(stderr,"Error: unknown parameter %s on .pss - ignored\n", word);
        }
    }
    return (0);
}
/* SP */

/* Enhancement-122: Periodic AC (PAC). Runs PSS then sweeps a small-signal input
 * frequency, solving the harmonic conversion matrix at each point. Reuses the PSS
 * analysis (the PAC sweep runs off the retained periodic operating point) with the
 * extra pac_* parameters set. */
static int
dot_pac(char *line, void *ckt, INPtables *tab, struct card *current,
        void *task, void *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *nname;		/* the oscNode name */
    CKTnode *nnode;		/* the oscNode node */
    int which;			/* which analysis we are performing */
    char *steptype;		/* pac sweep type: dec/oct/lin */

    NG_IGNORE(gnode);
    NG_IGNORE(current);

    /* .pac Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff
     *      <DEC|OCT|LIN> NumPts Fstart Fstop */
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic AC (PAC) analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Periodic AC Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "fguess", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "stabtime", parm));

    INPgetNetTok(&line, &nname, 0);
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));	/* OscNode given as string */

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    /* PAC sweep tail: <DEC|OCT|LIN> NumPts Fstart Fstop */
    INPgetTok(&line, &steptype, 1);
    ptemp.iValue = (strcmp(steptype, "dec") == 0) ? 1 :
                   (strcmp(steptype, "oct") == 0) ? 2 : 0;	/* default LIN */
    tfree(steptype);
    GCA(INPapName, (ckt, which, foo, "pac_step", &ptemp));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* number of points */
    GCA(INPapName, (ckt, which, foo, "pac_points", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstart */
    GCA(INPapName, (ckt, which, foo, "pac_fstart", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstop */
    GCA(INPapName, (ckt, which, foo, "pac_fstop", parm));

    /* card #465: optional trailing tokens, in any order --
     *   maxsideband          (bare integer, existing Enhancement-122 syntax)
     *   method shoot|driven  (which PSS algorithm the .pac card's OWN PSS
     *                         pass uses -- same keyword/value split as
     *                         dot_pss's tail, see that function's comment
     *                         on why INPgetTok() needs two separate
     *                         tokens for "method=shoot")
     *   solver td|fd         (card #465: time-domain vs the existing
     *                         frequency-domain PAC engine; td needs
     *                         method=shoot)
     * A bare integer is accepted only once, and only where the existing
     * Enhancement-122 maxsideband token could appear -- old decks with
     * just a trailing maxsideband are completely unaffected. */
    {
        int maxsb_seen = 0;
        while (*line) {
            char *word;
            char *save = line;
            INPgetTok(&line, &word, 1);
            if (!*word) { tfree(word); break; }
            if (!maxsb_seen && strcmp(word, "method") != 0 && strcmp(word, "solver") != 0 &&
                (isdigit((unsigned char) word[0]) ||
                 ((word[0] == '+' || word[0] == '-') && isdigit((unsigned char) word[1])))) {
                line = save;
                tfree(word);
                parm = INPgetValue(ckt, &line, IF_INTEGER, tab);
                GCA(INPapName, (ckt, which, foo, "pac_maxsb", parm));
                maxsb_seen = 1;
            } else if (strcmp(word, "method") == 0) {
                char *mword;
                tfree(word);
                INPgetTok(&line, &mword, 1);
                ptemp.sValue = mword;
                GCA(INPapName, (ckt, which, foo, "method", &ptemp));
            } else if (strcmp(word, "solver") == 0) {
                char *sword;
                tfree(word);
                INPgetTok(&line, &sword, 1);
                ptemp.sValue = sword;
                GCA(INPapName, (ckt, which, foo, "solver", &ptemp));
            } else {
                fprintf(stderr, "Error: unknown trailing token '%s' on .pac - ignored\n", word);
                tfree(word);
            }
        }
    }

    ptemp.iValue = 1;						/* enable the PAC sweep */
    GCA(INPapName, (ckt, which, foo, "pac", &ptemp));

    return (0);
}

/* card #465 (plan T5.1-T5.2): `.pstb` -- periodic loop gain. Runs its own PSS
 * (shooting-Newton only, method=shoot; the time-domain PAC/PSTB engine needs
 * the shooting orbit's per-step structure), then Tian two-injection loop-gain
 * sweep at a probe, on the sideband-0 LPTV response (pacshoot.c). Same PSS
 * params + sweep tail as .pac/.pss, plus a required probe=<Vname> (an
 * existing 0V V-source on the loop, same convention as native .stb's
 * probe=<Vname> loop-mode -- device-mode and probe=inst:terminal auto-insert
 * are NOT supported here, card #465 scope cut). */
static int
dot_pstb(char *line, void *ckt, INPtables *tab, struct card *current,
         void *task, void *gnode, JOB *foo)
{
    int error;
    IFvalue ptemp;
    IFvalue *parm;
    char *nname;
    CKTnode *nnode;
    int which;
    char *steptype;
    int probe_seen = 0;

    NG_IGNORE(gnode);
    NG_IGNORE(current);

    /* .pstb Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff
     *       <DEC|OCT|LIN> NumPts Fstart Fstop probe=<Vname> [method=shoot] */
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic loop-stability (PSTB) analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Periodic Loop Stability (PSTB) Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "fguess", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "stabtime", parm));

    INPgetNetTok(&line, &nname, 0);
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));	/* OscNode given as string */

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    /* sweep tail: <DEC|OCT|LIN> NumPts Fstart Fstop */
    INPgetTok(&line, &steptype, 1);
    ptemp.iValue = (strcmp(steptype, "dec") == 0) ? 1 :
                   (strcmp(steptype, "oct") == 0) ? 2 : 0;	/* default LIN */
    tfree(steptype);
    GCA(INPapName, (ckt, which, foo, "pac_step", &ptemp));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* number of points */
    GCA(INPapName, (ckt, which, foo, "pac_points", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstart */
    GCA(INPapName, (ckt, which, foo, "pac_fstart", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstop */
    GCA(INPapName, (ckt, which, foo, "pac_fstop", parm));

    /* trailing tokens: probe=<Vname> (required) and method shoot|driven
     * (optional, no '=' -- same INPgetTok split as dot_pss/dot_pac). */
    while (*line) {
        char *word;
        while (*line == ' ' || *line == '\t' || *line == '+')
            line++;
        if (!*line)
            break;
        if (strncmp(line, "probe=", 6) == 0) {
            char *name;
            line += 6;
            INPgetNetTok(&line, &name, 1);
            if (!*name) {
                LITERR("'.pstb' probe=<Vname> needs a non-empty probe name\n");
                tfree(name);
                return (0);
            }
            INPinsert(&name, tab);
            ptemp.uValue = name;
            GCA(INPapName, (ckt, which, foo, "pstb_probe", &ptemp));
            probe_seen = 1;
            continue;
        }
        INPgetTok(&line, &word, 1);
        if (!*word) { tfree(word); break; }
        if (strcmp(word, "method") == 0) {
            char *mword;
            tfree(word);
            INPgetTok(&line, &mword, 1);
            ptemp.sValue = mword;
            GCA(INPapName, (ckt, which, foo, "method", &ptemp));
        } else {
            fprintf(stderr, "Error: unknown trailing token '%s' on .pstb - ignored\n", word);
            tfree(word);
        }
    }

    if (!probe_seen) {
        LITERR("'.pstb' needs probe=<Vname> (an existing 0V V-source on the loop)\n");
        return (0);
    }

    ptemp.iValue = 1;						/* enable the PSTB sweep */
    GCA(INPapName, (ckt, which, foo, "pstb", &ptemp));

    return (0);
}

/* Enhancement-132: Periodic S-parameters (PSP). Runs PSS then, for each RF port,
 * injects through the harmonic conversion matrix and forms the periodic scattering
 * matrix vs input frequency. Same PSS params + sweep tail as .pac; the excitation is
 * the netlist's RF ports (portnum/z0), so there is no output node. */
static int
dot_psp(char *line, void *ckt, INPtables *tab, struct card *current,
        void *task, void *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *nname;		/* the oscNode name */
    CKTnode *nnode;		/* the oscNode node */
    int which;			/* which analysis we are performing */
    char *steptype;		/* psp sweep type: dec/oct/lin */

    NG_IGNORE(gnode);
    NG_IGNORE(current);

    /* .psp Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff
     *      <DEC|OCT|LIN> NumPts Fstart Fstop [maxsideband] */
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic S-parameter (PSP) analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Periodic S-parameter Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "fguess", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "stabtime", parm));

    INPgetNetTok(&line, &nname, 0);
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));	/* OscNode given as string */

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    /* PSP sweep tail: <DEC|OCT|LIN> NumPts Fstart Fstop */
    INPgetTok(&line, &steptype, 1);
    ptemp.iValue = (strcmp(steptype, "dec") == 0) ? 1 :
                   (strcmp(steptype, "oct") == 0) ? 2 : 0;	/* default LIN */
    tfree(steptype);
    GCA(INPapName, (ckt, which, foo, "pac_step", &ptemp));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* number of points */
    GCA(INPapName, (ckt, which, foo, "pac_points", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstart */
    GCA(INPapName, (ckt, which, foo, "pac_fstart", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstop */
    GCA(INPapName, (ckt, which, foo, "pac_fstop", parm));

    {   /* optional trailing maxsideband: output conversion sidebands each side */
        char *p = line;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p) {
            parm = INPgetValue(ckt, &line, IF_INTEGER, tab);
            GCA(INPapName, (ckt, which, foo, "pac_maxsb", parm));
        }
    }

    ptemp.iValue = 1;						/* enable the PSP sweep */
    GCA(INPapName, (ckt, which, foo, "psp", &ptemp));

    return (0);
}

/* Enhancement-124: Periodic noise (PNOISE). Runs PSS then folds each device's
 * noise through the conversion-matrix adjoint over all sidebands. Reuses the PSS
 * analysis (like .pac) with the pnoise output node, input source, and sweep set. */
static int
dot_pnoise(char *line, void *ckt, INPtables *tab, struct card *current,
           void *task, void *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *nname;		/* a node name */
    char *sname;		/* the input source name */
    CKTnode *nnode;		/* a node pointer */
    int which;			/* which analysis we are performing */
    char *steptype;		/* pnoise sweep type: dec/oct/lin */

    NG_IGNORE(gnode);
    NG_IGNORE(current);

    /* .pnoise Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff
     *         OutNode InSrc <DEC|OCT|LIN> NumPts Fstart Fstop */
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic noise (PNOISE) analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Periodic Noise Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "fguess", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "stabtime", parm));
    INPgetNetTok(&line, &nname, 0);				/* OscNode */
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    INPgetNetTok(&line, &nname, 0);				/* OutNode */
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "pnoise_out", &ptemp));

    INPgetTok(&line, &sname, 1);				/* InSrc */
    INPinsert(&sname, tab);
    ptemp.uValue = sname;
    GCA(INPapName, (ckt, which, foo, "pnoise_insrc", &ptemp));

    /* sweep tail: <DEC|OCT|LIN> NumPts Fstart Fstop */
    INPgetTok(&line, &steptype, 1);
    ptemp.iValue = (strcmp(steptype, "dec") == 0) ? 1 :
                   (strcmp(steptype, "oct") == 0) ? 2 : 0;
    tfree(steptype);
    GCA(INPapName, (ckt, which, foo, "pac_step", &ptemp));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* number of points */
    GCA(INPapName, (ckt, which, foo, "pac_points", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstart */
    GCA(INPapName, (ckt, which, foo, "pac_fstart", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstop */
    GCA(INPapName, (ckt, which, foo, "pac_fstop", parm));

    {   /* Enhancement-126: optional trailing "cyclo" keyword */
        char *p = line;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p) {
            char *word;
            INPgetTok(&line, &word, 1);
            if (strcmp(word, "cyclo") == 0) {
                ptemp.iValue = 1;
                GCA(INPapName, (ckt, which, foo, "pnoise_cyclo", &ptemp));
            } else {
                fprintf(stderr, "Error: unknown parameter %s on .pnoise - ignored\n", word);
            }
            tfree(word);
        }
    }

    ptemp.iValue = 1;						/* enable the pnoise sweep */
    GCA(INPapName, (ckt, which, foo, "pnoise", &ptemp));

    return (0);
}

/* Enhancement-125: Periodic transfer function (PXF). The adjoint of PAC: runs PSS
 * then solves Hᵀ Ψ = e_{out,0} and dots Ψ with the netlist AC-source pattern to get
 * the transfer from the input to a fixed output at each sideband. Reuses the PSS
 * analysis (like .pac). */
static int
dot_pxf(char *line, void *ckt, INPtables *tab, struct card *current,
        void *task, void *gnode, JOB *foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue *parm;		/* a pointer to a value struct for function returns */
    char *nname;		/* a node name */
    CKTnode *nnode;		/* a node pointer */
    int which;			/* which analysis we are performing */
    char *steptype;		/* pxf sweep type: dec/oct/lin */

    NG_IGNORE(gnode);
    NG_IGNORE(current);

    /* .pxf Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff
     *      OutNode <DEC|OCT|LIN> NumPts Fstart Fstop [maxsideband] */
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic transfer-function (PXF) analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Periodic Transfer Function Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "fguess", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "stabtime", parm));
    INPgetNetTok(&line, &nname, 0);				/* OscNode */
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    INPgetNetTok(&line, &nname, 0);				/* OutNode */
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "pxf_out", &ptemp));

    /* sweep tail: <DEC|OCT|LIN> NumPts Fstart Fstop [maxsideband] */
    INPgetTok(&line, &steptype, 1);
    ptemp.iValue = (strcmp(steptype, "dec") == 0) ? 1 :
                   (strcmp(steptype, "oct") == 0) ? 2 : 0;
    tfree(steptype);
    GCA(INPapName, (ckt, which, foo, "pac_step", &ptemp));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* number of points */
    GCA(INPapName, (ckt, which, foo, "pac_points", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstart */
    GCA(INPapName, (ckt, which, foo, "pac_fstart", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* fstop */
    GCA(INPapName, (ckt, which, foo, "pac_fstop", parm));

    {   /* optional trailing maxsideband */
        char *p = line;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p) {
            parm = INPgetValue(ckt, &line, IF_INTEGER, tab);
            GCA(INPapName, (ckt, which, foo, "pac_maxsb", parm));
        }
    }

    ptemp.iValue = 1;						/* enable the pxf sweep */
    GCA(INPapName, (ckt, which, foo, "pxf", &ptemp));

    return (0);
}
#endif


#ifdef RFSPICE
/* S-Parameter Analyis */
static int
dot_sp(char* line, void* ckt, INPtables* tab, struct card* current,
    void* task, void* gnode, JOB* foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue* parm;		/* a pointer to a value struct for function returns */
    int which;			/* which analysis we are performing */
    char* steptype;		/* ac analysis, type of stepping function */

    NG_IGNORE(gnode);

    /* .ac {DEC OCT LIN} NP FSTART FSTOP */
    which = ft_find_analysis("SP");
    if (which == -1) {
        LITERR("S-Params analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "SP Analysis", &foo, task));
    INPgetTok(&line, &steptype, 1);	/* get DEC, OCT, or LIN */
    ptemp.iValue = 1;
    GCA(INPapName, (ckt, which, foo, steptype, &ptemp));
    tfree(steptype);
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab); /* number of points */
    GCA(INPapName, (ckt, which, foo, "numsteps", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* fstart */
    GCA(INPapName, (ckt, which, foo, "start", parm));
    parm = INPgetValue(ckt, &line, IF_REAL, tab);	/* fstop */
    GCA(INPapName, (ckt, which, foo, "stop", parm));
    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);	/* fstop */
    GCA(INPapName, (ckt, which, foo, "donoise", parm));
    return (0);
}

#ifdef WITH_HB
/*SP: Steady State Analyis */
static int
dot_hb(char* line, void* ckt, INPtables* tab, struct card* current,
    void* task, void* gnode, JOB* foo)
{
    int error;			/* error code temporary */
    IFvalue ptemp;		/* a value structure to package resistance into */
    IFvalue* parm;		/* a pointer to a value struct for function returns */
    char* nname;		/* the oscNode name */
    CKTnode* nnode;		/* the oscNode node */
    int which;			/* which analysis we are performing */
    char* word;			/* something to stick a word of input into */

    NG_IGNORE(gnode);

    /* .pss Fguess StabTime OscNode <UIC>*/
    which = ft_find_analysis("PSS");
    if (which == -1) {
        LITERR("Periodic steady state analysis unsupported.\n");
        return (0);
    }
    IFC(newAnalysis, (ckt, which, "Harmonic Balance State Analysis", &foo, task));

    parm = INPgetValue(ckt, &line, IF_REALVEC, tab);		/* Fguess */
    GCA(INPapName, (ckt, which, foo, "freq", parm));

    parm = INPgetValue(ckt, &line, IF_INTVEC, tab);		/* StabTime */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));

    INPgetNetTok(&line, &nname, 0);
    INPtermInsert(ckt, &nname, tab, &nnode);
    ptemp.nValue = nnode;
    GCA(INPapName, (ckt, which, foo, "oscnode", &ptemp));	/* OscNode given as string */

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS points */
    GCA(INPapName, (ckt, which, foo, "points", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* PSS harmonics */
    GCA(INPapName, (ckt, which, foo, "harmonics", parm));

    parm = INPgetValue(ckt, &line, IF_INTEGER, tab);		/* SC iterations */
    GCA(INPapName, (ckt, which, foo, "sc_iter", parm));

    parm = INPgetValue(ckt, &line, IF_REAL, tab);		/* Steady coefficient */
    GCA(INPapName, (ckt, which, foo, "steady_coeff", parm));

    if (*line) {
        INPgetTok(&line, &word, 1);	/* uic? */
        if (strcmp(word, "uic") == 0) {
            ptemp.iValue = 1;
            GCA(INPapName, (ckt, which, foo, "uic", &ptemp));
        }
        else {
            fprintf(stderr, "Error: unknown parameter %s on .pss - ignored\n", word);
        }
    }
    return (0);
}
#endif

#endif


static int
dot_options(char *line, CKTcircuit *ckt, INPtables *tab, struct card *current,
            TSKtask *task, CKTnode *gnode, JOB *foo)
{
    NG_IGNORE(line);
    NG_IGNORE(gnode);
    NG_IGNORE(foo);

    /* .option - specify program options - rather complicated */
    /* use a subroutine to handle all of them to keep this    */
    /* subroutine managable.                                  */

    INPdoOpts(ckt, &(task->taskOptions), current, tab);
    return (0);
}


int
INP2dot(CKTcircuit *ckt, INPtables *tab, struct card *current, TSKtask *task, CKTnode *gnode)
{

    /* .<something> Many possibilities */
    char *token;		/* a token from the line, tmalloc'ed */
    JOB *foo = NULL;		/* pointer to analysis */
    /* the part of the current line left to parse */
    char *line = current->line;
    int rtn = 0;

    INPgetTok(&line, &token, 1);
    if (strcmp(token, ".model") == 0) {
        /* don't have to do anything, since models were all done in
         * pass 1 */
        goto quit;
    } else if (strcmp(token, ".param") == 0) {
        /* don't have to do anything, since params were all done
         * elsewhere */
        goto quit;
    } else if ((strcmp(token, ".width") == 0) ||
               strcmp(token, ".print") == 0 || strcmp(token, ".plot") == 0) {
        /* obsolete - ignore */
        char* token2 = tprintf(" obsolete dot command '%s' - ignored \n", token);
        LITERR(token2);
        tfree(token2);
        goto quit;
    } else if ((strcmp(token, ".temp") == 0)) {
        /* .temp temp1 temp2 temp3 temp4 ..... */
        /* not yet implemented - warn & ignore */
        /*
        LITERR(" Warning: .TEMP card obsolete - use .options TEMP and TNOM\n");
        */
        goto quit;
    } else if ((strcmp(token, ".op") == 0)) {
        rtn = dot_op(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".nodeset") == 0)) {
        goto quit;
    } else if ((strcmp(token, ".disto") == 0)) {
        rtn = dot_disto(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".noise") == 0)) {
        rtn = dot_noise(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".four") == 0)
               || (strcmp(token, ".fourier") == 0)) {
        /* .four */
        /* not implemented - warn & ignore */
        LITERR("Use fourier command to obtain fourier analysis\n");
        goto quit;
    } else if ((strcmp(token, ".ic") == 0)) {
        goto quit;
    } else if ((strcmp(token, ".ac") == 0)) {
        rtn = dot_ac(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".stb") == 0)) {
        rtn = dot_stb(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".pz") == 0)) {
        rtn = dot_pz(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".dc") == 0)) {
        rtn = dot_dc(line, ckt, tab, current, task, gnode, foo);
        if (rtn == 1) {
            current->error = copy("Bad syntax! ");
        }
        goto quit;
    } else if ((strcmp(token, ".tf") == 0)) {
        rtn = dot_tf(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    } else if ((strcmp(token, ".tran") == 0)) {
        rtn = dot_tran(line, ckt, tab, current, task, gnode, foo);
        goto quit;
#ifdef WITH_PSS
        /* SP: Steady State Analysis */
    } else if ((strcmp(token, ".pss") == 0)) {
        rtn = dot_pss(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* SP */
        /* Enhancement-122: Periodic AC */
    } else if ((strcmp(token, ".pac") == 0)) {
        rtn = dot_pac(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* card #465: Periodic loop gain (PSTB) */
    } else if ((strcmp(token, ".pstb") == 0)) {
        rtn = dot_pstb(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* Enhancement-124: Periodic noise */
    } else if ((strcmp(token, ".pnoise") == 0)) {
        rtn = dot_pnoise(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* Enhancement-125: Periodic transfer function */
    } else if ((strcmp(token, ".pxf") == 0)) {
        rtn = dot_pxf(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* Enhancement-132: Periodic S-parameters */
    } else if ((strcmp(token, ".psp") == 0)) {
        rtn = dot_psp(line, ckt, tab, current, task, gnode, foo);
        goto quit;
#endif
#ifdef RFSPICE
    }
    else if ((strcmp(token, ".sp") == 0)) {
        rtn = dot_sp(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* SP */
#ifdef WITH_HB
    }
    else if ((strcmp(token, ".hb") == 0)) {
        rtn = dot_hb(line, ckt, tab, current, task, gnode, foo);
        goto quit;
        /* SP */
#endif
#endif
    } else if ((strcmp(token, ".subckt") == 0) ||
               (strcmp(token, ".ends") == 0)) {
        /* not yet implemented - warn & ignore */
        LITERR(" Warning: Subcircuits not yet implemented - ignored \n");
        goto quit;
    } else if ((strcmp(token, ".end") == 0)) {
        /* .end - end of input */
        /* not allowed to pay attention to additional input - return */
        rtn = 1;
        goto quit;
    } else if (strcmp(token, ".sens") == 0) {
        rtn = dot_sens(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    }
#ifdef WANT_SENSE2
    else if ((strcmp(token, ".sens2") == 0)) {
        rtn = dot_sens2(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    }
#endif
    else if ((strcmp(token, ".probe") == 0)) {
        /* Maybe generate a "probe" format file in the future. */
        goto quit;
    } else if ((strcmp(token, ".options") == 0)||
               (strcmp(token,".option")==0) ||
               (strcmp(token,".opt")==0)) {
        rtn = dot_options(line, ckt, tab, current, task, gnode, foo);
        goto quit;
    }
    /* Added by H.Tanaka to find .global option */
    else if (strcmp(token, ".global") == 0) {
        rtn = 0;
        LITERR(" Warning: .global not yet implemented - ignored \n");
        goto quit;
    }
    /* ignore .meas statements -- these will be handled after analysis */
    /* also ignore .param statements */
    /* ignore .prot, .unprot */
    else if (strcmp(token, ".meas") == 0 || ciprefix(".para", token) || strcmp(token, ".measure") == 0 ||
             strcmp(token, ".prot") == 0 || strcmp(token, ".unprot") == 0) {
        rtn = 0;
        goto quit;
    }
    char *token2 = tprintf(" unimplemented dot command '%s'\n", token);
    LITERR(token2);
    tfree(token2);
quit:
    tfree(token);
    return rtn;
}
