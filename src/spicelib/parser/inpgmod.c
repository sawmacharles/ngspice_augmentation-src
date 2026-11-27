/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Copyright 2000 The ngspice team
3-Clause BSD license
(see COPYING or https://opensource.org/licenses/BSD-3-Clause)
Author: 1985 Thomas L. Quarles, 1991 David A. Gates
Modified: 2001 Paolo Nenzi (Cider Integration)
**********/

#include "ngspice/ngspice.h"
#include "ngspice/inpdefs.h"
#include "ngspice/ifsim.h"
#include "ngspice/cpstd.h"
#include "ngspice/fteext.h"
#include "ngspice/compatmode.h"
#include "ngspice/devdefs.h"
#include "inpxx.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>

#ifdef CIDER

#include "ngspice/numcards.h"
#include "ngspice/carddefs.h"
#include "ngspice/numgen.h"
#include "ngspice/suffix.h"

#define E_MISSING    -1
#define E_AMBIGUOUS  -2

extern IFcardInfo *INPcardTab[];
extern int INPnumCards;

static int INPparseNumMod(CKTcircuit *ckt, INPmodel *model, INPtables *tab, char **errMessage);
static int INPfindCard(char *name, IFcardInfo *table[], int numCards);
static int INPfindParm(char *name, IFparm *table, int numParms);

#endif

extern INPmodel *modtab;
extern NGHASHPTR modtabhash;


static IFparm *
find_model_parameter(const char *name, IFdevice *device)
{
    IFparm *p = device->modelParms;
    IFparm *p_end = p + *(device->numModelParms);

    for (; p < p_end; p++)
        if (strcmp(name, p->keyword) == 0)
            return p;

    return NULL;
}


static IFparm *
find_instance_parameter(const char *name, IFdevice *device)
{
    IFparm *p = device->instanceParms;
    IFparm *p_end = p + *(device->numInstanceParms);

    for (; p < p_end; p++)
        if (strcmp(name, p->keyword) == 0)
            return p;

    return NULL;
}


/*
 * code moved from INPgetMod
 */
static int
create_model(CKTcircuit *ckt, INPmodel *modtmp, INPtables *tab)
{
    char    *err = NULL, *line, *parm, *endptr;
    int     error;

    /* not already defined, so create & give parameters */
    error = ft_sim->newModel(ckt, modtmp->INPmodType, &(modtmp->INPmodfast), modtmp->INPmodName);
    if (error)
        return error;

#ifdef CIDER
    /* Handle Numerical Models Differently */
    if (modtmp->INPmodType == INPtypelook("NUMD") ||
        modtmp->INPmodType == INPtypelook("NBJT") ||
        modtmp->INPmodType == INPtypelook("NUMD2") ||
        modtmp->INPmodType == INPtypelook("NBJT2") ||
        modtmp->INPmodType == INPtypelook("NUMOS"))
    {
        error = INPparseNumMod(ckt, modtmp, tab, &err);
        if (error)
            return error;
        modtmp->INPmodLine->error = err;
        return 0;
    }
#endif

    IFdevice *device = ft_sim->devices[modtmp->INPmodType];

    /* parameter isolation, identification, binding */

    line = modtmp->INPmodLine->line;

#ifdef TRACE
    printf("In INPgetMod, inserting new model into table.  line = %s ...\n", line);
#endif

    INPgetTok(&line, &parm, 1);        /* throw away '.model' */
    tfree(parm);
    INPgetNetTok(&line, &parm, 1);        /* throw away 'modname' */
    tfree(parm);

#ifdef OSDI
    /* osdi models don't accept their device type as an argument */
    bool is_osdi = false;
    /* card of a model level bound to this module with osdi_bind */
    const struct osdi_bind *bind = NULL;
    const char *polarity = NULL;        /* value given to bind->polarity */
    if (device->registry_entry){
        /* osdi module name, or nmos/pmos/npn/pnp/d of a card bound with osdi_bind */
        INPgetNetTok(&line, &parm, 1);
        bool is_npn = strcmp(parm, "npn") == 0, is_pnp = strcmp(parm, "pnp") == 0;
        bool is_n = strcmp(parm, "nmos") == 0 || strcmp(parm, "nsoi") == 0 || is_npn;
        bool is_p = strcmp(parm, "pmos") == 0 || strcmp(parm, "psoi") == 0 || is_pnp;
        const char *kind = (is_npn || is_pnp) ? "bjt" :
                           (is_n || is_p) ? "mos" :
                           strcmp(parm, "d") == 0 ? "dio" : NULL;
        if (kind) {
            bind = INPosdiBindForCard(kind, modtmp->INPmodLine->line);
            if (bind && strcasecmp(bind->module, device->name) != 0)
                bind = NULL;
        }
        if (bind && bind->polarity && (is_n || is_p)) {
            IFparm *p = find_model_parameter(bind->polarity, device);
            if (p) {
                char *value = is_n ? bind->nval : bind->pval;
                polarity = value;
                error = ft_sim->setModelParm(ckt, modtmp->INPmodfast, p->id,
                                             INPgetValue(ckt, &value, p->dataType, tab), NULL);
                if (error)
                    return error;
            } else {
                err = INPerrCat(err, tprintf("module %s has no parameter %s - %s ignored"
                                             " (osdi_bind ... polarity=<par>|none)",
                                             device->name, bind->polarity, parm));
            }
        }
        tfree(parm);
        is_osdi = true;
    }
#endif

    while (*line) {
        INPgetTok(&line, &parm, 1);
        if (!*parm) {
            FREE(parm);
            continue;
        }

#ifdef OSDI
        if (bind) {
            /* level selects the binding, it is not a parameter of the module;
               nor is version, unless the module has this parameter */
            if (strcmp(parm, "level") == 0 ||
                (strcmp(parm, "version") == 0 && !find_model_parameter(parm, device))) {
                char *value;
                INPgetTok(&line, &value, 1);    /* a number or a string like "1.1" */
                tfree(value);
                FREE(parm);
                continue;
            }
            /* the polarity of a bound card comes from nmos/pmos */
            if (polarity && strcmp(parm, bind->polarity) == 0) {
                double v = INPgetValue(ckt, &line, IF_REAL, tab)->rValue;
                if (v != atof(polarity))
                    err = INPerrCat(err, tprintf("%s=%g contradicts the model type - ignored",
                                                 parm, v));
                FREE(parm);
                continue;
            }
            /* parameters of other simulators listed with 'osdi_bind ... ignore' */
            if (wl_find(parm, bind->ignore)) {
                INPgetValue(ckt, &line, IF_REAL, tab);
                FREE(parm);
                continue;
            }
        }
#endif

        IFparm *p = find_model_parameter(parm, device);

        if (p) {
#ifdef OSDI
            if (is_osdi && (p->dataType & IF_VECTOR)){  
                // we need to get rid if the leading [ in order to make sure 
                // that INPgetValue can parse the value properly
                // This is because, unlike other SPICEDev, OSDI models receive
                // array params in the syntax (param_name=[...])
                ++line;
            }
#endif
            IFvalue *val = INPgetValue(ckt, &line, p->dataType, tab);
            error = ft_sim->setModelParm(ckt, modtmp->INPmodfast, p->id, val, NULL);
            if (error)
                return error;
        } else if ((strcmp(parm, "level") == 0) || (strcmp(parm, "m") == 0)) {
            /* no instance parameter default for level and multiplier */
            /* just grab the number and throw away */
            /* since we already have that info from pass1 */
            INPgetValue(ckt, &line, IF_REAL, tab);
        
        } else {

            p = find_instance_parameter(parm, device);

            if (p) {
                char *value;

                INPgetTok(&line, &value, 1);
                if (p->dataType & IF_SET) {
                    modtmp->INPmodfast->defaults =
                        wl_cons(copy(parm),
                                wl_cons(value, modtmp->INPmodfast->defaults));
                } else {
                    fprintf(stderr,
                            "Ignoring attempt to set a default "
                            "for read-only instance parameter %s in:\n  %s\n",
                            p->keyword, modtmp->INPmodLine->line);
                }
            } else {

                double dval;

                /* want only the parameter names in output - not the values */
                errno = 0;    /* To distinguish success/failure after call */
                dval = strtod(parm, &endptr);
                /* Check for various possible errors */
                if ((errno == ERANGE && dval == HUGE_VAL) || errno != 0) {
                    perror("strtod");
                    controlled_exit(EXIT_FAILURE);
                }
                if (endptr == parm) /* it was no number - it is really a string */
                    err = INPerrCat(err,
                                    tprintf("unrecognized parameter (%s) - ignored",
                                            parm));
            }
        }
        FREE(parm);
    }

    modtmp->INPmodLine->error = err;
    return 0;
}


/* The next token of line as INPgetNetTok(&line, &token, 1) finds it, without copying it: *start and *len give
   the token, line moves past it and the separators after it. Card #2060: parse_line ran INPgetNetTok (one
   allocation and free per token) over every candidate card's whole parameter list for each instance. */
static void
next_net_tok(char **line, char **start, size_t *len)
{
    char *p = *line, *q;
    while (*p == ' ' || *p == '\t' || *p == '=' || *p == '(' || *p == ')' || *p == ',')
        p++;
    for (q = p; *q && *q != ' ' && *q != '\t' && *q != '\r' && *q != '=' && *q != ',' && *q != '(' && *q != ')'; q++)
        ;
    if (q == p && *q)
        q++;
    *start = p;
    *len = (size_t) (q - p);
    while (*q == ' ' || *q == '\t' || *q == '\r' || *q == '=' || *q == ',')
        q++;
    *line = q;
}


static int
binsel_nocache(void)
{
    static int v = -1;
    if (v < 0)
        v = getenv("BINSEL_NOCACHE") ? 1 : 0;
    return v;
}


static bool
parse_line(char *line, char *tokens[], int num_tokens, double values[], bool found[])
{
    int get_index = -1;
    int i;
    bool old = binsel_nocache();

    for (i = 0; i < num_tokens; i++)
        found[i] = FALSE;

    while (*line) {

        if (get_index != -1) {
            int error;
            values[get_index] = INPevaluate(&line, &error, 1);
            found[get_index] = TRUE;
            get_index = -1;
            continue;
        }

        if (old) {
            char *token = NULL;
            INPgetNetTok(&line, &token, 1);
            for (i = 0; i < num_tokens; i++)
                if (strcmp(tokens[i], token) == 0)
                    get_index = i;
            txfree(token);
        } else {
            char *t;
            size_t n;
            next_net_tok(&line, &t, &n);
            for (i = 0; i < num_tokens; i++)
                if (strncmp(tokens[i], t, n) == 0 && tokens[i][n] == '\0')
                    get_index = i;
        }
    }

    for (i = 0; i < num_tokens; i++)
        if (!found[i])
            return FALSE;

    return TRUE;
}


static bool
is_equal(double result, double expectedResult)
{
    return fabs(result - expectedResult) <= 1e-9 * fabs(expectedResult);
}


bool
INPbinInRange(double value, double min, double max)
{
    /* the standard binning rule is: min <= value < max */
    return (is_equal(value, min) || min < value) && value < max && !is_equal(value, max);
}


/* wnflag for the selection of a bin: 1 for w/nf, 0 for the total w.
   From set wnflag=0|1 or .option wnflag=0|1; set wnflag or .option wnflag
   without a value is 1, as in reference simulators.  If not given, 1 with ngbehavior hs
   or spe, else 0.  wnflag on a BSIM4 instance overrides this. */
int
INPgetWnflag(void)
{
    int wnflag;

    if (cp_getvar("wnflag", CP_NUM, &wnflag, 0))
        return wnflag ? 1 : 0;
    if (cp_getvar("wnflag", CP_BOOL, NULL, 0))
        return 1;
    return (newcompat.spe || newcompat.hs) ? 1 : 0;
}


/* Is the variable binoutermax set (set binoutermax, or .option binoutermax)?
   Then a size exactly on lmax or wmax is accepted when no bin holds it
   by the standard rule. Some PDKs (e.g. sky130) have outermost bin limits
   on round values, so their largest device size is on lmax or wmax. */
bool
INPbinOuterMax(void)
{
    return cp_getvar("binoutermax", CP_BOOL, NULL, 0);
}


/* Match the sizes l and w against the bin limits of one model card.
   0: the card does not hold the sizes.
   1: min <= size < max for l and for w (the standard rule).
   2, 3: only with outermax: one or both sizes are exactly on lmax, wmax.
   A lower value is the better match. */
int
INPbinMatch(double l, double lmin, double lmax,
            double w, double wmin, double wmax, bool outermax)
{
    int match = 1;

    if (!INPbinInRange(l, lmin, lmax)) {
        if (outermax && is_equal(l, lmax))
            match++;
        else
            return 0;
    }
    if (!INPbinInRange(w, wmin, wmax)) {
        if (outermax && is_equal(w, wmax))
            match++;
        else
            return 0;
    }
    return match;
}


/* The bin rule of ngspice before upstream commit 03ef19244 (2026-10-30, "Correct the limits for binning
   of MOS model"): a size within an absolute 1e-9 of lmin, lmax, wmin or wmax, or strictly inside, is in the
   bin, and the first card found wins. At an exact bin edge two cards hold the size and the card order
   decides; the standard rule (min <= size < max) takes the card that starts at the edge. */
static bool
legacy_in_range(double value, double min, double max)
{
    return fabs(value - min) < 1e-9 || fabs(value - max) < 1e-9 || (min < value && value < max);
}


int
INPbinMatchLegacy(double l, double lmin, double lmax, double w, double wmin, double wmax)
{
    return (legacy_in_range(l, lmin, lmax) && legacy_in_range(w, wmin, wmax)) ? 1 : 0;
}


/* Is binlegacy set (set binlegacy, or .option binlegacy)? Then INPgetModBin() bins by the old rule. */
bool
INPbinLegacy(void)
{
    return cp_getvar("binlegacy", CP_BOOL, NULL, 0);
}


static bool
binnable_type(int type)
{
    if (type == INPtypelook("BSIM3") || type == INPtypelook("BSIM3v32") || type == INPtypelook("BSIM3v0") ||
        type == INPtypelook("BSIM3v1") || type == INPtypelook("BSIM4") || type == INPtypelook("BSIM4v5") ||
        type == INPtypelook("BSIM4v6") || type == INPtypelook("BSIM4v7") || type == INPtypelook("HiSIM2") ||
        type == INPtypelook("HiSIMHV1") || type == INPtypelook("HiSIMHV2"))
        return TRUE;
#ifdef OSDI
    /* Verilog-A modules: bins by lmin, lmax, wmin, wmax of the card */
    if (type >= 0 && ft_sim->devices[type]->registry_entry)
        return TRUE;
#endif
    return FALSE;
}


/* #1191: index of the binned cards by base name (the name without its trailing ".<digits>"), in modtab order.
   Subcircuits that carry their own binned .model cards (sky130 and other PDKs) give every instance its own copy
   of the cards, so modtab grows with the circuit, and scanning all of it for every instance was quadratic: a
   5120-transistor sky130 netlist spent 12 minutes here. The index is rebuilt when modtab changes (models are
   prepended, so its head pointer changes) and lists exactly the cards model_name_match(name, card) == 2 selects,
   in the same order. BINSEL_SCAN=1 keeps the old full scan (red proof). */
/* lim[k] caches parse_line's lmin, lmax, wmin, wmax of v[k] (card #2060): state[k] 0 not read yet, 1 read,
   2 a size limit is missing; line[k] is the card text it was read from (read again if that changes). */
typedef double bin_lim[4];
typedef struct { INPmodel **v; bin_lim *lim; char *state; char **line; int n, cap; } bin_list;
static NGHASHPTR  bin_index = NULL;
static INPmodel  *bin_index_head = NULL;
static NGHASHPTR  bin_index_modhash = NULL;

static void
bin_list_free(void *p)
{
    bin_list *b = (bin_list *) p;
    tfree(b->v);
    tfree(b->lim);
    tfree(b->state);
    tfree(b->line);
    tfree(b);
}

static void
bin_key_free(void *p)
{
    tfree(p);
}

static void
bin_index_build(void)
{
    INPmodel *m;
    if (bin_index)
        nghash_free(bin_index, bin_list_free, bin_key_free);
    bin_index = nghash_init(NGHASH_MIN_SIZE);
    nghash_unique(bin_index, TRUE);
    for (m = modtab; m; m = m->INPnextModel) {
        const char *name = m->INPmodName, *dot = strrchr(name, '.'), *p;
        bin_list *b;
        char *key;
        if (!dot || dot == name || dot[1] == '\0')
            continue;
        for (p = dot + 1; *p; p++)
            if (!isdigit_c(*p))
                break;
        if (*p)
            continue;
        key = copy_substring(name, dot);
        b = (bin_list *) nghash_find(bin_index, key);
        if (!b) {
            b = TMALLOC(bin_list, 1);   /* zeroed */
            nghash_insert(bin_index, key, b);   /* the table keeps key */
        } else {
            tfree(key);
        }
        if (b->n == b->cap) {
            b->cap = b->cap ? 2 * b->cap : 8;
            b->v = TREALLOC(INPmodel *, b->v, b->cap);
            b->lim = TREALLOC(bin_lim, b->lim, b->cap);
            b->state = TREALLOC(char, b->state, b->cap);
            b->line = TREALLOC(char *, b->line, b->cap);
        }
        b->state[b->n] = 0;
        b->line[b->n] = NULL;
        b->v[b->n++] = m;
    }
    bin_index_head = modtab;
    bin_index_modhash = modtabhash;
}


/* The card for sizes l and w among the binned cards of model name: by the standard rule (with outermax as
   given), or by the legacy rule. *err is set for a card of unknown device type. */
static INPmodel *
bin_select(char *name, double l, double w, double w_nf, bool outermax, bool legacy, char **err)
{
    /* With binlegacy in hs or spe compatibility, the subcircuit expander kept only the cards of a
       subcircuit-local model (name x...:model) with min <= size < max in float on l and w/nf (subckt.c);
       the legacy rule then chose among those. */
    bool hs_pruned = legacy && (newcompat.hs || newcompat.spe) && strchr(name, ':');

    INPmodel    *modtmp, *best = NULL;
    int          best_match = 0;
    double       parse_values[4];
    bool         parse_found[4];
    static char *model_tokens[] = { "lmin", "lmax", "wmin", "wmax" };

    static int scan = -1;
    bin_list *cands = NULL;
    int ci = 0;

    if (scan < 0)
        scan = getenv("BINSEL_SCAN") ? 1 : 0;
    if (!scan) {
        if (!bin_index || bin_index_head != modtab || bin_index_modhash != modtabhash)
            bin_index_build();
        cands = (bin_list *) nghash_find(bin_index, name);
        if (!cands)
            return NULL;
    }

    for (modtmp = scan ? modtab : cands->v[0]; modtmp;
         modtmp = scan ? modtmp->INPnextModel : (++ci < cands->n ? cands->v[ci] : NULL)) {

        if (scan && model_name_match(name, modtmp->INPmodName) < 2)
            continue;

        /* skip if not binnable */
        if (!binnable_type(modtmp->INPmodType))
            continue;

        /* if illegal device type */
        if (modtmp->INPmodType < 0) {
            *err = tprintf("Unknown device type for model %s\n", name);
            return NULL;
        }

        if (cands && !binsel_nocache()) {
            char *ln = modtmp->INPmodLine->line;
            if (cands->state[ci] == 0 || cands->line[ci] != ln) {
                cands->state[ci] = parse_line(ln, model_tokens, 4, cands->lim[ci], parse_found) ? 1 : 2;
                cands->line[ci] = ln;
            }
            if (cands->state[ci] != 1)
                continue;
            memcpy(parse_values, cands->lim[ci], sizeof parse_values);
        } else if (!parse_line(modtmp->INPmodLine->line, model_tokens, 4, parse_values, parse_found)) {
            continue;
        }

        if (hs_pruned &&
            !((float) l >= (float) parse_values[0] && (float) l < (float) parse_values[1] &&
              (float) w_nf >= (float) parse_values[2] && (float) w_nf < (float) parse_values[3]))
            continue;

        int match = legacy
            ? INPbinMatchLegacy(l, parse_values[0], parse_values[1], w, parse_values[2], parse_values[3])
            : INPbinMatch(l, parse_values[0], parse_values[1], w, parse_values[2], parse_values[3], outermax);
        if (match == 0)
            continue;
        if (!best || match < best_match) {
            best = modtmp;
            best_match = match;
        }
        /* the standard rule holds: take the first card found, as before */
        if (match == 1)
            break;
    }
    return best;
}


char *
INPgetModBin(CKTcircuit *ckt, char *name, INPmodel **model, INPtables *tab, char *line)
{
    INPmodel    *best, *other;
    bool         outermax = INPbinOuterMax();
    bool         legacy = INPbinLegacy();
    char        *err = NULL;
    double       l, w, w_nf, nf_raw;
    double       parse_values[4];
    bool         parse_found[4];
    static char *instance_tokens[] = { "l", "w", "nf", "wnflag" };
    static CKTcircuit *warned_ckt = NULL;
    static int   n_edge_warned = 0;
    double       scale;
    int          wnflag = INPgetWnflag();

    if (!cp_getvar("scale", CP_REAL, &scale, 0))
        scale = 1;

    *model = NULL;

    /* read W and L. If not on the instance line, leave */
    if (!parse_line(line, instance_tokens, 2, parse_values, parse_found))
        return NULL;

    /* nf as given (1 if not), for the w/nf of the legacy subcircuit card filter */
    nf_raw = parse_line(line, instance_tokens, 3, parse_values, parse_found) ? parse_values[2] : 1.;
    if (nf_raw <= 0.)
        nf_raw = 1.;

    /* This is for reading nf. If nf is not available, set to 1 if in HSPICE or Spectre compatibility mode */
    if (!parse_line(line, instance_tokens, 3, parse_values, parse_found)) {
        parse_values[2] = 1.; /* divisor */
    }
    /* This is for reading wnflag from instance. If it is not available, no change.
       If instance wnflag == 0, set divisor to 1, else use instance nf */
    else if (parse_line(line, instance_tokens, 4, parse_values, parse_found)) {
        /* wnflag from instance overrules: no use of nf */
        if (parse_values[3] == 0) {
            parse_values[2] = 1.; /* divisor */
        }
    }
    /* We do have nf, but no wnflag on the instance. Now it depends on the default
       wnflag or on the .options wnflag */
    else {
        if (wnflag == 0)
            parse_values[2] = 1.; /* divisor */
    }


    l = parse_values[0] * scale;
    w = parse_values[1] / parse_values[2] * scale;
    w_nf = (float) scale * parse_values[1] / nf_raw;

    best = bin_select(name, l, w, w_nf, outermax, legacy, &err);
    if (err) {
        *model = NULL;
        return err;
    }

    /* #1766: a size exactly on a bin edge is held by two cards; the standard rule and the legacy rule
       (ngspice before 2026-10-30) may take different ones, and the two cards may simulate differently
       (a sync buck of sky130 5 V devices finished with one and stopped with "timestep too small" with the other).
       Say so, for the first instances of each circuit. */
    other = bin_select(name, l, w, w_nf, outermax, !legacy, &err);
    if (best && other && best != other) {
        if (warned_ckt != ckt) {
            warned_ckt = ckt;
            n_edge_warned = 0;
        }
        if (n_edge_warned < 5) {
            fprintf(stderr,
                    "Warning: model %s, l=%g w=%g%s: on a bin edge; card %s is used (%s rule), "
                    "the %s rule takes %s (%s)\n",
                    name, l, w, (parse_values[2] != 1.) ? " per finger" : "",
                    best->INPmodName, legacy ? "binlegacy" : "min <= size < max",
                    legacy ? "standard" : "pre-2026-10 ngspice", other->INPmodName,
                    legacy ? "without .option binlegacy" : ".option binlegacy");
        } else if (n_edge_warned == 5) {
            fprintf(stderr, "Warning: further instances on a bin edge are not listed\n");
        }
        n_edge_warned++;
    }

    if (best) {
        /* create unless model is already defined */
        if (!best->INPmodfast) {
            int error = create_model(ckt, best, tab);
            if (error)
                return NULL;
        }
        *model = best;
    }

    return NULL;
}


char *
INPgetMod(CKTcircuit *ckt, char *name, INPmodel **model, INPtables *tab)
{
    INPmodel *modtmp;

#ifdef TRACE
    printf("In INPgetMod, examining model %s ...\n", name);
#endif

    if (modtabhash) {
        modtmp = nghash_find(modtabhash, name);
        if (modtmp) {
            /* found the model in question - now instantiate if necessary */
            /* and return an appropriate pointer to it */

    /* if illegal device type */
            if (modtmp->INPmodType < 0) {
#ifdef TRACE
                printf("In INPgetMod, illegal device type for model %s ...\n", name);
#endif
                * model = NULL;
                return tprintf("Unknown device type for model %s\n", name);
            }

            /* create unless model is already defined */
            if (!modtmp->INPmodfast) {
                int error = create_model(ckt, modtmp, tab);
                if (error) {
                    *model = NULL;
                    return INPerror(error);
                }
            }

            *model = modtmp;
            return NULL;
        }
    }
#if (0)
    for (modtmp = modtab; modtmp; modtmp = modtmp->INPnextModel) {

#ifdef TRACE
        printf("In INPgetMod, comparing %s against stored model %s ...\n", name, modtmp->INPmodName);
#endif

        if (strcmp(modtmp->INPmodName, name) == 0) {
            /* found the model in question - now instantiate if necessary */
            /* and return an appropriate pointer to it */

            /* if illegal device type */
            if (modtmp->INPmodType < 0) {
#ifdef TRACE
                printf("In INPgetMod, illegal device type for model %s ...\n", name);
#endif
                *model = NULL;
                return tprintf("Unknown device type for model %s\n", name);
            }

            /* create unless model is already defined */
            if (!modtmp->INPmodfast) {
                int error = create_model(ckt, modtmp, tab);
                if (error) {
                    *model = NULL;
                    return INPerror(error);
                }
            }

            *model = modtmp;
            return NULL;
        }
    }
#endif
#ifdef TRACE
    printf("In INPgetMod, didn't find model for %s, using default ...\n", name);
#endif

    *model = NULL;
    return tprintf("Unable to find definition of model %s\n", name);
}


#ifdef CIDER
/*
 * Parse a numerical model by running through the list of original
 * input cards which make up the model
 * Given:
 * 1. First card looks like: .model modname modtype <level=val>
 * 2. Other cards look like: +<whitespace>? where ? tells us what
 * to do with the next card:
 *    '#$*' = comment card
 *    '+'   = continue previous card
 *    other = new card
 */
static int
INPparseNumMod(CKTcircuit *ckt, INPmodel *model, INPtables *tab, char **errMessage)
{
    struct card *txtCard;    /* Text description of a card */
    GENcard *tmpCard = NULL; /* Processed description of a card */
    IFcardInfo *info = NULL; /* Info about the type of card located */
    char *cardName = NULL;   /* name of a card */
    int cardNum = 0;         /* number of this card in the overall line */
    char *err = NULL;        /* Strings for error messages */
    int error;

    /* Chase down to the top of the list of actual cards */
    txtCard = model->INPmodLine->actualLine;

    /* Skip the first card if it exists since there's nothing interesting */
    /* txtCard will be empty if the numerical model is empty */
    if (txtCard)
        txtCard = txtCard->nextcard;

    /* Now parse each remaining card */
    for (; txtCard; txtCard = txtCard->nextcard) {
        char *line = txtCard->line;
        cardNum++;

        /* Skip the initial '+' and any whitespace. */
        line++;
        while (*line == ' ' || *line == '\t')
            line++;

        switch (*line) {
        case '*':
        case '$':
        case '#':
        case '\0':
        case '\n':
            /* comment or empty cards */
            info = NULL;
            continue;
        case '+':
            /* continuation card */
            if (!info) {
                err = INPerrCat(err,
                                tprintf("Error on card %d : illegal continuation \'+\' - ignored",
                                        cardNum));
                continue;
            }
            /* Skip leading '+'s */
            while (*line == '+')
                line++;
            break;
        default:
            info = NULL;
            break;
        }

        if (!info) {
            /* new command card */
            if (cardName)       /* get rid of old card name */
                FREE(cardName);
            INPgetTok(&line, &cardName, 1);        /* get new card name */
            if (*cardName) {                 /* Found a name? */
                int lastType = INPfindCard(cardName, INPcardTab, INPnumCards);
                if (lastType >= 0) {
                    /* Add card structure to model */
                    info = INPcardTab[lastType];
                    error = info->newCard(&tmpCard, model->INPmodfast);
                    if (error) {
                        FREE(cardName);
                        return error;
                    }
                    /* Handle parameter-less cards */
                } else if (cinprefix(cardName, "title", 3)) {
                    /* Do nothing */
                } else if (cinprefix(cardName, "comment", 3)) {
                    /* Do nothing */
                } else if (cinprefix(cardName, "end", 3)) {
                    /* Terminate parsing */
                    *errMessage = err;
                    FREE(cardName);
                    return 0;
                } else {
                    /* Error */
                    err = INPerrCat(err,
                                    tprintf("Error on card %d : unrecognized name (%s) - ignored",
                                            cardNum, cardName));
                }
                FREE(cardName);
            }
        }

        if (!info)
            continue;

        /* parse the rest of this line */
        while (*line) {

            int invert = FALSE;
            /* Strip leading carat from booleans */
            if (*line == '^') {
                invert = TRUE;
                line++;
            }

            char *parm;                /* name of a parameter */
            INPgetTok(&line, &parm, 1);
            if (!*parm) {
                FREE(parm);
                break;
            }

            int idx = INPfindParm(parm, info->cardParms, info->numParms);
            if (idx == E_MISSING) {
                /* parm not found */
                err = INPerrCat(err,
                                tprintf("Error on card %d : unrecognized parameter (%s) - ignored",
                                        cardNum, parm));
            } else if (idx == E_AMBIGUOUS) {
                /* parm ambiguous */
                err = INPerrCat(err,
                                tprintf("Error on card %d : ambiguous parameter (%s) - ignored",
                                        cardNum, parm));
            } else {
                IFvalue *value = INPgetValue(ckt, &line, info->cardParms[idx].dataType, tab);

                /* invert if this is a boolean entry */
                if (invert) {
                    if ((info->cardParms[idx].dataType & IF_VARTYPES) == IF_FLAG)
                        value->iValue = 0;
                    else
                        err = INPerrCat(err,
                                        tprintf("Error on card %d : non-boolean parameter (%s) - \'^\' ignored",
                                                cardNum, parm));
                }

                error = info->setCardParm(info->cardParms[idx].id, value, tmpCard);
                if (info->cardParms[idx].dataType & IF_STRING) {
                    FREE(value->sValue);
                } else if (info->cardParms[idx].dataType & IF_REALVEC) {
                    FREE(value->v.vec.rVec);
                } else if (info->cardParms[idx].dataType & IF_INTVEC) {
                    FREE(value->v.vec.iVec);
                }
                if (error)
                    return error;
            }
            FREE(parm);
        }
    }

    *errMessage = err;
    return 0;
}


/*
 * Locate the best match to a card name in an IFcardInfo table
 */
static int
INPfindCard(char *name, IFcardInfo *table[], int numCards)
{
    int length = (int) strlen(name);
    int best = E_MISSING;
    int bestMatch = 0;

    int test;

    /* compare all the names in the card table to this name */
    for (test = 0; test < numCards; test++) {
        int match = cimatch(name, table[test]->name);
        if ((match > 0) && (match == bestMatch)) {
            best = E_AMBIGUOUS;
        } else if ((match > bestMatch) && (match == length)) {
            best = test;
            bestMatch = match;
        }
    }

    return best;
}


/*
 * Locate the best match to a parameter name in an IFparm table
 */
static int
INPfindParm(char *name, IFparm *table, int numParms)
{
    int length = (int) strlen(name);
    int best = E_MISSING;
    int bestMatch = 0;
    int bestId = -1;

    int test;

    /* compare all the names in the parameter table to this name */
    for (test = 0; test < numParms; test++) {
        int match = cimatch(name, table[test].keyword);
        if ((match == length) && (match == (int) strlen(table[test].keyword))) {
            /* exact match */
            return test;
        }
        int id = table[test].id;
        if ((match > 0) && (match == bestMatch) && (id != bestId)) {
            best = E_AMBIGUOUS;
        } else if ((match > bestMatch) && (match == length)) {
            bestMatch = match;
            bestId = id;
            best = test;
        }
    }

    return best;
}

#endif /* CIDER */
