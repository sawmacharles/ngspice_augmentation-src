/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include <stdio.h>
#include <stdint.h>
#include "ngspice/ifsim.h"
#include "ngspice/inpmacs.h"
#include "ngspice/fteext.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/inpptree.h"
#include "inpxx.h"
#include "ngspice/cktdefs.h"

/* Uncomment to allow tracing */
/*#define TRACE*/


extern double PTfudge_factor;

/* accel candidate 0104 (round 3, card #483 T9.3): per-IFeval()-call node
 * cache, kept OUTSIDE INPparseNode -- see inpptree.h's INPparseTree
 * comment for why (round 2's 0101 grew INPparseNode itself and correlated
 * with an a_pex_buck regression). This is a small open-addressing hash
 * table (linear probe, no tombstones) keyed by INPparseNode pointer,
 * generation-stamped exactly like 0101's per-node fields were: a slot is
 * "live" only if its gen matches the table's current generation, which is
 * bumped once per IFeval() call (not once per PTeval() call), so the
 * value walk and every one of the numVars derivative walks inside the
 * SAME IFeval() invocation share one cache -- a node reached from more
 * than one of those walks (common: PTdifferentiate() deliberately shares
 * subtree pointers between the value tree and its derivative trees, see
 * its comment in inpptree.c) is computed once and reused.
 *
 * Correctness of the open-addressing scheme: lookup and insert use the
 * SAME deterministic probe sequence starting at hash(key) (linear
 * probing). A slot is "available" for a fresh insert as soon as its
 * stamped generation differs from the current one (no separate
 * empty/tombstone state needed, since nothing is ever deleted mid-
 * generation). Lookup stops and reports a HIT at the first slot with
 * (gen == current && key == tree); it stops and reports a MISS at the
 * first slot with (gen != current), because insertion always claims the
 * first available slot on its own probe from hash(key) -- if `tree` had
 * been inserted earlier this generation, lookup (using the identical
 * probe order) would have reached that slot no later than reaching any
 * later available slot. A slot with (gen == current && key != tree) is a
 * live entry for a DIFFERENT key (a hash collision this generation) and
 * must not stop either search; both continue probing past it. */
typedef struct PTcacheSlot {
    INPparseNode *key;
    unsigned long long gen;
    double val;
} PTcacheSlot;

struct PTprog;
typedef struct PTprog PTprog;
typedef struct PTctl {
    PTcacheSlot *slots;  /* hash cache of the tree walk (allocated only if needed) */
    PTprog *prog;
    int state;           /* 0 not tried, 1 flat program, 2 tree walk */
} PTctl;

/* Rough (non-deduplicated -- a safe over-count, not exact) upper bound on
 * the number of distinct INPparseNode pointers reachable from `tree`,
 * used only to size the cache table once, lazily, on the first IFeval()
 * call for a given INPparseTree. Cheap and one-time: not on the hot path
 * (never called again after the table exists), so a plain recursive walk
 * with no memoization of its own is fine even though it may revisit
 * shared subtrees once per parent that shares them. */
static int
PTcountNodes(INPparseNode *tree)
{
    if (!tree)
        return 0;
    switch (tree->type) {
    case PT_FUNCTION:
        switch (tree->funcnum) {
        case PTF_POW:
        case PTF_PWR:
        case PTF_MIN:
        case PTF_MAX:
            return 1 + PTcountNodes(tree->left->left) + PTcountNodes(tree->left->right);
        default:
            return 1 + PTcountNodes(tree->left);
        }
    case PT_TERN:
        return 1 + PTcountNodes(tree->left)
                 + PTcountNodes(tree->right->left)
                 + PTcountNodes(tree->right->right);
    case PT_PLUS:
    case PT_MINUS:
    case PT_TIMES:
    case PT_DIVIDE:
    case PT_POWER:
        return 1 + PTcountNodes(tree->left) + PTcountNodes(tree->right);
    default:
        return 1;
    }
}

/* Lazily allocate and size myTree's cache table on its first use. Called
 * once per INPparseTree, ever (myTree->cache is checked by the caller and
 * is non-NULL on every later call). Table size: next power of two >= 4x
 * the (over-counted) node total, minimum 16 slots, so load factor stays
 * low (<=25% at allocation time) and probe chains stay short even though
 * PTcountNodes() over-counts shared nodes. */
static void
PTcacheAlloc(INPparseTree *myTree)
{
    PTctl *ctl = (PTctl *) myTree->cache;
    int i;
    int n = PTcountNodes(myTree->tree);

    for (i = 0; i < myTree->p.numVars; i++)
        n += PTcountNodes(myTree->derivs[i]);

    unsigned int size = 16;
    while ((int) size < n * 4)
        size <<= 1;

    ctl->slots = TMALLOC(PTcacheSlot, size);   /* zero-initialized: gen==0 everywhere */
    myTree->cache_mask = size - 1;
    myTree->cache_gen = 0;   /* first real generation issued is 1, see below */
}

/* FNV-1a-style pointer hash-and-mix, cheap (a handful of ALU ops), no
 * division. Node pointers are tmalloc()'d (effectively 16-byte-or-better
 * aligned on any platform this fork targets) so the low bits are mostly
 * zero and are mixed away rather than relied on directly. */
static inline unsigned int
PThashPtr(const void *p, unsigned int mask)
{
    uintptr_t x = (uintptr_t) p;
    x ^= x >> 15;
    x *= (uintptr_t) 0x9E3779B97F4A7C15ULL;
    x ^= x >> 13;
    return (unsigned int) (x & mask);
}


/* card #532 (round 15): flat B-source evaluator.  The value tree and all the
 * derivative trees of one INPparseTree form a DAG (PTdifferentiate shares
 * subtrees).  PTeval walks it recursively and memoizes per node in a hash
 * table.  For ternary-free trees the whole DAG is compiled once into a linear
 * instruction list in the order in which the memoized walk first evaluates each
 * node (value tree, then derivs[0], derivs[1], ...: DFS post-order, left before
 * right), one register per distinct node.  Every instruction performs exactly
 * the operation PTeval performs for that node (same function pointers, same
 * HUGE checks, same error text), so values, derivatives and error behaviour
 * are bit-identical.  Trees with PT_TERN (lazy branches) or any node type not
 * handled here keep the tree walk.  NGSPICE_FLATB=0 disables.  */
enum { FB_CONST, FB_VAR, FB_TIME, FB_TEMP, FB_FREQ, FB_ADD, FB_SUB, FB_MUL, FB_DIV,
       FB_BIN, FB_UN, FB_UNP };
typedef struct FBinsn {
    int op, a, b;
    void (*fn)(void);
    INPparseNode *node;
} FBinsn;
struct PTprog {
    int n;               /* instructions == registers */
    int nroots;          /* 1 + numVars */
    int *root;           /* register of each root (value first) */
    FBinsn *code;
    double *reg;
};

typedef struct { INPparseNode **key; int *idx; unsigned int mask; int cap; PTprog *p; int bad; } FBcomp;

static int
FBcompile(FBcomp *c, INPparseNode *t)
{
    FBinsn in;
    unsigned int h;
    int a = -1, b = -1;
    if (c->bad || !t) { c->bad = 1; return -1; }
    h = PThashPtr(t, c->mask);
    while (c->key[h]) {
        if (c->key[h] == t) return c->idx[h];
        h = (h + 1) & c->mask;
    }
    memset(&in, 0, sizeof in);
    in.node = t;
    switch (t->type) {
    case PT_CONSTANT: in.op = FB_CONST; break;
    case PT_VAR: in.op = FB_VAR; in.a = t->valueIndex; break;
    case PT_TIME: in.op = FB_TIME; break;
    case PT_TEMPERATURE: in.op = FB_TEMP; break;
    case PT_FREQUENCY: in.op = FB_FREQ; break;
    case PT_FUNCTION:
        switch (t->funcnum) {
        case PTF_POW: case PTF_PWR: case PTF_MIN: case PTF_MAX:
            a = FBcompile(c, t->left->left);
            b = FBcompile(c, t->left->right);
            in.op = FB_BIN; break;
        default:
            a = FBcompile(c, t->left);
            in.op = t->data ? FB_UNP : FB_UN; break;
        }
        in.fn = t->function;
        break;
    case PT_PLUS: case PT_MINUS: case PT_TIMES: case PT_DIVIDE: case PT_POWER:
        a = FBcompile(c, t->left);
        b = FBcompile(c, t->right);
        in.op = t->type == PT_PLUS ? FB_ADD : t->type == PT_MINUS ? FB_SUB :
                t->type == PT_TIMES ? FB_MUL : t->type == PT_DIVIDE ? FB_DIV : FB_BIN;
        in.fn = t->function;
        break;
    default:   /* PT_TERN (lazy), PT_PLACEHOLDER, ...: not handled */
        c->bad = 1;
        return -1;
    }
    if (c->bad) return -1;
    if (c->p->n >= c->cap) { c->bad = 1; return -1; }
    in.a = (in.op == FB_VAR) ? in.a : a;
    in.b = b;
    c->p->code[c->p->n] = in;
    h = PThashPtr(t, c->mask);
    while (c->key[h]) h = (h + 1) & c->mask;
    c->key[h] = t; c->idx[h] = c->p->n;
    return c->p->n++;
}

static PTprog *
FBbuild(INPparseTree *myTree)
{
    int i, n = PTcountNodes(myTree->tree);
    FBcomp c;
    PTprog *p;
    unsigned int size = 16;
    for (i = 0; i < myTree->p.numVars; i++)
        n += PTcountNodes(myTree->derivs[i]);
    while ((int) size < n * 4) size <<= 1;
    p = TMALLOC(PTprog, 1);
    p->code = TMALLOC(FBinsn, n + 1);
    p->reg = TMALLOC(double, n + 1);
    p->nroots = 1 + myTree->p.numVars;
    p->root = TMALLOC(int, p->nroots);
    p->n = 0;
    c.key = TMALLOC(INPparseNode *, size);
    c.idx = TMALLOC(int, size);
    c.mask = size - 1; c.cap = n + 1; c.p = p; c.bad = 0;
    p->root[0] = FBcompile(&c, myTree->tree);
    for (i = 0; i < myTree->p.numVars && !c.bad; i++)
        p->root[1 + i] = FBcompile(&c, myTree->derivs[i]);
    free(c.key); free(c.idx);
    if (c.bad) {
        free(p->code); free(p->reg); free(p->root); free(p);
        return NULL;
    }
    for (i = 0; i < p->n; i++)
        if (p->code[i].op == FB_CONST) p->reg[i] = p->code[i].node->constant;
    return p;
}

/* returns OK or the error code of the first failing instruction; *fail is the
 * index of that instruction (so that roots already computed can be stored) */
static int
FBrun(PTprog *p, double gmin, const double *vals, double *r, int *fail)
{
    double *R = p->reg;
    const FBinsn *code = p->code;
    int i, n = p->n;
    PTfudge_factor = gmin * 1.0e-20;
    for (i = 0; i < n; i++) {
        const FBinsn *in = &code[i];
        double r1, r2, res;
        switch (in->op) {
        case FB_CONST: continue;
        case FB_VAR: R[i] = vals[in->a]; continue;
        case FB_TIME: R[i] = ((CKTcircuit *) in->node->data)->CKTtime; continue;
        case FB_TEMP: R[i] = ((CKTcircuit *) in->node->data)->CKTtemp - CONSTCtoK; continue;
        case FB_FREQ: R[i] = (((CKTcircuit *) in->node->data)->CKTomega) / 2. / M_PI; continue;
        case FB_UN:
            r1 = R[in->a];
            res = PTunary(in->fn)(r1);
            if (res == HUGE) {
                fprintf(stderr, "Error: %g out of range for %s\n", r1, in->node->funcname);
                *fail = i; return E_PARMVAL;
            }
            R[i] = res; continue;
        case FB_UNP:
            r1 = R[in->a];
            res = PTunary_with_private(in->fn)(r1, in->node->data);
            if (res == HUGE) {
                fprintf(stderr, "Error: %g out of range for %s\n", r1, in->node->funcname);
                *fail = i; return E_PARMVAL;
            }
            R[i] = res; continue;
        default:
            r1 = R[in->a]; r2 = R[in->b];
            switch (in->op) {
            case FB_ADD: res = r1 + r2; break;
            case FB_SUB: res = r1 - r2; break;
            case FB_MUL: res = r1 * r2; break;
            case FB_DIV: res = PTbinary(in->fn)(r1, r2); break;
            default:     res = PTbinary(in->fn)(r1, r2); break;
            }
            if (res == HUGE) {
                if (in->node->type == PT_FUNCTION)
                    fprintf(stderr, "Error: %g, %g out of range for %s\n", r1, r2, in->node->funcname);
                else
                    fprintf(stderr, "\nError: %g, %g out of range for %s\n", r1, r2, in->node->funcname);
                *fail = i; return E_PARMVAL;
            }
            R[i] = res; continue;
        }
    }
    (void) r;
    return OK;
}

static int fb_enabled = -1;
static int fb_perturb = 0;   /* test-only red-proof: NGSPICE_FLATB=2 scales the flat result by 1+1e-9 */

static int PTeval(INPparseNode * tree, double gmin, double *res,
		  double *vals, PTcacheSlot *cache, unsigned int cache_mask,
		  unsigned long long cache_gen);



int
IFeval(IFparseTree * tree, double gmin, double *result, double *vals,
       double *derivs)
{
    int i, err;
    INPparseTree *myTree = (INPparseTree *) tree;
    PTcacheSlot *cache;

    if (!myTree) {
        fprintf(stderr, "\nInternal error: No tree to evaluate.\n");
        controlled_exit(EXIT_BAD);
    }

    if (!myTree->cache)
        myTree->cache = TMALLOC(PTctl, 1);
    {
        PTctl *ctl = (PTctl *) myTree->cache;
        if (ctl->state == 0) {
            if (fb_enabled < 0) {
                const char *e = getenv("NGSPICE_FLATB");
                fb_enabled = !(e && atoi(e) == 0);
                fb_perturb = (e && atoi(e) == 2);
            }
            ctl->prog = fb_enabled ? FBbuild(myTree) : NULL;
            ctl->state = ctl->prog ? 1 : 2;
        }
        if (ctl->state == 1) {
            PTprog *p = ctl->prog;
            int fail = p->n;
            err = FBrun(p, gmin, vals, NULL, &fail);
            if (err == OK) {
                *result = p->reg[p->root[0]];
                for (i = 0; i < myTree->p.numVars; i++)
                    derivs[i] = p->reg[p->root[1 + i]];
                if (fb_perturb)
                    *result *= 1.0 + 1e-9;
#ifdef TRACE
                printf("results: function = %lg\n", *result);
#endif
                return OK;
            }
            /* same partial outputs as the tree walk: roots that were complete before the failing node */
            if (p->root[0] < fail)
                *result = p->reg[p->root[0]];
            for (i = 0; i < myTree->p.numVars; i++)
                if (p->root[1 + i] < fail)
                    derivs[i] = p->reg[p->root[1 + i]];
            if (ft_ngdebug) {
                INPptPrint("calling PTeval, tree = ", tree);
                printf("values:");
                for (i = 0; i < myTree->p.numVars; i++)
                    printf("\tvar%d = %lg\n", i, vals[i]);
            }
            if (ft_stricterror)
                controlled_exit(EXIT_BAD);
            return err;
        }
        if (!ctl->slots)
            PTcacheAlloc(myTree);
        cache = ctl->slots;
    }
    myTree->cache_gen++;   /* generation 0 is reserved (never issued), so a
                             * freshly zero-initialized slot never reads
                             * back as a spurious hit. */

#ifdef TRACE
    INPptPrint("calling PTeval, tree = ", tree);
    printf("values:");
    for (i = 0; i < myTree->p.numVars; i++)
	printf("\tvar%d = %lg\n", i, vals[i]);
#endif

    if ((err = PTeval(myTree->tree, gmin, result, vals,
                       cache, myTree->cache_mask, myTree->cache_gen)) != OK) {
        if (ft_ngdebug) {
            INPptPrint("calling PTeval, tree = ", tree);
            printf("values:");
            for (i = 0; i < myTree->p.numVars; i++)
	        printf("\tvar%d = %lg\n", i, vals[i]);
        }
        if (ft_stricterror)
            controlled_exit(EXIT_BAD);
        return err;
    }

    for (i = 0; i < myTree->p.numVars; i++)
        if ((err = PTeval(myTree->derivs[i], gmin, &derivs[i], vals,
                           cache, myTree->cache_mask, myTree->cache_gen)) != OK) {
            if (ft_ngdebug) {
                INPptPrint("calling PTeval, tree = ", tree);
                printf("results: function = %lg\n", *result);
                for (i = 0; i < myTree->p.numVars; i++)
	                printf("\td / d var%d = %lg\n", i, derivs[i]);
            }
            if (ft_stricterror)
                controlled_exit(EXIT_BAD);
            return err;
        }

#ifdef TRACE
    printf("results: function = %lg\n", *result);
    for (i = 0; i < myTree->p.numVars; i++)
	printf("\td / d var%d = %lg\n", i, derivs[i]);
#endif

    return (OK);
}

static int
PTeval(INPparseNode * tree, double gmin, double *res, double *vals,
       PTcacheSlot *cache, unsigned int cache_mask, unsigned long long cache_gen)
{
    double r1, r2;
    int err;
    /* accel candidate 0104: only the node types below that actually
     * recurse (PT_FUNCTION/PT_TERN/the binary arithmetic ops) are worth
     * caching -- PT_CONSTANT/PT_VAR/PT_TIME/PT_TEMPERATURE/PT_FREQUENCY
     * are already O(1) leaf reads, and a hash probe costs more than just
     * doing them again. `slot` is set (non-NULL) only when this call did
     * a cache lookup and found a miss, so the tail of the function knows
     * where to memoize the result; it stays NULL (no memoization) for the
     * cheap leaf types. */
    PTcacheSlot *slot = NULL;

    switch (tree->type) {
    case PT_FUNCTION:
    case PT_TERN:
    case PT_PLUS:
    case PT_MINUS:
    case PT_TIMES:
    case PT_DIVIDE:
    case PT_POWER:
        if (cache) {
            unsigned int h = PThashPtr(tree, cache_mask);
            for (;;) {
                PTcacheSlot *s = &cache[h];
                if (s->gen == cache_gen) {
                    if (s->key == tree) {
                        *res = s->val;
                        return (OK);
                    }
                    /* live entry this generation, but a different key
                     * (hash collision) -- keep probing, see the
                     * correctness note above PTcacheSlot. */
                } else {
                    /* first available slot on this probe sequence: `tree`
                     * was not inserted earlier this generation (see the
                     * correctness note) -- this is where it goes once
                     * computed below. */
                    slot = s;
                    break;
                }
                h = (h + 1) & cache_mask;
            }
        }
        break;
    default:
        break;
    }

    PTfudge_factor = gmin * 1.0e-20; /* defaults to 1e-32, should be small enough */
    switch (tree->type) {
    case PT_CONSTANT:
	*res = tree->constant;
	break;

    case PT_VAR:
	*res = vals[tree->valueIndex];
	break;

    case PT_FUNCTION:
        switch(tree->funcnum) {

        case PTF_POW:
        case PTF_PWR:
        case PTF_MIN:
        case PTF_MAX:
            err = PTeval(tree->left->left, gmin, &r1, vals, cache, cache_mask, cache_gen);
            if (err != OK)
                return (err);
            err = PTeval(tree->left->right, gmin, &r2, vals, cache, cache_mask, cache_gen);
            if (err != OK)
                return (err);
            *res = PTbinary(tree -> function) (r1, r2);
            if (*res == HUGE) {
                fprintf(stderr, "Error: %g, %g out of range for %s\n",
                r1, r2, tree->funcname);
                return (E_PARMVAL);
            }
        break;
        /* fcns with single argument */
        default:
            err = PTeval(tree->left, gmin, &r1, vals, cache, cache_mask, cache_gen);
            if (err != OK)
                return (err);
            if(tree->data == NULL)
                *res = PTunary(tree -> function) (r1);
            else
                *res = PTunary_with_private(tree -> function) (r1, tree->data);
            if (*res == HUGE) {
                fprintf(stderr, "Error: %g out of range for %s\n",
                r1, tree->funcname);
                return (E_PARMVAL);
            }
        break;
        }
	break;

    case PT_TERN:
      {
        INPparseNode *arg1 = tree->left;
        INPparseNode *arg2 = tree->right->left;
        INPparseNode *arg3 = tree->right->right;

        err = PTeval(arg1, gmin, &r1, vals, cache, cache_mask, cache_gen);
        if (err != OK)
          return (err);

        /*FIXME > 0.0, >= 0.5, != 0.0 or what ? */
        err = PTeval((r1 != 0.0) ? arg2 : arg3, gmin, &r2, vals, cache, cache_mask, cache_gen);
        if (err != OK)
           return (err);

        *res = r2;
        break;
      }

    case PT_PLUS:
    case PT_MINUS:
    case PT_TIMES:
    case PT_DIVIDE:
    case PT_POWER:
	err = PTeval(tree->left, gmin, &r1, vals, cache, cache_mask, cache_gen);
	if (err != OK)
	    return (err);
	err = PTeval(tree->right, gmin, &r2, vals, cache, cache_mask, cache_gen);
	if (err != OK)
	    return (err);
	*res = PTbinary(tree -> function) (r1, r2);
	if (*res == HUGE) {
	    fprintf(stderr, "\nError: %g, %g out of range for %s\n",
		    r1, r2, tree->funcname);
	    return (E_PARMVAL);
	}
	break;

    case PT_TIME:
        *res = ((CKTcircuit*) tree->data) -> CKTtime;
	break;

    case PT_TEMPERATURE:
        *res = ((CKTcircuit*) tree->data) -> CKTtemp - CONSTCtoK;
	break;

    case PT_FREQUENCY:
        *res = (((CKTcircuit*) tree->data) -> CKTomega)/2./M_PI;
	break;
	
    default:
	fprintf(stderr, "Internal Error: bad node type %d\n", tree->type);
	return (E_PANIC);
    }

    /* accel candidate 0104: only reached on a successful evaluation of a
     * cacheable node type (every error path above returns early) -- if
     * this call found a miss and claimed an available slot up top,
     * memoize the result there now for the rest of this IFeval() call. */
    if (slot) {
        slot->key = tree;
        slot->gen = cache_gen;
        slot->val = *res;
    }

    return (OK);
}
