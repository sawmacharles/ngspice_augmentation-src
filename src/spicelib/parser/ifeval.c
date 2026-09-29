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
    int i;
    int n = PTcountNodes(myTree->tree);

    for (i = 0; i < myTree->p.numVars; i++)
        n += PTcountNodes(myTree->derivs[i]);

    unsigned int size = 16;
    while ((int) size < n * 4)
        size <<= 1;

    myTree->cache = TMALLOC(PTcacheSlot, size);   /* zero-initialized: gen==0 everywhere */
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
        PTcacheAlloc(myTree);
    myTree->cache_gen++;   /* generation 0 is reserved (never issued), so a
                             * freshly zero-initialized slot never reads
                             * back as a spurious hit. */
    cache = (PTcacheSlot *) myTree->cache;

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
