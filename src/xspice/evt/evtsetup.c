/*============================================================================
FILE    EVTsetup.c

MEMBER OF process XSPICE

Public Domain

Georgia Tech Research Corporation
Atlanta, Georgia 30332
PROJECT A-8503

AUTHORS

    9/12/91  Bill Kuhn

MODIFICATIONS

    <date> <person name> <nature of modifications>

SUMMARY

    This file contains function EVTsetup which clears/allocates the
    event-driven queues and data structures immediately prior to a new
    analysis.  In addition, it places entries in the job list so that
    results data from multiple analysis can be retrieved similar to
    SPICE3C1 saving multiple 'plots'.

INTERFACES

    int EVTsetup(CKTcircuit *ckt)

REFERENCED FILES

    None.

NON-STANDARD FEATURES

    None.

============================================================================*/

#include "ngspice/ngspice.h"
//#include "misc.h"

#include "ngspice/cktdefs.h"
#include "ngspice/sperror.h"
//#include "util.h"

#include "ngspice/mif.h"
#include "ngspice/evt.h"
#include "ngspice/evtudn.h"
#include "ngspice/mifproto.h"
#include "ngspice/evtproto.h"
#include "ngspice/devdefs.h"
#include "ngspice/cmtypes.h"   /* card #714: Digital_State_t */
#include "ngspice/mifparse.h"   /* card #714: Mif_Param_Info_t */
#include <math.h>




static int EVTsetup_queues(CKTcircuit *ckt);
static int EVTsetup_data(CKTcircuit *ckt);
static int EVTsetup_jobs(CKTcircuit *ckt);
static int EVTsetup_load_ptrs(CKTcircuit *ckt);

int EVTsetup_plot(CKTcircuit* ckt, char* plottypename);
int EVTswitch_plot(CKTcircuit* ckt, const char* plottypename);


/* Allocation macros with built-in check for out-of-memory */
/* Adapted from SPICE 3C1 code in CKTsetup.c */

#define CKALLOC(var,size,type) \
    if(size) { \
        if((var = TMALLOC(type, size)) == NULL) \
            return(E_NOMEM); \
    }

#define CKREALLOC(var,size,type) \
    if((size) == 1) { \
        if((var = TMALLOC(type, size)) == NULL) \
            return(E_NOMEM); \
    } else if((size) > 1) { \
        if((var = TREALLOC(type, (var), size)) == NULL) \
            return(E_NOMEM); \
    }


/*
EVTsetup

This function clears/allocates the event-driven queues and data structures
immediately prior to a new analysis.  In addition, it places entries
in the job list so that results data from multiple analysis can be
retrieved similar to SPICE3C1 saving multiple 'plots'.
*/


int EVTsetup(
    CKTcircuit *ckt)   /* The circuit structure */
{

    int err;


    /* Exit immediately if no event-driven instances in circuit */
    if(ckt->evt->counts.num_insts == 0)
        return(OK);

    /* Clear the inst, node, and output queues, and initialize the to_call */
    /* elements in the instance queue to call all event-driven instances */
    err = EVTsetup_queues(ckt);
    if(err)
        return(err);

    /* Allocate and initialize the node, state, message, and statistics data */
    err = EVTsetup_data(ckt);
    if(err)
        return(err);

    /* Set the job pointers to the allocated results, states, messages, */
    /* and statistics so that data will be accessable after run */
    err = EVTsetup_jobs(ckt);
    if(err)
        return(err);

    /* Setup the pointers in the MIFinstance structure for inputs, outputs, */
    /* and total loads */
    err = EVTsetup_load_ptrs(ckt);
    if(err)
        return(err);

    /* Initialize additional event data */
    g_mif_info.circuit.evt_step = 0.0;

    /* Return OK */
    return(OK);
}


int EVTunsetup(
    CKTcircuit* ckt)   /* The circuit structure */
{
    int err;

    /* Exit immediately if no event-driven instances in circuit */
    if (ckt->evt->counts.num_insts == 0)
        return(OK);

    /* Clear the inst, node, and output queues, and initialize the to_call */
    /* elements in the instance queue to call all event-driven instances */
    err = EVTsetup_queues(ckt);
    if (err)
        return(err);

    /* Initialize additional event data */
    g_mif_info.circuit.evt_step = 0.0;

    /* Return OK */
    return(OK);
}




/*
EVTsetup_queues

This function clears the event-driven queues in preparation for
a new simulation.
*/


static int EVTsetup_queues(
    CKTcircuit *ckt)         /* The circuit structure */
{

    int i;
    int num_insts;
    int num_nodes;
    int num_outputs;

    Evt_Inst_Queue_t    *inst_queue;
    Evt_Node_Queue_t    *node_queue;
    Evt_Output_Queue_t  *output_queue;

    Evt_Inst_Event_t    *inst_event;
    Evt_Output_Event_t  *output_event;
    void                *ptr;

    /* ************************ */
    /* Clear the instance queue */
    /* ************************ */

    num_insts = ckt->evt->counts.num_insts;
    inst_queue = &(ckt->evt->queue.inst);

    for(i = 0; i < num_insts; i++) {
        inst_event = inst_queue->head[i];
        while(inst_event) {
            ptr = inst_event;
            inst_event = inst_event->next;
            FREE(ptr);
        }
        inst_event = inst_queue->free[i];
        while(inst_event) {
            ptr = inst_event;
            inst_event = inst_event->next;
            FREE(ptr);
        }
        inst_queue->head[i] = NULL;
        inst_queue->current[i] = &(inst_queue->head[i]);
        inst_queue->last_step[i] = &(inst_queue->head[i]);
        inst_queue->free[i] = NULL;
    }

    inst_queue->next_time = 0.0;
    inst_queue->last_time = 0.0;

    inst_queue->num_modified = 0;
    inst_queue->num_pending = 0;
    inst_queue->num_to_call = 0;

    for(i = 0; i < num_insts; i++) {
        inst_queue->modified[i] = MIF_FALSE;
        inst_queue->pending[i] = MIF_FALSE;
        inst_queue->to_call[i] = MIF_FALSE;
    }


    /* ******************** */
    /* Clear the node queue */
    /* ******************** */

    num_nodes = ckt->evt->counts.num_nodes;
    node_queue = &(ckt->evt->queue.node);

    node_queue->num_changed = 0;
    node_queue->num_to_eval = 0;

    for(i = 0; i < num_nodes; i++) {
        node_queue->changed[i] = MIF_FALSE;
        node_queue->to_eval[i] = MIF_FALSE;
    }

    /* ********************** */
    /* Clear the output queue */
    /* ********************** */

    num_outputs = ckt->evt->counts.num_outputs;
    output_queue = &(ckt->evt->queue.output);

    for(i = 0; i < num_outputs; i++) {
        output_event = output_queue->head[i];
        while(output_event) {
            ptr = output_event;
            output_event = output_event->next;
            FREE(ptr);
        }
        output_queue->head[i] = NULL;
        output_queue->current[i] = &(output_queue->head[i]);
        output_queue->last_step[i] = &(output_queue->head[i]);
    }

    output_queue->next_time = 0.0;
    output_queue->last_time = 0.0;

    output_queue->num_modified = 0;
    output_queue->num_pending = 0;
    output_queue->num_changed = 0;

    if (num_outputs > 0) {
        for (i = 0; i < num_outputs; i++) {
            output_queue->modified[i] = MIF_FALSE;
            output_queue->pending[i] = MIF_FALSE;
            output_queue->changed[i] = MIF_FALSE;
        }

        if (output_queue->free_list[0]) {
            Evt_purge_free_outputs();
        } else {
            Evt_Output_Info_t *output_info;
            Evt_Node_Info_t   *node;

            /* ********************************************************* *
             * On first call for this circuit, set the free-list pointer
             * for each output queue.
             * ********************************************************* */

            output_info = ckt->evt->info.output_list;
            for (i = 0; i < num_outputs; i++) {
                node =  ckt->evt->info.node_table[output_info->node_index];
                output_queue->free_list[i] =
                    &g_evt_udn_info[node->udn_index]->free_list;
                output_info = output_info->next;
            }
        }
    }
    return OK;
}

void Evt_purge_free_outputs(void)
{
    Evt_Output_Event_t  *output_event, *next;
    int                  i;

    for (i = 0; i < g_evt_num_udn_types; ++i) {
        output_event = g_evt_udn_info[i]->free_list;
        g_evt_udn_info[i]->free_list = NULL;
        while (output_event) {
            next = output_event->next;
            tfree(output_event->value);
            tfree(output_event);
            output_event = next;
        }
    }
}


/*
EVTsetup_data

This function sets up the event-driven node, state, and
message data runtime structures in preparation for
a new simulation.
*/


static int EVTsetup_data(
    CKTcircuit *ckt)       /* The circuit structure */
{

    Evt_Data_t  *data;

    int         i;
    int         j;

    int         num_insts;
    int         num_ports;
    int         num_nodes;

    int         udn_index;
    int         num_outputs;

    Mif_Boolean_t  invert;

    Evt_Node_Data_t     *node_data;
    Evt_State_Data_t    *state_data;
    Evt_Msg_Data_t      *msg_data;
    /*    Evt_Statistic_t     *statistics_data;*/

    Evt_Node_t          *rhs;
    Evt_Node_t          *rhsold;

    Evt_Node_Info_t     *node_info;


    /* Allocate main substructures of data */
    /* Note that we don't free any old structures */
    /* since they are pointed to by jobs and need */
    /* to be maintained so that results from multiple */
    /* jobs are kept around like SPICE does */

    data = &(ckt->evt->data);
    CKALLOC(data->node, 1, Evt_Node_Data_t)
    CKALLOC(data->state, 1, Evt_State_Data_t)
    CKALLOC(data->msg, 1, Evt_Msg_Data_t)
    CKALLOC(data->statistics, 1, Evt_Statistic_t)

    /* Allocate node data */

    num_nodes = ckt->evt->counts.num_nodes;
    node_data = data->node;

    CKALLOC(node_data->head, num_nodes, Evt_Node_t *)
    CKALLOC(node_data->tail, num_nodes, Evt_Node_t **)
    CKALLOC(node_data->last_step, num_nodes, Evt_Node_t **)
    CKALLOC(node_data->free, num_nodes, Evt_Node_t *)
    CKALLOC(node_data->modified_index, num_nodes, int)
    CKALLOC(node_data->modified, num_nodes, Mif_Boolean_t)
    CKALLOC(node_data->rhs, num_nodes, Evt_Node_t)
    CKALLOC(node_data->rhsold, num_nodes, Evt_Node_t)
    CKALLOC(node_data->total_load, num_nodes, double)

    /* Initialize the node data */

    for(i = 0; i < num_nodes; i++) {
        node_data->tail[i] = &(node_data->head[i]);
        node_data->last_step[i] = &(node_data->head[i]);
    }

    for(i = 0; i < num_nodes; i++) {

        /* Get pointers to rhs & rhsold, the user-defined node type index, */
        /* the number of outputs on the node and the invert flag */
        rhs = &(node_data->rhs[i]);
        rhsold = &(node_data->rhsold[i]);
        node_info = ckt->evt->info.node_table[i];
        udn_index = node_info->udn_index;
        num_outputs = node_info->num_outputs;
        invert = node_info->invert;

        /* Initialize the elements within rhs and rhsold */
        rhs->step = 0.0;
        rhsold->step = 0.0;
        if(num_outputs > 1) {
            CKALLOC(rhs->output_value, num_outputs, void *)
            CKALLOC(rhsold->output_value, num_outputs, void *)
            for(j = 0; j < num_outputs; j++) {
                g_evt_udn_info[udn_index]->create (&(rhs->output_value[j]));
                g_evt_udn_info[udn_index]->initialize (rhs->output_value[j]);
                g_evt_udn_info[udn_index]->create (&(rhsold->output_value[j]));
                g_evt_udn_info[udn_index]->initialize (rhsold->output_value[j]);
            }
        }
        g_evt_udn_info[udn_index]->create (&(rhs->node_value));
        g_evt_udn_info[udn_index]->initialize (rhs->node_value);
        g_evt_udn_info[udn_index]->create (&(rhsold->node_value));
        g_evt_udn_info[udn_index]->initialize (rhsold->node_value);
        if(invert) {
            g_evt_udn_info[udn_index]->create (&(rhs->inverted_value));
            g_evt_udn_info[udn_index]->initialize (rhs->inverted_value);
            g_evt_udn_info[udn_index]->create (&(rhsold->inverted_value));
            g_evt_udn_info[udn_index]->initialize (rhsold->inverted_value);
        }

        /* Initialize the total load value to zero */
        node_data->total_load[i] = 0.0;
    }


    /* Allocate and initialize state data */

    num_insts = ckt->evt->counts.num_insts;
    state_data = data->state;

    CKALLOC(state_data->head, num_insts, Evt_State_t *)
    CKALLOC(state_data->tail, num_insts, Evt_State_t **)
    CKALLOC(state_data->last_step, num_insts, Evt_State_t **)
    CKALLOC(state_data->free, num_insts, Evt_State_t *)
    CKALLOC(state_data->modified_index, num_insts, int)
    CKALLOC(state_data->modified, num_insts, Mif_Boolean_t)
    CKALLOC(state_data->total_size, num_insts, int)
    CKALLOC(state_data->desc, num_insts, Evt_State_Desc_t *)

    for(i = 0; i < num_insts; i++) {
        state_data->tail[i] = &(state_data->head[i]);
        state_data->last_step[i] = &(state_data->head[i]);
    }


    /* Allocate and initialize msg data */

    num_ports = ckt->evt->counts.num_ports;
    msg_data = data->msg;

    CKALLOC(msg_data->head, num_ports, Evt_Msg_t *)
    CKALLOC(msg_data->tail, num_ports, Evt_Msg_t **)
    CKALLOC(msg_data->last_step, num_ports, Evt_Msg_t **)
    CKALLOC(msg_data->free, num_ports, Evt_Msg_t *)
    CKALLOC(msg_data->modified_index, num_ports, int)
    CKALLOC(msg_data->modified, num_ports, Mif_Boolean_t)

    for(i = 0; i < num_ports; i++) {
        msg_data->tail[i] = &(msg_data->head[i]);
        msg_data->last_step[i] = &(msg_data->head[i]);
    }

    /* Don't need to initialize statistics since they were */
    /* calloc'ed above */

    return(OK);
}




/*
EVTsetup_jobs

This function prepares the jobs data for a new simulation.
*/


static int EVTsetup_jobs(
    CKTcircuit *ckt)       /* The circuit structure */
{

    int  i;
    int  num_jobs;

    Evt_Job_t  *jobs;
    Evt_Data_t *data;


    jobs = &(ckt->evt->jobs);
    data = &(ckt->evt->data);

    /* Increment the number of jobs */
    num_jobs = ++(jobs->num_jobs);

    /* Allocate/reallocate necessary pointers */
    CKREALLOC(jobs->job_name, num_jobs, char *)
    CKREALLOC(jobs->job_plot, num_jobs, char *)
    CKREALLOC(jobs->node_data, num_jobs, Evt_Node_Data_t *)
    CKREALLOC(jobs->state_data, num_jobs, Evt_State_Data_t *)
    CKREALLOC(jobs->msg_data, num_jobs, Evt_Msg_Data_t *)
    CKREALLOC(jobs->statistics, num_jobs, Evt_Statistic_t *)

    /* Fill in the pointers, etc. for this new job */
    i = num_jobs - 1;
    jobs->job_name[i] = MIFcopy(ckt->CKTcurJob->JOBname);
    jobs->job_plot[i] = NULL; /* fill in later */
    jobs->node_data[i] = data->node;
    jobs->state_data[i] = data->state;
    jobs->msg_data[i] = data->msg;
    jobs->statistics[i] = data->statistics;

    return(OK);
}



/*
EVTsetup_load_ptrs

This function setups up the required data in the MIFinstance
structure of event-driven and hybrid instances.
*/


static int EVTsetup_load_ptrs(
    CKTcircuit *ckt)            /* The circuit structure */
{

    int         i;
    int         j;
    int         k;

    int         num_insts;
    int         num_conn;
    int         num_port;
    int         num_outputs;

    int         node_index;
    int         output_subindex;

    MIFinstance         *fast;

    Mif_Conn_Data_t     *conn;

    Mif_Port_Type_t     type;
    Mif_Port_Data_t     *port;

    Evt_Node_Data_t     *node_data;


    /* This function setups up the required data in the MIFinstance */
    /* structure of event-driven and hybrid instances */

    /* Loop through all event-driven and hybrid instances */
    num_insts = ckt->evt->counts.num_insts;
    for(i = 0; i < num_insts; i++) {

        /* Get the MIFinstance pointer */
        fast = ckt->evt->info.inst_table[i]->inst_ptr;

        /* Reset init flag, required when any run is called a second time */
        fast->initialized = FALSE;

        /* Loop through all connections */
        num_conn = fast->num_conn;
        for(j = 0; j < num_conn; j++) {

            /* Skip if connection is null */
            if(fast->conn[j]->is_null)
                continue;

            conn = fast->conn[j];

            /* Loop through all ports */
            num_port = conn->size;
            for(k = 0; k < num_port; k++) {

                /* Get port data pointer for quick access */
                port = conn->port[k];

                if(port->is_null)
                    continue;

                /* Skip if port is not digital or user-defined type */
                type = port->type;
                if((type != MIF_DIGITAL) && (type != MIF_USER_DEFINED))
                    continue;

                /* Set input.pvalue to point to rhsold.node_value or to */
                /* rhsold.inverted_value as appropriate */
                node_index = port->evt_data.node_index;
                node_data  = ckt->evt->data.node;
                if(conn->is_input) {
                    if(port->invert) {
                        port->input.pvalue = node_data->rhsold[node_index].
                                             inverted_value;
                    }
                    else {
                        port->input.pvalue = node_data->rhsold[node_index].
                                             node_value;
                    }
                }

                /* Set output.pvalue to point to rhs.node_value or rhs.output_value[i] */
                /* where i is given by the output_subindex in output info */
                /* depending on whether more than one output is connected to the node. */
                /* Note that this is only for the DCOP analysis.  During a transient */
                /* analysis, new structures will be created and the pointers will */
                /* be set by EVTload */
                if(conn->is_output) {
                    num_outputs = ckt->evt->info.node_table[node_index]->num_outputs;
                    if(num_outputs <= 1) {
                        port->output.pvalue = node_data->rhs[node_index].
                                              node_value;
                    }
                    else {
                        output_subindex = port->evt_data.output_subindex;
                        port->output.pvalue = node_data->rhs[node_index].
                                              output_value[output_subindex];
                    }
                }

            } /* end for number of ports */
        } /* end for number of connections */
    } /* end for number of insts */

    return(OK);
}

/* get the analog plot name and store it into the current event job */
int EVTsetup_plot(CKTcircuit* ckt, char *plotname) {
    if (ckt->evt->counts.num_insts == 0)
        return(OK);
    
    Evt_Job_t* jobs = &(ckt->evt->jobs);
    if (jobs) {
        jobs->job_plot[jobs->num_jobs - 1] = copy(plotname);
        jobs->cur_job = jobs->num_jobs - 1;
        return OK;
    }
    return 1;
}

/* If command 'setplot' is called, we switch to the corresponding event data.
   Their pointers have been stored in the jobs structure. The circuit must
   be still available!
*/
int EVTswitch_plot(CKTcircuit* ckt, const char* plottypename) {
    int i;
    bool found = FALSE;

    Evt_Job_t* jobs;
    Evt_Data_t* data;

    if (!ckt)
        return (E_NOTFOUND);

    if (ckt->evt->counts.num_insts == 0)
        return(E_NOTFOUND);

    jobs = &(ckt->evt->jobs);
    data = &(ckt->evt->data);

    if (jobs) {
        /* check for the job with current plot type name , e.g. tran2 */
        for (i = 0; i < jobs->num_jobs; i++) {
            if (jobs->job_plot[i] && eq(jobs->job_plot[i], plottypename)) {
                found = TRUE;
                jobs->cur_job = i;
                break;
            }
        }
        if (found) {
            data->node = jobs->node_data[i];
            data->state = jobs->state_data[i];
            data->msg = jobs->msg_data[i];
            data->statistics = jobs->statistics[i];
            return OK;
        }
    }
    return 1;
}


/* ======================================================================================================================
 * card #551: restartable XSPICE state for `.pss method=shoot`.
 *
 * A shooting period run is a fresh DCtran(restart=1) from t = 0 with the native C/L states put back.  For XSPICE circuits that
 * lost (a) the analog code-model state (CKTstates entries allocated through cm_analog_alloc / cm_analog_integrate / the old
 * input of every port), which was left at whatever the last run (trial, perturbation) ended with, and (b) the whole event
 * system, which is only (re)initialised by EVTsetup() once per job: a second DCtran() restart crashed in the first hybrid
 * model (output pointers still aimed at the previous run's event structs) and, had it not, would have run on stale queues.
 *
 * Here the logical state at the end of a run (analog state values, digital node values, the tail state block of every event
 * instance, the pending instance and output events with their times relative to the end of the run) is copied into an XSsnap
 * and put back before the next run.  Time: the models must see a CONTINUING time (g_xs_toff, see mif.h): their stored change
 * times / phases stay valid and the `TIME == 0.0` initialisation branches run only in the very first run.
 * ====================================================================================================================== */

typedef struct { void *nv, *inv; void **out; } XSnodeVal;

struct XSsnap {
    int        valid;
    double     toff;               /* g_xs_toff of a run that starts from this snapshot */
    int        n_an;               /* analog: CKTstate indices and their values */
    int       *an_idx;
    double    *an_val;
    int        num_nodes, num_insts;
    XSnodeVal *rhs, *rhsold;       /* [num_nodes] */
    void     **st_blk;             /* [num_insts] tail state block (or NULL) */
    int       *st_size;
    int        n_ie;               /* pending instance events */
    int       *ie_inst;
    double    *ie_t;
    int        n_oe;               /* pending output events */
    int       *oe_out;
    double    *oe_t;
    void     **oe_val;
};

static void
xs_for_each_mif_instance(CKTcircuit *ckt, void (*fn)(MIFinstance *, void *), void *ctx)
{
    int i;
    for (i = 0; i < DEVmaxnum; i++) {
        GENmodel *m;
        if (!DEVices[i] || !DEVices[i]->DEVpublic.cm_func || !ckt->CKThead[i])
            continue;
        for (m = ckt->CKThead[i]; m; m = m->GENnextModel) {
            GENinstance *g;
            for (g = m->GENinstances; g; g = g->GENnextInstance)
                fn((MIFinstance *) g, ctx);
        }
    }
}

int
XSPICEnumInstances(CKTcircuit *ckt)
{
    int i, n = 0;
    for (i = 0; i < DEVmaxnum; i++) {
        GENmodel *m;
        if (!DEVices[i] || !DEVices[i]->DEVpublic.cm_func || !ckt->CKThead[i])
            continue;
        for (m = ckt->CKThead[i]; m; m = m->GENnextModel) {
            GENinstance *g;
            for (g = m->GENinstances; g; g = g->GENnextInstance)
                n++;
        }
    }
    return n;
}

typedef struct { int n, cap; int *idx; } XSidxList;

static void
xs_idx_add(XSidxList *l, int idx)
{
    if (l->n == l->cap) {
        l->cap = l->cap ? 2 * l->cap : 64;
        l->idx = TREALLOC(int, l->idx, l->cap);
    }
    l->idx[l->n++] = idx;
}

static void
xs_collect_analog(MIFinstance *here, void *ctx)
{
    XSidxList *l = (XSidxList *) ctx;
    int i, j, k, d;
    for (k = 0; k < here->num_state; k++)
        for (d = 0; d < here->state[k].doubles; d++)
            xs_idx_add(l, here->state[k].index + d);
    if (!here->analog)
        return;
    for (i = 0; i < here->num_conn; i++) {
        if (here->conn[i]->is_null || !here->conn[i]->is_input)
            continue;
        for (j = 0; j < here->conn[i]->size; j++)
            if (!here->conn[i]->port[j]->is_null)
                xs_idx_add(l, here->conn[i]->port[j]->old_input);
    }
}

static void *
xs_udn_clone(int udn, void *src)
{
    void *p = NULL;
    g_evt_udn_info[udn]->create(&p);
    g_evt_udn_info[udn]->copy(src, p);
    return p;
}

static void
xs_udn_free(int udn, void *p)
{
    if (!p)
        return;
    if (g_evt_udn_info[udn]->dismantle)
        g_evt_udn_info[udn]->dismantle(p);
    tfree(p);
}

static void
xs_nodeval_clone(CKTcircuit *ckt, int i, const Evt_Node_t *src, XSnodeVal *dst)
{
    Evt_Node_Info_t *ni = ckt->evt->info.node_table[i];
    int j;
    dst->nv = xs_udn_clone(ni->udn_index, src->node_value);
    dst->inv = (ni->invert && src->inverted_value) ? xs_udn_clone(ni->udn_index, src->inverted_value) : NULL;
    dst->out = NULL;
    if (ni->num_outputs > 1 && src->output_value) {
        dst->out = TMALLOC(void *, ni->num_outputs);
        for (j = 0; j < ni->num_outputs; j++)
            dst->out[j] = xs_udn_clone(ni->udn_index, src->output_value[j]);
    }
}

static void
xs_nodeval_free(CKTcircuit *ckt, int i, XSnodeVal *v)
{
    Evt_Node_Info_t *ni = ckt->evt->info.node_table[i];
    int j;
    xs_udn_free(ni->udn_index, v->nv);
    xs_udn_free(ni->udn_index, v->inv);
    if (v->out) {
        for (j = 0; j < ni->num_outputs; j++)
            xs_udn_free(ni->udn_index, v->out[j]);
        tfree(v->out);
    }
    v->nv = v->inv = NULL;
    v->out = NULL;
}

static void
xs_nodeval_put(CKTcircuit *ckt, int i, const XSnodeVal *src, Evt_Node_t *dst)
{
    Evt_Node_Info_t *ni = ckt->evt->info.node_table[i];
    int j;
    g_evt_udn_info[ni->udn_index]->copy(src->nv, dst->node_value);
    if (src->inv && dst->inverted_value)
        g_evt_udn_info[ni->udn_index]->copy(src->inv, dst->inverted_value);
    if (src->out && dst->output_value)
        for (j = 0; j < ni->num_outputs; j++)
            g_evt_udn_info[ni->udn_index]->copy(src->out[j], dst->output_value[j]);
}

XSsnap *
XSsnapNew(void)
{
    return TMALLOC(XSsnap, 1);
}

static void
xs_snap_clear(CKTcircuit *ckt, XSsnap *s)
{
    int i;
    if (s->rhs)
        for (i = 0; i < s->num_nodes; i++)
            xs_nodeval_free(ckt, i, &s->rhs[i]);
    if (s->rhsold)
        for (i = 0; i < s->num_nodes; i++)
            xs_nodeval_free(ckt, i, &s->rhsold[i]);
    if (s->st_blk)
        for (i = 0; i < s->num_insts; i++)
            if (s->st_blk[i])
                tfree(s->st_blk[i]);
    if (s->oe_val)
        for (i = 0; i < s->n_oe; i++) {
            int out = s->oe_out[i];
            int node = ckt->evt->info.output_table[out]->node_index;
            xs_udn_free(ckt->evt->info.node_table[node]->udn_index, s->oe_val[i]);
        }
    tfree(s->an_idx); tfree(s->an_val); tfree(s->rhs); tfree(s->rhsold); tfree(s->st_blk); tfree(s->st_size);
    tfree(s->ie_inst); tfree(s->ie_t); tfree(s->oe_out); tfree(s->oe_t); tfree(s->oe_val);
    memset(s, 0, sizeof *s);
}

void
XSsnapFree(CKTcircuit *ckt, XSsnap *s)
{
    if (!s)
        return;
    xs_snap_clear(ckt, s);
    tfree(s);
}

int
XSsnapValid(const XSsnap *s)
{
    return s && s->valid;
}

/* Take the logical XSPICE state at the end of a run of length L (CKTtime == L).  base_toff is the g_xs_toff the run used. */
void
XSsnapTake(CKTcircuit *ckt, XSsnap *s, double L, double base_toff)
{
    XSidxList l = { 0, 0, NULL };
    int i, k;

    xs_snap_clear(ckt, s);
    s->toff = base_toff + L;

    xs_for_each_mif_instance(ckt, xs_collect_analog, &l);
    if (l.n && ckt->CKTstate0) {
        s->n_an = l.n;
        s->an_idx = l.idx;
        s->an_val = TMALLOC(double, l.n);
        for (k = 0; k < l.n; k++)
            s->an_val[k] = ckt->CKTstate0[l.idx[k]];
    } else {
        tfree(l.idx);
    }

    if (ckt->evt && ckt->evt->counts.num_insts > 0 && ckt->evt->data.node) {
        Evt_Node_Data_t *nd = ckt->evt->data.node;
        Evt_State_Data_t *sd = ckt->evt->data.state;
        Evt_Inst_Queue_t *iq = &ckt->evt->queue.inst;
        Evt_Output_Queue_t *oq = &ckt->evt->queue.output;
        int ni = ckt->evt->counts.num_nodes, ns = ckt->evt->counts.num_insts, no = ckt->evt->counts.num_outputs;
        int cap;

        s->num_nodes = ni;
        s->num_insts = ns;
        s->rhs = TMALLOC(XSnodeVal, ni);
        s->rhsold = TMALLOC(XSnodeVal, ni);
        for (i = 0; i < ni; i++) {
            xs_nodeval_clone(ckt, i, &nd->rhs[i], &s->rhs[i]);
            xs_nodeval_clone(ckt, i, &nd->rhsold[i], &s->rhsold[i]);
        }
        s->st_blk = TMALLOC(void *, ns);
        s->st_size = TMALLOC(int, ns);
        for (i = 0; i < ns; i++) {
            Evt_State_t *st = sd->tail[i] ? *(sd->tail[i]) : NULL;
            s->st_size[i] = sd->total_size[i];
            if (st && st->block && sd->total_size[i] > 0) {
                s->st_blk[i] = tmalloc((size_t) sd->total_size[i]);
                memcpy(s->st_blk[i], st->block, (size_t) sd->total_size[i]);
            }
        }
        cap = 0;
        for (i = 0; i < ns; i++) {
            Evt_Inst_Event_t *ev;
            if (!iq->current[i])
                continue;
            for (ev = *(iq->current[i]); ev; ev = ev->next) {
                if (s->n_ie == cap) {
                    cap = cap ? 2 * cap : 16;
                    s->ie_inst = TREALLOC(int, s->ie_inst, cap);
                    s->ie_t = TREALLOC(double, s->ie_t, cap);
                }
                s->ie_inst[s->n_ie] = i;
                s->ie_t[s->n_ie] = ev->event_time - L;
                s->n_ie++;
            }
        }
        cap = 0;
        for (i = 0; i < no; i++) {
            Evt_Output_Event_t *ev;
            int node = ckt->evt->info.output_table[i]->node_index;
            int udn = ckt->evt->info.node_table[node]->udn_index;
            if (!oq->current[i])
                continue;
            for (ev = *(oq->current[i]); ev; ev = ev->next) {
                if (ev->removed)
                    continue;
                if (s->n_oe == cap) {
                    cap = cap ? 2 * cap : 16;
                    s->oe_out = TREALLOC(int, s->oe_out, cap);
                    s->oe_t = TREALLOC(double, s->oe_t, cap);
                    s->oe_val = TREALLOC(void *, s->oe_val, cap);
                }
                s->oe_out[s->n_oe] = i;
                s->oe_t[s->n_oe] = ev->event_time - L;
                s->oe_val[s->n_oe] = xs_udn_clone(udn, ev->value);
                s->n_oe++;
            }
        }
    }
    s->valid = 1;
    if (getenv("PSS_XS_TRACE")) {
        fprintf(stderr, "XS-TRACE take: toff_next=%.9g n_an=%d", s->toff, s->n_an);
        for (k = 0; k < s->n_an && k < 8; k++)
            fprintf(stderr, " [%d]=%.6g", s->an_idx[k], s->an_val[k]);
        fprintf(stderr, " digital: nodes=%d insts=%d events inst=%d out=%d", s->num_nodes, s->num_insts, s->n_ie, s->n_oe);
        for (k = 0; k < s->num_nodes && k < 12; k++) {
            double rv = -1.0;
            Evt_Node_Info_t *ni = ckt->evt->info.node_table[k];
            if (g_evt_udn_info[ni->udn_index]->plot_val)
                g_evt_udn_info[ni->udn_index]->plot_val(s->rhsold[k].nv, "state", &rv);
            fprintf(stderr, " %s=%g", ni->name, rv);
        }
        for (k = 0; k < s->n_ie; k++) fprintf(stderr, " ie[%d]@%.6e", s->ie_inst[k], s->ie_t[k]);
        for (k = 0; k < s->n_oe; k++) fprintf(stderr, " oe[%d]@%.6e", s->oe_out[k], s->oe_t[k]);
        fprintf(stderr, "\n");
    }
}

/* Put the snapshot back before a restarted run (DCtran restart = 1, CKTtime = 0). */
void
XSsnapPut(CKTcircuit *ckt, const XSsnap *s, int values)
{
    int i, k;
    if (!s || !s->valid)
        return;
    /* values == 0 (PSS_XS_NO_STATE=1, red proof only): restart the event framework but carry no state (the pre-#551 restart) */
    g_xs_toff = values ? s->toff : 0.0;
    if (values && s->n_an && ckt->CKTstate0) {
        for (k = 0; k < s->n_an; k++) {
            ckt->CKTstate0[s->an_idx[k]] = s->an_val[k];
            if (ckt->CKTstate1)
                ckt->CKTstate1[s->an_idx[k]] = s->an_val[k];
        }
    }
    if (ckt->evt && ckt->evt->counts.num_insts > 0 && ckt->evt->data.node && s->rhs) {
        Evt_Node_Data_t *nd = ckt->evt->data.node;
        Evt_State_Data_t *sd = ckt->evt->data.state;
        int ni = ckt->evt->counts.num_nodes, ns = ckt->evt->counts.num_insts;
        char *was_init = TMALLOC(char, ns);

        for (i = 0; i < ns; i++)
            was_init[i] = ckt->evt->info.inst_table[i]->inst_ptr->initialized ? 1 : 0;

        EVTsetup_queues(ckt);
        EVTsetup_load_ptrs(ckt);
        for (i = 0; i < ns; i++)
            ckt->evt->info.inst_table[i]->inst_ptr->initialized = was_init[i] ? MIF_TRUE : MIF_FALSE;
        tfree(was_init);

        for (i = 0; i < ni; i++) {
            /* history list -> free list; EVTop_save() reseeds it from rhsold */
            if (nd->head[i]) {
                Evt_Node_t *last = nd->head[i];
                while (last->next)
                    last = last->next;
                last->next = nd->free[i];
                nd->free[i] = nd->head[i];
                nd->head[i] = NULL;
            }
            nd->tail[i] = &(nd->head[i]);
            nd->last_step[i] = &(nd->head[i]);
            nd->modified[i] = MIF_FALSE;
            if (values) {
                xs_nodeval_put(ckt, i, &s->rhs[i], &nd->rhs[i]);
                xs_nodeval_put(ckt, i, &s->rhsold[i], &nd->rhsold[i]);
            }
        }
        nd->num_modified = 0;

        for (i = 0; i < ns; i++) {
            Evt_State_t *st = sd->tail[i] ? *(sd->tail[i]) : NULL;
            if (!st || !st->block || !s->st_blk[i] || s->st_size[i] != sd->total_size[i])
                continue;
            if (values)
                memcpy(st->block, s->st_blk[i], (size_t) s->st_size[i]);
            if (sd->head[i] != st) {              /* history still linked: older steps go to the free list, as EVTaccept() does */
                if (st->prev) {
                    st->prev->next = sd->free[i];
                    sd->free[i] = sd->head[i];
                }
                sd->head[i] = st;
            }
            st->prev = NULL;                      /* (after EVTaccept() prev dangles into the free list; head == st already) */
            st->next = NULL;
            st->step = 0.0;
            sd->tail[i] = &(sd->head[i]);
            sd->last_step[i] = &(sd->head[i]);
            sd->modified[i] = MIF_FALSE;
        }
        sd->num_modified = 0;

        for (k = 0; values && k < s->n_ie; k++)
            EVTqueue_inst(ckt, s->ie_inst[k], 0.0, s->ie_t[k]);
        for (k = 0; values && k < s->n_oe; k++) {
            int out = s->oe_out[k];
            int node = ckt->evt->info.output_table[out]->node_index;
            int udn = ckt->evt->info.node_table[node]->udn_index;
            Evt_Output_Event_t *ev = TMALLOC(Evt_Output_Event_t, 1);
            g_evt_udn_info[udn]->create(&ev->value);
            g_evt_udn_info[udn]->copy(s->oe_val[k], ev->value);
            g_xs_dig_cur_tag = (g_xs_dig_on && g_xs_dig_carry_tag && out < g_xs_dig_ncarry) ? g_xs_dig_carry_tag[out] : 0;   /* card #714 */
            EVTqueue_output(ckt, out, udn, ev, 0.0, s->oe_t[k]);
            g_xs_dig_cur_tag = 0;
        }
        g_mif_info.circuit.evt_step = 0.0;
        EVTop_save(ckt, MIF_FALSE, 0.0);   /* the (emptied) node history is seeded from rhsold exactly as after EVTop() */
    }
}

/* The value a digital node has once the output events that fall due within XS_EDGE_EPS of the snapshot instant have been applied.  A clock edge
 * that lands on the period boundary is processed at the end of one run or at the start of the next depending on rounding of the continuing model
 * time (10^-20 s); both are the same orbit, so the periodicity comparison must not tell them apart (card #712). */
#define XS_EDGE_EPS 1e-15
static const void *
xs_snap_node_effective(CKTcircuit *ckt, const XSsnap *s, int node)
{
    const void *v = s->rhsold[node].nv;
    double tbest = -1e300;
    int k;
    for (k = 0; k < s->n_oe; k++) {
        int out = s->oe_out[k];
        if (ckt->evt->info.output_table[out]->node_index != node)
            continue;
        if (s->oe_t[k] <= XS_EDGE_EPS && s->oe_t[k] >= tbest) {
            tbest = s->oe_t[k];
            v = s->oe_val[k];
        }
    }
    return v;
}

/* Number of digital nodes whose resolved value differs between two snapshots (hidden-state periodicity check). */
int
XSsnapDigitalDiff(CKTcircuit *ckt, const XSsnap *a, const XSsnap *b)
{
    int i, n = 0;
    if (!a || !b || !a->valid || !b->valid || !a->rhsold || !b->rhsold || a->num_nodes != b->num_nodes)
        return 0;
    for (i = 0; i < a->num_nodes; i++) {
        Mif_Boolean_t eq = MIF_TRUE;
        const void *va = xs_snap_node_effective(ckt, a, i), *vb = xs_snap_node_effective(ckt, b, i);
        g_evt_udn_info[ckt->evt->info.node_table[i]->udn_index]->compare((void *) va, (void *) vb, &eq);
        if (!eq)
            n++;
    }
    return n;
}

/* Largest |a-b| over the analog state values of two snapshots, relative to max(|a|,|b|,floor) (diagnostic only). */
double
XSsnapAnalogDiff(const XSsnap *a, const XSsnap *b)
{
    int k;
    double m = 0.0;
    if (!a || !b || !a->n_an || a->n_an != b->n_an)
        return 0.0;
    for (k = 0; k < a->n_an; k++) {
        double x = a->an_val[k], y = b->an_val[k];
        double sc = fmax(fmax(fabs(x), fabs(y)), 1e-9);
        double d = fabs(x - y) / sc;
        if (d > m)
            m = d;
    }
    return m;
}

/* ------------------------------------------------------------------------------------------------------------------------------- */
/* card #714: EXACT LOCATION OF adc_bridge THRESHOLD CROSSINGS.
 *
 * The stock adc_bridge compares its input with in_low / in_high at every analog time point and, when the digital value it would give
 * differs from its state, schedules itself at that time point: the digital edge is quantised to the adaptive time grid (up to one step
 * late).  In a closed loop that makes the shooting map x0 -> x(T) piecewise constant in the edge time; the finite-difference monodromy
 * then sees 0 or a whole step, and the Newton iteration of `.pss` stalls (card #712: residual 0.75 after 150 iterations on a comparator
 * loop) or "converges" with a wrong Floquet radius.  Here a converged trial step in which the input of a single-bit adc_bridge crossed
 * the level its state changes at is rejected and redone to end at the linearly interpolated crossing (plus half the tolerance), until
 * the overshoot is below XS_ADC_TTOL (default 1e-16 s), exactly like the ideal-switch location of pssshoot.c (card #498).  The digital
 * edge then is a smooth function of the analog state.
 *
 * Active in the shooting and td recording runs (XSperiodicRun; PSS_ADC_LOCATE=0 turns it off) and in any transient with
 * XS_ADC_LOCATE=1 (so a settled .tran / FRA reference can use the same edge times).  Not done: multi-bit adc_bridge (one input, several
 * outputs), the input types vnam; their edges stay on the grid. */
int XSadcLocateActive = 0;
int XSperiodicRun = 0;

typedef struct {
    MIFinstance *inst;
    int inst_index, port, type, pos, neg, ibr;
    double lo, hi;
    double t_acc, v_acc;   /* input at the last accepted time point */
} XSadcPort;
static XSadcPort *g_adc = NULL;
static int g_adc_n = -1;
static CKTcircuit *g_adc_ckt = NULL;
static double g_adc_tacc = -1.0;
static int g_adc_iter = 0;

int
XSadcLocateWanted(void)
{
    const char *e = getenv("XS_ADC_LOCATE"), *p = getenv("PSS_ADC_LOCATE");
    if (e && atoi(e) != 0)
        return 1;
    if (XSperiodicRun && !(p && atoi(p) == 0))
        return 1;
    return 0;
}

static double
xs_param_real(MIFinstance *inst, const char *name, double dflt)
{
    int type = MIFmodPtr(inst)->MIFmodType, k;
    IFdevice *pub = &DEVices[type]->DEVpublic;
    for (k = 0; k < pub->num_param && k < inst->num_param; k++)
        if (pub->param[k].name && !strcmp(pub->param[k].name, name) && inst->param[k] && inst->param[k]->element)
            return inst->param[k]->element[0].rvalue;
    return dflt;
}

static void
xs_adc_collect(CKTcircuit *ckt)
{
    int i, cap = 0;
    tfree(g_adc);
    g_adc = NULL;
    g_adc_n = 0;
    g_adc_ckt = ckt;
    if (!ckt->evt)
        return;
    for (i = 0; i < ckt->evt->counts.num_insts; i++) {
        MIFinstance *inst = ckt->evt->info.inst_table[i]->inst_ptr;
        const char *nm = DEVices[MIFmodPtr(inst)->MIFmodType]->DEVpublic.name;
        Mif_Conn_Data_t *cin, *cout;
        int j;
        if (!nm || strcmp(nm, "adc_bridge") || inst->num_conn < 2)
            continue;
        cin = inst->conn[0];
        cout = inst->conn[1];
        if (cin->is_null || cout->is_null || cin->size != cout->size)
            continue;   /* multi-bit (one input, several outputs): not located */
        for (j = 0; j < cin->size; j++) {
            Mif_Port_Data_t *p = cin->port[j];
            XSadcPort *a;
            if (p->is_null || !(p->type == MIF_VOLTAGE || p->type == MIF_DIFF_VOLTAGE || p->type == MIF_CURRENT || p->type == MIF_DIFF_CURRENT))
                continue;
            if (g_adc_n == cap) {
                cap = cap ? 2 * cap : 8;
                g_adc = TREALLOC(XSadcPort, g_adc, cap);
            }
            a = &g_adc[g_adc_n++];
            a->inst = inst;
            a->inst_index = inst->inst_index;
            a->port = j;
            a->type = p->type;
            a->pos = p->smp_data.pos_node;
            a->neg = p->smp_data.neg_node;
            a->ibr = p->smp_data.ibranch;
            a->lo = xs_param_real(inst, "in_low", 0.1);
            a->hi = xs_param_real(inst, "in_high", 0.9);
            a->t_acc = -1e300;
            a->v_acc = 0.0;
        }
    }
}

static double
xs_adc_value(const XSadcPort *a, const double *x)
{
    if (a->type == MIF_CURRENT || a->type == MIF_DIFF_CURRENT)
        return x[a->ibr];
    return x[a->pos] - x[a->neg];
}

/* the digital state the adc_bridge holds for output j (its event state block: one Digital_State_t per output) */
static int
xs_adc_state(CKTcircuit *ckt, const XSadcPort *a)
{
    Evt_State_Data_t *sd = ckt->evt->data.state;
    Evt_State_t *st;
    if (a->inst_index < 0 || !sd->tail || !sd->tail[a->inst_index])
        return -1;
    st = *(sd->tail[a->inst_index]);
    if (!st || !st->block)
        return -1;
    return (int) ((Digital_State_t *) st->block)[a->port];
}

/* digital value for input v and the level at which state o changes to n (-1: no change / not located) */
static int
xs_adc_level(const XSadcPort *a, int o, double v, double *L)
{
    int n;
    if (o < 0 || o > 2)
        return -1;
    if (a->lo <= a->hi) {
        n = (v >= a->hi) ? 1 : (v <= a->lo) ? 0 : 2;
        if (n == o) return -1;
        *L = (n == 1) ? a->hi : (n == 0) ? a->lo : (o == 0 ? a->lo : a->hi);
    } else {   /* Schmitt trigger: the band keeps the state */
        n = (v >= a->lo) ? 1 : (v <= a->hi) ? 0 : o;
        if (n == o) return -1;
        *L = (n == 1) ? a->lo : a->hi;
    }
    return n;
}

void
XSadcAccept(void *vckt)
{
    CKTcircuit *ckt = (CKTcircuit *) vckt;
    int k;
    if (!XSadcLocateActive)
        return;
    if (g_adc_n < 0 || g_adc_ckt != ckt)
        xs_adc_collect(ckt);
    for (k = 0; k < g_adc_n; k++) {
        g_adc[k].t_acc = ckt->CKTtime;
        g_adc[k].v_acc = xs_adc_value(&g_adc[k], ckt->CKTrhs);
    }
}

/* 1 = reject this converged step and retry with *newdelta (to the estimated crossing); 0 = carry on */
int
XSadcLocate(void *vckt, double *newdelta)
{
    CKTcircuit *ckt = (CKTcircuit *) vckt;
    static double ttol = -1.0;
    double t = ckt->CKTtime, tacc = ckt->CKTtime - ckt->CKTdelta, tbest = 1e300;
    int k, rej = 0;
    if (ttol < 0.0)
        ttol = getenv("XS_ADC_TTOL") ? atof(getenv("XS_ADC_TTOL")) : 1e-16;
    if (g_adc_n < 0 || g_adc_ckt != ckt)
        xs_adc_collect(ckt);
    if (g_adc_n == 0)
        return 0;
    if (tacc != g_adc_tacc) {
        g_adc_tacc = tacc;
        g_adc_iter = 0;
    }
    for (k = 0; k < g_adc_n; k++) {
        XSadcPort *a = &g_adc[k];
        double v, L, tc;
        int o;
        if (fabs(a->t_acc - tacc) > 1e-12 * fabs(tacc) + 1e-30)
            continue;   /* no input value at the last accepted point */
        o = xs_adc_state(ckt, a);
        v = xs_adc_value(a, ckt->CKTrhs);   /* NIiter leaves the converged x_k in CKTrhs (see pacshoot.c) */
        if (xs_adc_level(a, o, v, &L) < 0)
            continue;
        if ((a->v_acc - L) * (v - L) > 0.0 || v == a->v_acc)
            continue;   /* not a crossing of this step */
        tc = tacc + (t - tacc) * (L - a->v_acc) / (v - a->v_acc);
        if (t - tc > ttol && tc < tbest)
            tbest = tc;
    }
    if (tbest < 1e300 && g_adc_iter < 60) {
        double nd = tbest + 0.5 * ttol - tacc;
        if (nd < 2.0 * ckt->CKTdelmin)
            nd = 2.0 * ckt->CKTdelmin;
        if (nd < ckt->CKTdelta) {
            *newdelta = nd;
            g_adc_iter++;
            rej = 1;
        }
    }
    return rej;
}

/* card #714: the first non-zero timing tag among the output events still queued (0: none) */
int
XSdigPendingTag(void *vckt, int (*allowed)(int tag))
{
    CKTcircuit *ckt = (CKTcircuit *) vckt;
    Evt_Output_Queue_t *oq;
    int i;
    if (!ckt->evt)
        return 0;
    oq = &(ckt->evt->queue.output);
    for (i = 0; i < ckt->evt->counts.num_outputs; i++) {
        Evt_Output_Event_t *ev;
        if (!oq->current || !oq->current[i])
            continue;
        for (ev = *(oq->current[i]); ev; ev = ev->next)
            if (!ev->removed && ev->tag != 0 && !(allowed && allowed(ev->tag)))
                return ev->tag;
    }
    return 0;
}
