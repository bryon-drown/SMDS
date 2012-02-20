/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    SEDH Tabulator Helper Class
 *    $Id: Tabulator.c,v 1.1.1.1 2006/05/18 21:12:43 mculbert Exp $
 *
 *******************************************/
 
#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL smds_array_symbol
#include <numpy/arrayobject.h>
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE { Py_INCREF(Py_None); return Py_None; }
#endif

typedef long Int32;
typedef short Int16;

typedef struct _smds_SEDH_Tabulator_extent extent;
struct _smds_SEDH_Tabulator_extent { 
  int num; 
  Int32 *data; 
  extent *next;
};

typedef struct {
    PyObject_HEAD
    double binWidthIn, binWidthHist;
    int threshold;			// per bin
    int reBin;
    extent *head;
    int in_event;
    Int32 event_start;
    Int16 reBinPool;        /* If the data stop before we've finished an */
    int reBinRemaining;     /*   aggregate bin, store the current photon
                                 bin counts here */
} smds_SEDH_Tabulator;

// takes three args: binWidthIn (ms), threshold (phot/ms), binWidthHist (ms)
static int smds_SEDH_Tabulator_new(PyObject *self, 
					PyObject *args, PyObject *kwds);
static void smds_SEDH_Tabulator_free(smds_SEDH_Tabulator *);
// takes an smds.Task.RTR
static PyObject * smds_SEDH_Tabulator_tabulate(PyObject *self, PyObject *arg);
static PyObject * smds_SEDH_Tabulator_fetch(PyObject *self);

static PyTypeObject smds_TabulatorType = {
  PyObject_HEAD_INIT(NULL)
  0,				// ob_size
  "SEDH Tabulator",		// tp_name
  sizeof(smds_SEDH_Tabulator),	// tp_basicsize
  0,				// tp_itemsize
  (destructor)smds_SEDH_Tabulator_free,	// tp_dealloc
  0,				// tp_print
  0,				// tp_getattr
  0,				// tp_setattr
  0,				// tp_compare
  0,				// tp_repr
  0,				// tp_as_number
  0,				// tp_as_sequence
  0,				// tp_as_mapping
  0,				// tp_hash
  0,				// tp_call
  0,				// tp_str
  0,				// tp_getattro
  0,				// tp_setattro
  0,				// tp_as_buffer
  Py_TPFLAGS_DEFAULT,		// tp_flags
  "Multi-tau Auto-Correlator",	// tp_doc
};

static PyObject *SEDHType = NULL;
static PyObject *RTRType = NULL;

static PyMethodDef smds_SEDH_Tabulator_methods[] = {
  {"tabulate", (PyCFunction)smds_SEDH_Tabulator_tabulate, METH_VARARGS,
   "Append more data to the SEDH.\n  "
   "Argument is of type RTR."},
  {"fetch", (PyCFunction)smds_SEDH_Tabulator_fetch, METH_NOARGS,
   "Returns the smds.Analysis.SEDH corresponding the the tabulated data."},
  {NULL}	// Sentinel
};

static PyMethodDef module_methods[] = {
  {NULL}	// Sentinel
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initTabulator(void)
{
  PyObject *m;
  smds_TabulatorType.tp_new = PyType_GenericNew;
  smds_TabulatorType.tp_init = smds_SEDH_Tabulator_new;
  smds_TabulatorType.tp_methods = smds_SEDH_Tabulator_methods;
  if (PyType_Ready(&smds_TabulatorType) < 0)
    return;
  m = Py_InitModule3("Tabulator", module_methods,
  			"A Single-Event Duration Histogram Tabulator");
  Py_INCREF(&smds_TabulatorType);
  PyModule_AddObject(m, "Tabulator", (PyObject*)&smds_TabulatorType);
  import_array();

  m = PyImport_ImportModule("smds.Analysis");
  if (m)
  {
    SEDHType = PyObject_GetAttrString(m, "SEDH");
    Py_DECREF(m);
  }
  m = PyImport_ImportModule("smds.Task");
  if (m)
  {
    RTRType = PyObject_GetAttrString(m, "RTR");
    Py_DECREF(m);
  }
}


// Num must be > 0
static extent * smds_SEDH_Tabulator_extent_new(int num)
{
  extent *r;
  r = (extent*)malloc(sizeof(extent));
  r->num = num;
  r->data = (Int32*)malloc(sizeof(Int32)*num);
  for (num--;  num >= 0; num--)
    r->data[num] = 0;
  r->next = NULL;
  return r;
}

static void smds_SEDH_Tabulator_extent_free(extent *p)
{
  if (!p) return;
  smds_SEDH_Tabulator_extent_free(p->next);
  free(p->data);
  free(p);
}

/***
   Usual constructor:
     Begin tabulating an SEDH of the given threshold 
     from data with the given bin Width, reBinning as needed.
   Arguements: bw (ms), thr (phot/ms), bwH (ms)
***/
int smds_SEDH_Tabulator_new(PyObject *self, PyObject *args, PyObject *kwds)
{
  static  char *kwlist[] = { "binWidthData", "threshold", "binWidthHist",
  				NULL };
  smds_SEDH_Tabulator *t = (smds_SEDH_Tabulator*)self;
  double bw, thr, bwH;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddd", kwlist, 
  					&bw, &thr, &bwH))
    return -1;

  smds_SEDH_Tabulator_extent_free(t->head);
  t->head = smds_SEDH_Tabulator_extent_new(500);
  if (!t->head)
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    return -1;
  }
  t->binWidthIn = bw;
  if (bwH < bw) bwH = bw;		// Cannot increase resolution
  t->reBin = (int)(bwH/bw+0.5);
  t->binWidthHist = bw * (double)(t->reBin);	// Account for rounding
  t->threshold = (int)(thr*t->binWidthHist+0.5);
  t->in_event = 0;
  t->event_start = 0;
  t->reBinPool = 0;
  t->reBinRemaining = t->reBin;
  return 0;
}

/***
   Destructor
***/
void smds_SEDH_Tabulator_free(smds_SEDH_Tabulator *t)
{
  if (!t) return;
  smds_SEDH_Tabulator_extent_free(t->head);
  t->ob_type->tp_free((PyObject*)t);
}

/***
   Tabulate the next set of data.
     Takes an smds.Task.RTR
***/
PyObject * smds_SEDH_Tabulator_tabulate(PyObject *self, PyObject *args)
{
  smds_SEDH_Tabulator *t = (smds_SEDH_Tabulator*)self;
  PyObject *o;
  PyArrayObject *a;
  int i, pos = 0;		// _i_ is input pos, _pos_ is rebinned pos
  int numBins;
  Int16 *data;
  extent *p;

  if (!PyArg_ParseTuple(args, "O!", &PyInstance_Type, &o))  return NULL;
  // o is borrowed
  {
    PyObject *c = PyObject_GetAttrString(o, "__class__");
    if (!c) return NULL;
    i = PyObject_Compare(c, RTRType);
    Py_DECREF(c);
    if (PyErr_Occurred()) return NULL;
    if (i != 0)
    {
      PyErr_SetString(PyExc_TypeError, 
      		"Argument must be an instance of smds.Task.RTR.");
      return NULL;
    }
  }
  o = PyObject_GetAttrString(o, "data");
  // o is owned
  a = (PyArrayObject*)PyArray_ContiguousFromObject(o, PyArray_SHORT, 1, 1);
  if (!a)
  {
    Py_DECREF(o);
    return NULL;
  }

  numBins = a->dimensions[0];
  data = (Int16*)a->data;

  for (i=0; i < numBins; i++)
  {
    t->reBinPool += data[i];
    if (--t->reBinRemaining) continue;
    
    if (!t->in_event && t->reBinPool > t->threshold)   // Start event
    {
      t->in_event = 1;
      t->event_start = pos;
    }
    
    if (t->in_event && t->reBinPool <= t->threshold)  // End event
    {
      t->in_event = 0;
      t->event_start = pos - t->event_start;		// Event length
      p = t->head;
      while (t->event_start >= p->num)
      {
        if (p->next) p = p->next;
        else p->next = smds_SEDH_Tabulator_extent_new(500);
        t->event_start -= p->num;
      }
      p->data[t->event_start]++;
    }
    
    t->reBinPool = 0;				// Increment rebin pos
    t->reBinRemaining = t->reBin;
    pos++;
  }
  
  if (t->in_event) t->event_start -= pos;
  				// A neg starting pos for next add.

  Py_DECREF(o);
  Py_RETURN_NONE;
}

PyObject * smds_SEDH_Tabulator_fetch(PyObject *self)
{
  smds_SEDH_Tabulator *t = (smds_SEDH_Tabulator*)self;
  PyObject *o = NULL, *args;
  PyArrayObject *a;
  Int32 *data;
  int i, pos;
  extent *e;
  
  args = Py_BuildValue("(dd)", t->binWidthHist, 
  				(double)t->threshold/t->binWidthHist);
  if (!args) return NULL;
  if (SEDHType) o = PyInstance_New(SEDHType, args, NULL);
  else PyErr_SetString(PyExc_SystemError, "Unable to create SEDH object.");
  if (!o) return NULL;

  i = 0;
  for(e = t->head; e; e = e->next)
    i += e->num;
  if (t->in_event && -t->event_start > i) i = -t->event_start;
  		// Insure enough space for last event
  a = (PyArrayObject*)PyArray_FromDims(1, &i, PyArray_LONG);
  if (!a) { Py_DECREF(o);  return NULL; }
  if (PyObject_SetAttrString(o, "data", (PyObject*)a) == -1)
  {
    Py_DECREF(o);
    Py_DECREF(a);
    return NULL;
  }
  Py_DECREF(a);			// transfered reference to the SEDH object
  data = (Int32*)a->data;

  pos = 0;
  for (e = t->head; e; e = e->next)
    for(i=0; i < e->num; i++, pos++)
      data[pos] = e->data[i];
  
  for (i=pos; i < a->dimensions[0]; i++)
    data[i] = 0;
  if (t->in_event) data[-t->event_start]++;

  return o;
}

