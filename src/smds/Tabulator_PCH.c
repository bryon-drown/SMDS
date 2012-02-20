/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    PCH Tabulator Helper Class
 *    $Id: Tabulator_PCH.c,v 1.1 2008/11/19 15:48:11 mculbert Exp $
 *
 *******************************************/
 
#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL smds_array_symbol
#include <numpy/arrayobject.h>
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE { Py_INCREF(Py_None); return Py_None; }
#endif

typedef unsigned long Int32;
typedef unsigned short Int16;

typedef struct _smds_PCH_Tabulator_extent extent;
struct _smds_PCH_Tabulator_extent { 
  int num; 
  Int32 *data; 
  extent *next;
};

typedef struct {
    PyObject_HEAD
    double binWidthIn, binWidthHist;
    int reBin;
    extent *head;
    Int16 reBinPool;        /* If the data stop before we've finished an */
    int reBinRemaining;     /*   aggregate bin, store the current photon
                                 bin counts here */
} smds_PCH_Tabulator;

// takes two args: binWidthIn (ms), binWidthHist (ms)
static int smds_PCH_Tabulator_new(PyObject *self, 
					PyObject *args, PyObject *kwds);
static void smds_PCH_Tabulator_free(smds_PCH_Tabulator *);
// takes an smds.Task.RTR
static PyObject * smds_PCH_Tabulator_tabulate(PyObject *self, PyObject *arg);
static PyObject * smds_PCH_Tabulator_fetch(PyObject *self);

static PyTypeObject smds_TabulatorType = {
  PyObject_HEAD_INIT(NULL)
  0,				// ob_size
  "PCH Tabulator",		// tp_name
  sizeof(smds_PCH_Tabulator),	// tp_basicsize
  0,				// tp_itemsize
  (destructor)smds_PCH_Tabulator_free,	// tp_dealloc
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
  "Photon Counting Histogram Tabulator",	// tp_doc
};

static PyObject *PCHType = NULL;
static PyObject *RTRType = NULL;

static PyMethodDef smds_PCH_Tabulator_methods[] = {
  {"tabulate", (PyCFunction)smds_PCH_Tabulator_tabulate, METH_VARARGS,
   "Append more data to the PCH.\n  "
   "Argument is of type RTR."},
  {"fetch", (PyCFunction)smds_PCH_Tabulator_fetch, METH_NOARGS,
   "Returns the smds.Analysis.PCH corresponding the the tabulated data."},
  {NULL}	// Sentinel
};

static PyMethodDef module_methods[] = {
  {NULL}	// Sentinel
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initTabulator_PCH(void)
{
  PyObject *m;
  smds_TabulatorType.tp_new = PyType_GenericNew;
  smds_TabulatorType.tp_init = smds_PCH_Tabulator_new;
  smds_TabulatorType.tp_methods = smds_PCH_Tabulator_methods;
  if (PyType_Ready(&smds_TabulatorType) < 0)
    return;
  m = Py_InitModule3("Tabulator_PCH", module_methods,
  			"A Photon Counting Histogram Tabulator");
  Py_INCREF(&smds_TabulatorType);
  PyModule_AddObject(m, "Tabulator_PCH", (PyObject*)&smds_TabulatorType);
  import_array();

  m = PyImport_ImportModule("smds.Analysis");
  if (m)
  {
    PCHType = PyObject_GetAttrString(m, "PCH");
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
static extent * smds_PCH_Tabulator_extent_new(int num)
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

static void smds_PCH_Tabulator_extent_free(extent *p)
{
  if (!p) return;
  smds_PCH_Tabulator_extent_free(p->next);
  free(p->data);
  free(p);
}

/***
   Usual constructor:
     Begin tabulating an PCH from data with the given bin Width,
     reBinning as needed.
   Arguements: bw (ms), bwH (ms)
***/
int smds_PCH_Tabulator_new(PyObject *self, PyObject *args, PyObject *kwds)
{
  static  char *kwlist[] = { "binWidthData", "binWidthHist", NULL };
  smds_PCH_Tabulator *t = (smds_PCH_Tabulator*)self;
  double bw, bwH;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "dd", kwlist, 
  					&bw, &bwH))
    return -1;

  smds_PCH_Tabulator_extent_free(t->head);
  t->head = smds_PCH_Tabulator_extent_new(500);
  if (!t->head)
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    return -1;
  }
  t->binWidthIn = bw;
  if (bwH < bw) bwH = bw;		// Cannot increase resolution
  t->reBin = (int)(bwH/bw+0.5);
  t->binWidthHist = bw * (double)(t->reBin);	// Account for rounding
  t->reBinPool = 0;
  t->reBinRemaining = t->reBin;
  return 0;
}

/***
   Destructor
***/
void smds_PCH_Tabulator_free(smds_PCH_Tabulator *t)
{
  if (!t) return;
  smds_PCH_Tabulator_extent_free(t->head);
  t->ob_type->tp_free((PyObject*)t);
}

/***
   Tabulate the next set of data.
     Takes an smds.Task.RTR
***/
PyObject * smds_PCH_Tabulator_tabulate(PyObject *self, PyObject *args)
{
  smds_PCH_Tabulator *t = (smds_PCH_Tabulator*)self;
  PyObject *o;
  PyArrayObject *a;
  int i, numBins;
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
    
    p = t->head;
    while (t->reBinPool >= p->num)
    {
      if (p->next) p = p->next;
      else p->next = smds_PCH_Tabulator_extent_new(500);
      t->reBinPool -= p->num;
    }
    p->data[t->reBinPool]++;
    
    t->reBinPool = 0;				// Increment rebin pos
    t->reBinRemaining = t->reBin;
  }
  
  Py_DECREF(o);
  Py_RETURN_NONE;
}

PyObject * smds_PCH_Tabulator_fetch(PyObject *self)
{
  smds_PCH_Tabulator *t = (smds_PCH_Tabulator*)self;
  PyObject *o = NULL, *args;
  PyArrayObject *a;
  Int32 *data;
  int i, pos;
  extent *e;
  
  args = Py_BuildValue("(d)", t->binWidthHist);
  if (!args) return NULL;
  if (PCHType) o = PyInstance_New(PCHType, args, NULL);
  else PyErr_SetString(PyExc_SystemError, "Unable to create PCH object.");
  if (!o) return NULL;

  i = 0;
  for(e = t->head; e; e = e->next)
    i += e->num;
  a = (PyArrayObject*)PyArray_FromDims(1, &i, PyArray_LONG);
  if (!a) { Py_DECREF(o);  return NULL; }
  if (PyObject_SetAttrString(o, "data", (PyObject*)a) == -1)
  {
    Py_DECREF(o);
    Py_DECREF(a);
    return NULL;
  }
  Py_DECREF(a);			// transfered reference to the PCH object
  data = (Int32*)a->data;

  pos = 0;
  for (e = t->head; e; e = e->next)
    for(i=0; i < e->num; i++, pos++)
      data[pos] = e->data[i];
  
  return o;
}

