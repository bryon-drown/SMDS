/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Alpha-Hemolysin Channel Simulation Cores
 *    $Id: base.c,v 1.11 2009/06/16 22:36:11 mculbert Exp $
 *
 *  You must define CORE_NAME as an identifier.
 *  Optionally define CORE_DOC to a double-quoted doc-string.
 *  For variable step size, define CORE_FP.
 *  Note that these cores require that their corresponding Params
 *    and Results have the same name (CORE_NAME), but this isn't
 *    strictly necessary.  I chose this constraint for these cores
 *    to facilitate the generation of multiple cores from a single
 *    source file based on different defined macros.
 *
 *******************************************/

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL smds_array_symbol
#include <numpy/arrayobject.h>
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE { Py_INCREF(Py_None); return Py_None; }
#endif
#ifndef CORE_DOC
#define CORE_DOC "Single Molecule Diffusion Core"
#endif
#define xinitname(s) init ## s
#define initname(s) xinitname(s)
#define xstr(s) #s
#define str(s) xstr(s)

#include "random.h"
#include <math.h>

typedef unsigned long long Int64;
typedef int Int32;
typedef short Int16;
#ifdef CORE_FP
typedef double PosType;
#define ROUND 0
#else // lattice
typedef Int32 PosType;
#define ROUND 0.5
#endif

typedef struct _smds_AHL_extent extent;
struct _smds_AHL_extent { 
  int num; 
  Int32 *data; 
  extent *next;
};

typedef struct {
  PyObject_HEAD
  int ready;			// Initialization succeeded.
  double bw;
  Int64 stepsToMiss;
  Int32 stepsPerBin;
  PosType capHeight, capRadius, channelRadius, initHeight, hitPoint;
#ifdef CORE_FP
  double stepsize;
#endif
  /*** WORKING SPACE ***/
  RandomStream rs;
  /*** RESULTS ***/
  Int32 num, hits, misses;
  extent *timeToHitDist;
} smds_core;

// takes one arg: Params
static int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds);
// takes three arguments: numbins, array, pos
//    PERFORMS NO CHECKING OF PARAMETERS!
static PyObject * smds_core_run(PyObject *self, PyObject *args);
static PyObject * smds_core_getResults(PyObject *self);
static void smds_core_free(smds_core *);
static PyObject * smds_core_getParamsType(PyObject *cls);
static PyObject * smds_core_getResultsType(PyObject *cls);
static PyObject * smds_core_getName(PyObject *cls);
static PyObject * smds_core_sRnd(PyObject *self, PyObject *args);

static PyTypeObject smds_coreType = {
  PyObject_HEAD_INIT(NULL)
  0,				// ob_size
  str(CORE_NAME),			// tp_name
  sizeof(smds_core),		// tp_basicsize
  0,				// tp_itemsize
  (destructor)smds_core_free,	// tp_dealloc
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
  CORE_DOC,			// tp_doc
};

static PyMethodDef smds_core_methods[] = {
  {"run", (PyCFunction)smds_core_run, METH_VARARGS,
   "Run the simulator"},
  {"getResults", (PyCFunction)smds_core_getResults, METH_NOARGS,
   "Returns the smds.Results for the simulation."},
  {"getParamsType", (PyCFunction)smds_core_getParamsType,
    METH_NOARGS|METH_CLASS, "Returns the corresponding smds.Params class."},
  {"getResultsType", (PyCFunction)smds_core_getResultsType,
    METH_NOARGS|METH_CLASS, "Returns the corresponding smds.Results class."},
  {"getName", (PyCFunction)smds_core_getName, METH_NOARGS|METH_CLASS,
   "Returns the name of the core."},
  {NULL}	// Sentinel
};

static PyMethodDef module_methods[] = {
  {"sRnd", (PyCFunction)smds_core_sRnd, METH_VARARGS,
   "Seeds the random number generator."},
  {NULL}	// Sentinel
};

static PyObject *ParamsType = NULL;
static PyObject *ResultsType = NULL;

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initname(CORE_NAME) (void)
{
  PyObject *m;
  smds_coreType.tp_new = PyType_GenericNew;
  smds_coreType.tp_init = smds_core_init;
  smds_coreType.tp_methods = smds_core_methods;
  if (PyType_Ready(&smds_coreType) < 0)
    return;
  m = Py_InitModule3( str(CORE_NAME), module_methods,
  			CORE_DOC);
  Py_INCREF(&smds_coreType);
  PyModule_AddObject(m, str(CORE_NAME), (PyObject*)&smds_coreType);
  import_array();

//  initRnd();

  // Link core with Params and Results; Register core in smds.Cores.CoreTypes
  m = PyImport_ImportModule("smds.Results");
  if (m)
  {
    ResultsType = PyObject_GetAttrString(m, str(CORE_NAME));
    Py_DECREF(m);
  }
  
  m = PyImport_ImportModule("smds.Params");
  if (m)
  {
    ParamsType = PyObject_GetAttrString(m, str(CORE_NAME));
    Py_DECREF(m);
  }

  m = PyImport_ImportModule("smds.Cores");
  if (m)
  {
    PyObject *o = PyObject_GetAttrString(m, "CoreTypes");
    if (o)
    {
      PyDict_SetItemString(o, str(CORE_NAME), (PyObject*)&smds_coreType);
      Py_DECREF(o);
    }
    o = PyObject_GetAttrString(m, "paramsToCore");
    if (o)
    {
      if (ParamsType)
        PyDict_SetItem(o, ParamsType, (PyObject*)&smds_coreType);
      Py_DECREF(o);
    }
    Py_DECREF(m);
  }
}

static extent * smds_AHL_extent_new(int num)
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

static void smds_AHL_extent_free(extent *p)
{
  if (!p) return;
  smds_AHL_extent_free(p->next);
  free(p->data);
  free(p);
}

// takes one arg: Params
int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds)
{
  static  char *kwlist[] = { "p", NULL };
  smds_core *c = (smds_core*)self;
  PyObject *p = NULL, *o = NULL;
  double dT, dX, channelR, capH, capW;
  int i;
  
  if (!ParamsType)
  {
    PyErr_SetString(PyExc_SystemError, "Unable to find Params type.");
    return -1;
  }

  seedRnd(Random(NULL), &c->rs);

  // p is borrowed
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, 
  					&PyInstance_Type, &p))
    return -1;
  {
    PyObject *cls = PyObject_GetAttrString(p, "__class__");
    if (!cls) return -1;
    i = PyObject_Compare(cls, ParamsType);
    Py_DECREF(cls);
    if (PyErr_Occurred()) return -1;
    if (i != 0)
    {
      PyErr_SetString(PyExc_TypeError,
           "Argument must be an instance of the correct smds.Params type.");
      return -1;
    }
  }

  c->ready = 0;
  c->num = c->hits = c->misses = 0;
  c->timeToHitDist = smds_AHL_extent_new(500);
  if (c->timeToHitDist == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    return -1;
  }

  o = PyObject_GetAttrString(p, "dT");
  if (!o) return -1;
  dT = PyFloat_AsDouble(o);			// ns
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "binWidth");
  if (!o) return -1;
  c->bw = PyFloat_AsDouble(o);			// ns
  c->stepsPerBin = (int)(c->bw/dT+0.5);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "D");		// cm^2/s
  if (!o) return -1;
  dX = sqrt(6*PyFloat_AsDouble(o)*dT*1e-9)*1e7;	// nm
#ifdef CORE_FP
  c->stepsize = dX;
#endif
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "channelRadius");
  if (!o) return -1;
  channelR = PyFloat_AsDouble(o);
  c->channelRadius = (PosType)(channelR*channelR/dX/dX + ROUND);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "capWidth");
  if (!o) return -1;
  capW = PyFloat_AsDouble(o);
  c->capRadius = (PosType)((channelR + capW)*(channelR + capW)/dX/dX + ROUND);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "capHeight");
  if (!o) return -1;
  capH = PyFloat_AsDouble(o);
  c->capHeight = (PosType)(capH/dX+ROUND);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "initHeight");
  if (!o) return -1;
  c->initHeight = (PosType)((capH+PyFloat_AsDouble(o))/dX+ROUND);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "hitPoint");
  if (!o) return -1;
  c->hitPoint = (PosType)((capH-PyFloat_AsDouble(o))/dX+ROUND);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "timeToMiss");		// ms
  if (!o) return -1;
  c->stepsToMiss = (Int64)(PyFloat_AsDouble(o)*1e6/dT+0.5);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  c->ready = 1;
  return 0;
}

// takes three arguments: numbins, array[short], pos
//    PERFORMS NO CHECKING OF PARAMETERS!
PyObject * smds_core_run(PyObject *self, PyObject *args)
{
  smds_core *s = (smds_core*)self;
  PyArrayObject *o;
  extent *e;
  PosType x, y, z, r, R;
  int dur, pos, i, hit = 0;
  unsigned long tt;
  int dir, rndCounter = 0;
  unsigned int rnd = 0;
  Int64 t, stepsToMiss;
  PosType capHeight, capRadius, channelRadius, hitPoint;
#ifdef CORE_FP
  double stepsize, step;
#define STEP step
#else
#define STEP 1
#endif

  if (!PyArg_ParseTuple(args, "iO!i", &dur, &PyArray_Type, &o, &pos))
    return NULL;

  stepsToMiss = s->stepsToMiss;
  capHeight = s->capHeight;
  capRadius = s->capRadius;
  channelRadius = s->channelRadius;
  hitPoint = s->hitPoint;
#ifdef CORE_FP
  stepsize = s->stepsize;
#endif

  s->num += dur;

  Py_BEGIN_ALLOW_THREADS
  for (i=0; i < dur; i++)
  {
    x = y = r = 0;
    z = s->initHeight;
    for (t=0; t < stepsToMiss; t++)
    {
      /* Reposition */
#ifdef CORE_FP
    step = fabs(Normal(0, 1, &s->rs));
#endif
      do {
        if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
        dir = rnd & 0x7;  rndCounter--;
        rnd >>= 3;
      } while (dir == 6 || dir == 7);
      switch (dir)
      {
        case 0: // x++;
          R = (x+STEP)*(x+STEP)+y*y;
          if (!((z <= capHeight && r <= channelRadius && R > channelRadius) ||
                (z <= capHeight && r >= capRadius && R < capRadius)
             )) { x += STEP; r = R; }
          break;
        case 1: // x--;
          R = (x-STEP)*(x-STEP)+y*y;
          if (!((z <= capHeight && r <= channelRadius && R > channelRadius) ||
                (z <= capHeight && r >= capRadius && R < capRadius)
             )) { x -= STEP; r = R; }
          break;
        case 2: // y++;
          R = (y+STEP)*(y+STEP)+x*x;
          if (!((z <= capHeight && r <= channelRadius && R > channelRadius) ||
                (z <= capHeight && r >= capRadius && R < capRadius)
             )) { y += STEP; r = R; }
          break;
        case 3: // y--;
          R = (y-STEP)*(y-STEP)+x*x;
          if (!((z <= capHeight && r <= channelRadius && R > channelRadius) ||
                (z <= capHeight && r >= capRadius && R < capRadius)
             )) { y -= STEP; r = R; }
          break;
        case 4: // z++;
          z += STEP;
          break;
        case 5: // z--;
          if (!((r >= capRadius && z-STEP < 0) ||
                (r < capRadius && r > channelRadius && z-STEP < capHeight)
             )) z -= STEP;
          break;
      }
      
      /* Check for hit */
      if (r <= channelRadius && z <= hitPoint)
      {
        hit = 1;
        break;
      }
    } // each step

    if (hit)
    {
      hit = 0;
      s->hits++;
      tt = t/s->stepsPerBin;
      e = s->timeToHitDist;
      while (tt > e->num)
      {
        if (e->next == NULL)
        {
          e->next = smds_AHL_extent_new(500);
          if (e->next == NULL)
          {
            PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
            return NULL;
          } // error check
        } // create new extent
        tt -= e->num;
        e = e->next;
      } // find correct extent
      e->data[tt]++;
    } else  s->misses++;
  } // each rep

  Py_END_ALLOW_THREADS
  Py_RETURN_NONE;
}

PyObject * smds_core_getResults(PyObject *self)
{
  smds_core *c = (smds_core*)self;
  PyObject *r, *o, *ts = NULL, *counts = NULL;
  extent *e;
  int i, j, base, num;
  
  if (!ResultsType)
  {
    PyErr_SetString(PyExc_SystemError, "Couldn't find ResultsType.");
    return NULL;
  }
  r = PyInstance_New(ResultsType, NULL, NULL);
  if (!r) return NULL;

  o = PyInt_FromLong(c->num);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "num", o) == -1) goto bail;
  Py_DECREF(o);

  o = PyInt_FromLong(c->hits);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "hits", o) == -1) goto bail;
  Py_DECREF(o);

  o = PyInt_FromLong(c->misses);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "misses", o) == -1) goto bail;
  Py_DECREF(o);
  o = NULL;

  j = base = num = 0;
  e = c->timeToHitDist;
  while (e != NULL)
  {
    for (i=0; i < e->num; i++)
      if (e->data[i] > 0) num++;
    e = e->next;
  }

  ts = PyList_New(num);
  if (!ts) goto bail;
  counts = PyList_New(num);
  if (!counts) goto bail;
  e = c->timeToHitDist;
  while (e != NULL)
  {
    for (i=0; i < e->num; i++)
      if (e->data[i] > 0)
      {
        o = PyFloat_FromDouble((base+i)*c->bw*1e-9);
        if (!o) goto bail;
        if (PyList_SetItem(ts, j, o)) goto bail;	// steals ref
        
        o = PyInt_FromLong(e->data[i]);
        if (!o) goto bail;
        if (PyList_SetItem(counts, j, o)) goto bail;	// steals ref
        
        j++;
      }
    base += e->num;
    e = e->next;
  }
  o = NULL;
  if (PyObject_SetAttrString(r, "hitTimeDist_t", ts) == -1) goto bail;
  Py_DECREF(ts);
  ts = NULL;
  if (PyObject_SetAttrString(r, "hitTimeDist_count", counts) == -1) goto bail;
  Py_DECREF(counts);
  counts = NULL;

  return r;
bail:
  Py_DECREF(r);
  if (o) { Py_DECREF(o); }
  if (ts) { Py_DECREF(ts); }
  if (counts) { Py_DECREF(counts); }
  return NULL;    
}

void smds_core_free(smds_core *p)
{
  if (!p) return;
  smds_AHL_extent_free(p->timeToHitDist);
  p->ob_type->tp_free((PyObject*)p);
}

PyObject * smds_core_getParamsType(PyObject *cls)
{
  if (ParamsType) { Py_INCREF(ParamsType); }
  return ParamsType;
}

PyObject * smds_core_getResultsType(PyObject *cls)
{
  if (ResultsType) { Py_INCREF(ResultsType); }
  return ResultsType;
}

PyObject * smds_core_getName(PyObject *cls)
{
  return PyString_FromString(str(CORE_NAME));
}

PyObject * smds_core_sRnd(PyObject *self, PyObject *args)
{
  if (PyTuple_Size(args) != 1)
  {
    initRnd();
  }
  else
  {
    int seed;
    if (!PyArg_ParseTuple(args, "k", &seed))
      return NULL;
    seedRnd(seed, NULL);
  }

  Py_RETURN_NONE;
}
