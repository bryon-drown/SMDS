/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Alpha-Hemolysin Channel Simulation Cores
 *    $Id$
 *
 *  You must define CORE_NAME as an identifier.
 *  Optionally define CORE_DOC to a double-quoted doc-string.
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
typedef double PosType;

typedef struct { PosType x, y, z; } smds_molec;

typedef struct {
  PyObject_HEAD
  int ready;			// Initialization succeeded.
  double bw;
  Int32 stepsPerBin, microSteps;
  int numMolecs;
  PosType size, microStepSize;
  PosType capHeight, capRadius, channelRadius, hitPoint;
  PosType vestibuleRadius, vestibuleTop, vestibuleDepth;
//  double potA1, potB0, potB1, potC0, potC1, potExtent;
  double potExtent;
  PosType vA, vB, vC;
  /*** WORKING SPACE ***/
  RandomStream rs;
  smds_molec *m;
  /*** RESULTS ***/
  Int32 dur, hits;
} smds_core;

// takes one arg: Params
static int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds);
// takes three arguments: numbins, array, pos
//    PERFORMS NO CHECKING OF PARAMETERS!
static void smds_core_create_molec(smds_molec *m, smds_core *s, int edge);
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

void smds_core_create_molec(smds_molec *m, smds_core *s, int edge)
{
  /* Place the new molecule at the edge */
  if (edge && DRandom(&s->rs) < 0.5)
  {
    m->x = (PosType)((DRandom(&s->rs) - 0.5) * s->size * 2.0);
    m->y = -s->size;
    m->z = (PosType)(DRandom(&s->rs) * s->size * 2.0);
  }
  else
  {
    m->x = edge ? -s->size :
           (PosType)((DRandom(&s->rs) - 0.5) * s->size * 2.0);
    m->y = (PosType)((DRandom(&s->rs) - 0.5) * s->size * 2.0);
    m->z = (PosType)(DRandom(&s->rs) * s->size * 2.0);
  }

  if (!edge && m->z < s->capHeight &&
      (m->x*m->x + m->y*m->y) < s->capRadius)
    m->z += s->capHeight;
}

// takes one arg: Params
int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds)
{
  static  char *kwlist[] = { "p", NULL };
  smds_core *c = (smds_core*)self;
  PyObject *p = NULL, *o = NULL;
  double dT, dX, channelR, capH, capW, potA, potB, potC, mobility;
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
  c->dur = 0;

  o = PyObject_GetAttrString(p, "macro_dT");
  if (!o) return -1;
  dT = PyFloat_AsDouble(o);			// ns
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "binWidth");
  if (!o) return -1;
  c->bw = PyFloat_AsDouble(o);			// us
  c->stepsPerBin = (int)(c->bw*1e3/dT+0.5);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "numMolecs");
  if (!o) return -1;
  c->numMolecs = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "D");		// cm^2/s
  if (!o) return -1;
  dX = sqrt(6*PyFloat_AsDouble(o)*dT*1e-9)*1e7;	// nm
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "concentration");
  if (!o) return -1;
  c->size = (PosType)(pow(c->numMolecs/PyFloat_AsDouble(o), 1./3.) * 1e3
                      / dX) / 2;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "micro_dT");		// ns
  if (!o) return -1;
  c->microSteps = (int)(dT / PyFloat_AsDouble(o) + 0.5);
  dT /= c->microSteps;
  c->microStepSize = 1/sqrt(c->microSteps);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  free(c->m);  c->m = (smds_molec*)malloc(sizeof(smds_molec)*c->numMolecs);
  if (!c->m)
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    return -1;
  }
  for (i=0; i < c->numMolecs; i++)
    smds_core_create_molec(c->m+i, c, 0);

  o = PyObject_GetAttrString(p, "channelRadius");
  if (!o) return -1;
  channelR = PyFloat_AsDouble(o);
  c->channelRadius = (PosType)(channelR*channelR/dX/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "capWidth");
  if (!o) return -1;
  capW = PyFloat_AsDouble(o);
  c->capRadius = (PosType)((channelR + capW)*(channelR + capW)/dX/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "capHeight");
  if (!o) return -1;
  capH = PyFloat_AsDouble(o);
  c->capHeight = (PosType)(capH/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "hitPoint");
  if (!o) return -1;
  c->hitPoint = (PosType)((capH-PyFloat_AsDouble(o))/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "vestibuleRadius");
  if (!o) return -1;
  c->vestibuleRadius = (PosType)(PyFloat_AsDouble(o)/dX);
  c->vestibuleRadius *= c->vestibuleRadius;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "vestibulePos");
  if (!o) return -1;
  c->vestibuleTop = c->capHeight - (PosType)(PyFloat_AsDouble(o)/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "vestibuleHeight");
  if (!o) return -1;
  c->vestibuleDepth = c->vestibuleTop - (PosType)(PyFloat_AsDouble(o)/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "potExtent");
  if (!o) return -1;
  c->potExtent = PyFloat_AsDouble(o)/dX;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "potA");
  if (!o) return -1;
  potA = PyFloat_AsDouble(o) * 1e-3;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "potB");
  if (!o) return -1;
  potB = PyFloat_AsDouble(o) * 1e-3;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "potC");
  if (!o) return -1;
  potC = PyFloat_AsDouble(o) * 1e-3;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "mobility");
  if (!o) return -1;
  mobility = PyFloat_AsDouble(o) / (dX*dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

/*
  c->potA1 = potA / c->potExtent;
  c->potB1 = (potB - potA) / (c->vestibuleDepth - c->capHeight);
  c->potC1 = (potC - potB) / (c->hitPoint - c->vestibuleDepth);
  c->potB0 = potA - c->potB1 * c->capHeight;
  c->potC0 = potB - c->potC1 * c->vestibuleDepth;
*/
  c->vA = mobility * dT * 1e-9 * (potA-0) / (0-c->potExtent);
  c->vB = mobility * dT * 1e-9 * (potB-potA) /
                                     (c->vestibuleDepth - c->capHeight);
  c->vC = mobility * dT * 1e-9 * (potC-potB) /
                                     (c->hitPoint - c->vestibuleDepth);
  
  c->ready = 1;
  return 0;
}

// takes three arguments: numbins, array[short], pos
//    PERFORMS NO CHECKING OF PARAMETERS!
PyObject * smds_core_run(PyObject *self, PyObject *args)
{
  smds_core *s = (smds_core*)self;
  PyArrayObject *o;
  Int16 *iData;
  PosType x, y, z, r, Z;
  double step;
  int dur, pos, i, ii, m, micro;
  Int64 t;
  int dir, rndCounter = 0;
  unsigned int rnd = 0;
  PosType size, microStepSize;
  PosType capHeight, capRadius, channelRadius, hitPoint;
  PosType vestibuleRadius, vestibuleTop, vestibuleDepth;
//  double potA1, potB0, potB1, potC0, potC1, potExtent;
  double potExtent;
  PosType vA, vB, vC;
  
  inline int legalPos(PosType x, PosType y, PosType z)
  {
    PosType r = x*x + y*y;
    if (z > capHeight) return 1;
    if (z >= 0 && r >= capRadius) return 1;
    if (z <= vestibuleTop && z >= vestibuleDepth && r <= vestibuleRadius)
      return 1;
    if (r <= channelRadius) return 1;
    return 0;
  }

  if (!PyArg_ParseTuple(args, "iO!i", &dur, &PyArray_Type, &o, &pos))
    return NULL;

  size = s->size;
  microStepSize = s->microStepSize;
  capHeight = s->capHeight;
  capRadius = s->capRadius;
  channelRadius = s->channelRadius;
  hitPoint = s->hitPoint;
  vestibuleRadius = s->vestibuleRadius;
  vestibuleTop = s->vestibuleTop;
  vestibuleDepth = s->vestibuleDepth;
/*
  potA1 = s->potA1;
  potB0 = s->potB0;
  potB1 = s->potB1;
  potC0 = s->potC0;
  potC1 = s->potC1;
*/  
  potExtent = s->potExtent;
  vA = s->vA;
  vB = s->vB;
  vC = s->vC;

  s->dur += dur;
  iData = (Int16*)o->data + pos;
  for (i=0; i < dur; i++)  iData[i] = 0;

  Py_BEGIN_ALLOW_THREADS
  for (m=0; m < s->numMolecs; m++)
  {
    x = s->m[m].x;
    y = s->m[m].y;
    z = s->m[m].z;
    
    for (t=0; t < dur; t++)
    {
      for (i=s->stepsPerBin; i; i--)
      {
        micro = (x*x+y*y+z*z <= 1);
        for (ii=(micro?s->microSteps:1); ii; ii--)
        {
          /* Reposition */
          step = fabs(Normal(0, 1, &s->rs));
          if (micro) step *= microStepSize;
          do {
            if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
            dir = rnd & 0x7;  rndCounter--;
            rnd >>= 3;
          } while (dir == 6 || dir == 7);
          switch (dir)
          {
            case 0: // x++;
              if (legalPos(x+step, y, z) && x+step < size)  x += step;
              break;
            case 1: // x--;
              if (legalPos(x-step, y, z) && x-step > -size) x -= step;
              break;
            case 2: // y++;
              if (legalPos(x, y+step, z) && y+step < size)  y += step;
              break;
            case 3: // y--;
              if (legalPos(x, y-step, z) && y-step > -size) y -= step;
              break;
            case 4: // z++;
              if (legalPos(x, y, z+step) && z+step < 2*size) z += step;
              break;
            case 5: // z--;
              if (legalPos(x, y, z-step)) z -= step;
              break;
          }
          
          /* Electrophoretic Flow */
  /*
          if (z > capHeight)
          {
            Z = z-capHeight;
            r = sqrt(x*x + y*y + Z*Z);
            if (r < potExtent && r > 0)		// Zone A
            {
              delta = (r - potExtent) * potA1 / r;
              x += x * delta;
              y += y * delta;
              z += Z * delta;
              if (z < capHeight)		// overshot
              { x = y = 0; z = capHeight; }
            }
          } else if (z >= vestibuleDepth) {	// Zone B
            z += potB0 + potB1 * z;
            if (z < vestibuleDepth && x*x + y*y > channelRadius)
            {  z = vestibuleDepth; }
          } else {				// Zone C
            z += potC0 + potC1 * z;
          }
  */
          if (z > capHeight)
          {
            Z = z-capHeight;
            r = sqrt(x*x + y*y + Z*Z);
            if (r < potExtent && r > 0)		// Zone A
            {
              x += x / r * vA;
              y += y / r * vA;
              z += Z / r * vA;
              if (z < capHeight)		// overshot
              { x = y = 0; z = capHeight; }
            }
          } else if (z >= vestibuleDepth) {	// Zone B
            z += vB;
            if (z < vestibuleDepth && x*x + y*y > channelRadius)
            {  z = vestibuleDepth; }
          } else {				// Zone C
            z += vC;
          }
          
          /* Check for hit */
          if (z <= hitPoint && x*x+y*y <= channelRadius)
          {
            iData[t]++;
            s->hits++;
            smds_core_create_molec(s->m+m, s, 1);
            x = s->m[m].x;
            y = s->m[m].y;
            z = s->m[m].z;
            break;		// terminate micro loop
          }
        } // each microstep
      } // each step
    } // each bin

    s->m[m].x = x;
    s->m[m].y = y;
    s->m[m].z = z;
    
  } // each molec     

  Py_END_ALLOW_THREADS
  Py_RETURN_NONE;
}

PyObject * smds_core_getResults(PyObject *self)
{
  smds_core *c = (smds_core*)self;
  PyObject *r, *o;
  
  if (!ResultsType)
  {
    PyErr_SetString(PyExc_SystemError, "Couldn't find ResultsType.");
    return NULL;
  }
  r = PyInstance_New(ResultsType, NULL, NULL);
  if (!r) return NULL;

  o = PyFloat_FromDouble(c->dur*c->bw*1e-6);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "length", o) == -1) goto bail;
  Py_DECREF(o);

  o = PyInt_FromLong(c->hits);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "hits", o) == -1) goto bail;
  Py_DECREF(o);

  return r;
bail:
  Py_DECREF(r);
  if (o) { Py_DECREF(o); }
  return NULL;    
}

void smds_core_free(smds_core *p)
{
  if (!p) return;
  free(p->m);
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
