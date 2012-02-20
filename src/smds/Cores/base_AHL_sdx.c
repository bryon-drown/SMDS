/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Alpha-Hemolysin Channel Simulation Cores
 *
 *  You must define CORE_NAME as an identifier.
 *  Optionally define CORE_DOC to a double-quoted doc-string.
 *  For numerical potential profiles, define CORE_POT.
 *  For electroosmotic flow, define CORE_OSM
 *  For spherical particles, define CORE_SPH
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
#define PI (3.141592653589793)

#include "random.h"
#include <math.h>

typedef unsigned long long Int64;
typedef int Int32;
typedef short Int16;
typedef double PosType;

typedef struct _smds_AHL_extent extent;
struct _smds_AHL_extent {
  int num;
  Int32 *data;
  extent *next;
};

#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
typedef struct _smds_AHL_coordinate coordinate;
struct _smds_AHL_coordinate {    // Linked list of positions of particles
  PosType t, x, y, z;
  coordinate *next;
};
#endif

typedef struct {
  PyObject_HEAD
  int ready;			// Initialization succeeded.
  double bw;
  Int64 timeToMiss;
  Int32 stepsPerBin;
  PosType capHeight, capRadius, channelRadius, initHeight, hitPoint, offset;
  PosType vestibuleRadius, vestibuleTop, vestibuleDepth;
#ifdef CORE_SPH
  PosType rP;
#endif
#ifdef CORE_OSM
  double eoDrift;
#endif
#ifdef CORE_POT
  PyArrayObject *potA, *potB;		// Potential profile
  int potNAz, potNAr, potNBz, potNBr;	// Number of bins in PotProfile
  double potSC_Az, potSC_Ar, potSC_Bz, potSC_Br;  // Dim conv for PotProfile
  PosType vA, vB;
#else
  double potExtent;
  PosType vA, vB, vC;
#endif
#ifdef CORE_1TAPE
  int stepsPerTape;
#endif
#ifdef CORE_2TAPE
  int stepsPerTape, tapeSelect;
#endif
  double dT, dX;
  /*** WORKING SPACE ***/
  RandomStream rs;
#ifdef CORE_1TAPE
  coordinate *tapeTail;
#endif
  /*** RESULTS ***/
  Int32 num, hits, misses;
  extent *timeToHitDist;
#ifdef CORE_1TAPE
  coordinate *tape;
#endif
#ifdef CORE_2TAPE
  coordinate *hitTape;
  coordinate *missTape;
#endif
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
#if (defined(CORE_1TAPE) || defined(CORE_2TAPE))
static coordinate * smds_AHL_coordinate_new(PosType t, PosType x, PosType y, PosType z);
static void smds_AHL_coordinate_free(coordinate *);
#endif

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
#ifdef CORE_POT
static PyObject *PotProfileType = NULL;
#endif

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

#ifdef CORE_POT
  m = PyImport_ImportModule("smds.Task");
  if (m)
  {
    PotProfileType = PyObject_GetAttrString(m, "PotProfile");
    Py_DECREF(m);
  }
#endif

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
  r->data = (Int32*)calloc(sizeof(Int32),num);
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

#if(defined(CORE_1TAPE) || defined(CORE_2TAPE))
static coordinate * smds_AHL_coordinate_new(PosType t, PosType x, PosType y, PosType z)
{
    coordinate *new;
    new = (coordinate*)malloc(sizeof(coordinate));
    new->t = t;
    new->x = x;
    new->y = y;
    new->z = z;
    new->next = NULL;
    return new;
}

static void smds_AHL_coordinate_free(coordinate *toFree)
{
    if (!toFree) return;
    smds_AHL_coordinate_free(toFree->next);
    free(toFree);
}
#endif

// takes one arg: Params
int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds)
{
  static  char *kwlist[] = { "p", NULL };
  smds_core *c = (smds_core*)self;
  PyObject *p = NULL, *o = NULL;
  double dT, dX, channelR, capH, capW, potA, potB, mobility;
  double vestibuleHeight, vestibuleRadius;
#ifdef CORE_OSM
  double g, sel, Nw, totPot, J, V;
#endif
#ifndef CORE_POT
  double potC;
#endif
  int i;

  if (!ParamsType)
  {
    PyErr_SetString(PyExc_SystemError, "Unable to find Params type.");
    return -1;
  }

#ifdef CORE_POT
  if (!PotProfileType)
  {
    PyErr_SetString(PyExc_SystemError, "Unable to find PotProfile type.");
    return -1;
  }
#endif

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
  c->dT = dT;
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
  c->dX = dX;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "channelRadius"); // nm
  if (!o) return -1;
  channelR = PyFloat_AsDouble(o);
  c->channelRadius = (PosType)(channelR/dX); // nm
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "capWidth");
  if (!o) return -1;
  capW = PyFloat_AsDouble(o);
  c->capRadius = (PosType)((channelR + capW)/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "capHeight");
  if (!o) return -1;
  capH = PyFloat_AsDouble(o);
  c->capHeight = (PosType)(capH/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "initHeight");
  if (!o) return -1;
  c->initHeight = (PosType)((capH+PyFloat_AsDouble(o))/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "hitPoint");
  if (!o) return -1;
  c->hitPoint = (PosType)((capH-PyFloat_AsDouble(o))/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "offset");
  if (!o) return -1;
  c->offset = (PosType)((PyFloat_AsDouble(o))/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "timeToMiss");		// ms
  if (!o) return -1;
  c->timeToMiss = (Int64)(PyFloat_AsDouble(o)*1e6/dT+0.5);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "vestibuleRadius");
  if (!o) return -1;
  vestibuleRadius = (PosType)((PyFloat_AsDouble(o)- channelR)/dX);
  c->vestibuleRadius = vestibuleRadius;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "vestibulePos");
  if (!o) return -1;
  c->vestibuleTop = c->capHeight - (PosType)(PyFloat_AsDouble(o)/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "vestibuleHeight");
  if (!o) return -1;
  vestibuleHeight = (PosType)(PyFloat_AsDouble(o)/dX);
  c->vestibuleDepth = c->vestibuleTop - vestibuleHeight;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

#ifdef CORE_SPH
  o = PyObject_GetAttrString(p, "particleRadius");
  if (!o) return -1;
  c->rP = (PosType)(PyFloat_AsDouble(o)/dX);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  // Correction for spherical particle
  c->vestibuleDepth += c->rP;
  c->vestibuleRadius = vestibuleRadius + c->rP - c->rP*sqrt(vestibuleRadius*vestibuleRadius + vestibuleHeight*vestibuleHeight)/vestibuleHeight;
  c->vestibuleTop -= (c->rP/vestibuleRadius)*(sqrt(vestibuleRadius*vestibuleRadius + vestibuleHeight*vestibuleHeight) - vestibuleHeight);
  c->channelRadius -= c->rP;
#endif
#ifdef CORE_OSM
  o = PyObject_GetAttrString(p, "conductance");		// pS
  if (!o) return -1;
  g = PyFloat_AsDouble(o) * 1e-12;                      // S
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "ionicSelectivity");    // K/Cl
  if (!o) return -1;
  sel = PyFloat_AsDouble(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "Nw");                  // particles/ion
  if (!o) return -1;
  Nw = PyFloat_AsDouble(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
#endif

  o = PyObject_GetAttrString(p, "mobility");		// nm^2 / Vs
  if (!o) return -1;
  mobility = PyFloat_AsDouble(o) / (dX*dX);
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

#ifdef CORE_POT
{
  PyObject *s = PyObject_GetAttrString(p, "pot");
  i = 1;
  if (!s || !PyInstance_Check(s) ||
      !(o = PyObject_GetAttrString(s, "__class__")) ||
      PyObject_Compare(PotProfileType, o) != 0 || (i = (int)PyErr_Occurred()))
  {
    if (!i)  PyErr_SetString(PyExc_TypeError,
             "Expected an instance of type smds.Task.PotProfile for pot.");
    if (o) Py_DECREF(o);
    if (s) Py_DECREF(s);
    return -1;
  }
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "dataA");
  if (!o)
  {
    PyErr_SetString(PyExc_ValueError, "Potential Profile lacks data for Zone A!");
    Py_DECREF(s);
    return -1;
  }
  c->potA = (PyArrayObject*)PyArray_ContiguousFromObject(o,
  				PyArray_DOUBLE, 1, 1);
			  	// new ref
  Py_DECREF(o);
  if (PyErr_Occurred())
  {
    PyErr_SetString(PyExc_ValueError, "Unintelligible Potential profile.");
    Py_DECREF(s);
    return -1;
  }

  o = PyObject_GetAttrString(s, "dataB");
  if (!o)
  {
    PyErr_SetString(PyExc_ValueError, "Potential Profile lacks data for Zone B!");
    Py_DECREF(s);
    return -1;
  }
  c->potB = (PyArrayObject*)PyArray_ContiguousFromObject(o,
  				PyArray_DOUBLE, 1, 1);
			  	// new ref
  Py_DECREF(o);
  if (PyErr_Occurred())
  {
    PyErr_SetString(PyExc_ValueError, "Unintelligible Potential profile.");
    Py_DECREF(s);
    return -1;
  }

  o = PyObject_GetAttrString(s, "bwAz");	// nm
  if (o) {
    c->potSC_Az = dX / PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "bwAr");
  if (o) {
    c->potSC_Ar = dX / PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "bwBz");
  if (o) {
    c->potSC_Bz = dX / PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "bwBr");
  if (o) {
    c->potSC_Br = dX / PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "NAz");
  if (o) {
    c->potNAz = PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "NAr");
  if (o) {
    c->potNAr = PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "NBz");
  if (o) {
    c->potNBz = PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  o = PyObject_GetAttrString(s, "NBr");
  if (o) {
    c->potNBr = PyFloat_AsDouble(o);
    Py_DECREF(o); }
  if (!o || PyErr_Occurred())  {  Py_DECREF(s);  return -1;  }

  Py_DECREF(s);
}

  c->vA = mobility * dT * 1e-9 * potA * dX;
  c->vB = mobility * dT * 1e-9 * potB * dX;
#ifdef CORE_OSM
  totPot = potA + potB;
#endif
#else
  o = PyObject_GetAttrString(p, "potC");
  if (!o) return -1;
  potC = PyFloat_AsDouble(o) * 1e-3;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "potExtent");
  if (!o) return -1;
  c->potExtent = -PyFloat_AsDouble(o)/dX;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  /* Calculation of stepSize due to electrophoresis.
   * drift speed = mobility * electric field
   */
  c->vA = mobility * dT * 1e-9 * potA / c->potExtent;	// < 0
  c->vB = mobility * dT * 1e-9 * (potB) /
                                     (c->vestibuleDepth - c->capHeight);
  c->vC = mobility * dT * 1e-9 * (potC) /
                                     (c->vestibuleDepth - 10/dX);   // pore is 10 nm tall
#ifdef CORE_OSM
  totPot = potA + potB + potC;
#endif
#endif
#ifdef CORE_OSM
  /* Calculation of electroosmotic drift as determined by Equation 7 in
   * Gu et al., PNAS, vol. 100, no. 26, pp 15498-15503
   */
  J = fabs(Nw*g*totPot/1.60217646e-19*((sel-1)/(sel+1))); // particles/s
  V = PI*pow(channelR,2)*10*1e-24*55*6.022e23; // particles/channel pi*r^2*h*conversion*55 mol/L*Avogadro's number
  c->eoDrift = J/V*10*1e-7;    // cm/s
  c->eoDrift = c->eoDrift * 1e-2 * dT /dX;  // step size due to EOF
#endif

#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
  o = PyObject_GetAttrString(p, "stepsPerTape");
  if(!o) return -1;
  c->stepsPerTape = (int)(PyFloat_AsDouble(o)/dT + 0.5);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

#ifdef CORE_2TAPE
  o = PyObject_GetAttrString(p, "tapeSelect");
  if(!o) return -1;
  c->tapeSelect = (int)(PyInt_AsLong(o));
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  c->missTape = c->hitTape = NULL;
#else
  c->tapeTail = c->tape = smds_AHL_coordinate_new(0,0,0,0);
#endif
#endif

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
    PosType x, y, z, r;
    double step;
    int dur, pos, i, hit = 0;
    unsigned long tt;
    int dir, rndCounter = 0;
    unsigned int rnd = 0;
    Int64 t, timeToMiss;
    PosType capHeight, capRadius, channelRadius, hitPoint;
    PosType vestibuleRadius, vestibuleTop, vestibuleDepth;
#ifdef CORE_SPH
    PosType rP = s->rP; // radius of particle
#endif
#ifdef CORE_POT
    PosType vA = s->vA, vB = s->vB;
    double *potA = (double*)(s->potA->data);
    double *potB = (double*)(s->potB->data);
    double scAz = s->potSC_Az, scAr = s->potSC_Ar,
         scBz = s->potSC_Bz, scBr = s->potSC_Br;
    double dZ, dR, dX, dY;
    int NAz = s->potNAz, NAr = s->potNAr, NBz = s->potNBz, NBr = s->potNBr;
#else
    double potExtent = s->potExtent;
    PosType Z, vA = s->vA, vB = s->vB, vC = s->vC;
#endif
#ifdef CORE_OSM
    double eoDrift = s->eoDrift;
#endif
#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
    int tapeCheck;
    int stepsPerTape = s->stepsPerTape;
#ifdef CORE_1TAPE
    coordinate *tail = s->tapeTail;
#else
    coordinate *tail, *tempTape;
    tail = tempTape = smds_AHL_coordinate_new(0,0,0,0); // create empty head for temp list
#endif
#endif

    inline int legalPos(PosType x, PosType y, PosType z)
    {
#ifdef CORE_SPH
    if (z > capHeight + rP) return 1;    // above cap
        r = sqrt(x*x + y*y);
    return (z >= rP && r >= capRadius + rP) ||    // sides of cap
            (r <= channelRadius) || // inside channel
            (z <= vestibuleTop && z >= vestibuleDepth && r <= (vestibuleTop-z)*
            vestibuleRadius / (vestibuleTop-vestibuleDepth) + channelRadius);   // in triangular vestibule
#else
    if (z > capHeight) return 1;    // above cap
    r = sqrt(x*x + y*y);    
    return (z >= 0 && r >= capRadius) ||    // sides of cap
            (r <= channelRadius) || // inside channel
            (z <= vestibuleTop && z >= vestibuleDepth && r <= (vestibuleTop-z)*
            vestibuleRadius / (vestibuleTop-vestibuleDepth) + channelRadius);   // in triangular vestibule
#endif
    }

    if (!PyArg_ParseTuple(args, "iO!i", &dur, &PyArray_Type, &o, &pos))
        return NULL;

    timeToMiss = s->timeToMiss;
    capHeight = s->capHeight;
    capRadius = s->capRadius;
    channelRadius = s->channelRadius;
    hitPoint = s->hitPoint;
    vestibuleRadius = s->vestibuleRadius;
    vestibuleTop = s->vestibuleTop;
    vestibuleDepth = s->vestibuleDepth;

    s->num += dur;

    Py_BEGIN_ALLOW_THREADS
    for (i=0; i < dur; i++)   // iterations
    {
        x = 0;
        y = s->offset;
        z = s->initHeight;
#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
        tapeCheck = stepsPerTape;
#endif
        for (t=0; t < timeToMiss; t++)  // time since launch
        {
            /* Diffusion */
            step = fabs(Normal(0, 1, &s->rs));    //random step size
            do {
                if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
                dir = rnd & 0x7;  rndCounter--;
                rnd >>= 3;
            } while (dir == 6 || dir == 7);
            switch (dir)
            {
                case 0: // x++;
                    if (legalPos(x+step, y, z))  x += step;
                    break;
                case 1: // x--;
                    if (legalPos(x-step, y, z)) x -= step;
                    break;
                case 2: // y++;
                    if (legalPos(x, y+step, z))  y += step;
                    break;
                case 3: // y--;
                    if (legalPos(x, y-step, z)) y -= step;
                    break;
                case 4: // z++;
                    if (legalPos(x, y, z+step)) z += step;
                    break;
                case 5: // z--;
                    if (legalPos(x, y, z-step)) z -= step;
                    break;
            }

            /* Electrophoretic Flow */
            r = sqrt(x*x + y*y);
#ifdef CORE_POT
            dZ = dR = 0;
            if (z > capHeight)			// Zone A (out of channel)
            {
                int iZ = (int)(z * scAz);
                int iR = (int)(r * scAr);
                if (iZ < NAz && iR < NAr)
                {
                    dZ = vA * potA[iZ * NAr + iR];
                    dR = vA * potA[iZ * NAr + iR + NAz*NAr];
                }
            } else {				// Zone B (in channel)
                int iZ = (int)((capHeight-z) * scBz);  // iZ axis points down!
                int iR = (int)(r * scBr);
                if (iZ < NBz && iR < NBr)
                {
                    dZ = vB * potB[iZ * NBr + iR];
                    dR = vB * potB[iZ * NBr + iR + NBz*NBr];
                }
            }

            if (r > 0)
            {
                dX = dR * x / r;
                dY = dR * y / r;
            } else
                dX = dY = 0;

            // euf... this really isn't the best way to handle bumping into a
            // wall, but for well-condiitoned potential profiles, provided the
            // step size is small enough and the potential isn't that large,
            // we probably shouldn't need to invoke the safety mechanism.
            if (legalPos(x+dX, y+dY, z+dZ))
            {
                x += dX;
                y += dY;
                z += dZ;
            }
#else
            if (z > capHeight)			// Zone A
            {
                double norm, u, v;
                Z = z-capHeight;  // Z > 0

                if (r < channelRadius)			// Down to mouth
                {
                    z += vA * exp(Z/potExtent);  // Note: vA, potExtent < 0
                } else {				// Toward mouth edge
                    u = (r-channelRadius);
                    norm = sqrt(u*u + Z*Z);
                    v = vA * exp(norm/potExtent);	// Note: vA, potExtent < 0

                    x += x * u / (r * norm) * v;
                    y += y * u / (r * norm) * v;
                    z += Z / norm * v;
                }
            } else if (r <= vestibuleRadius + channelRadius)  // if not above cap but inside pore
            {
                if (z >= vestibuleDepth)	// Zone B
                {
                    z += vB;
                    if (z < vestibuleDepth && r > channelRadius) // if particle goes through vestibule bottom
                        {  z = vestibuleDepth; }    // set particle to right on vestibule bottom
                } else {				// Zone C
                    z += vC;
                }
            } // End Electrophoretic Flow
#endif

#ifdef CORE_OSM
            /* Electroosmotic Flow */
            r = sqrt(x*x + y*y);
            if(z > capHeight)
            {
                double norm, u, v;
                Z = z-capHeight;  // Z > 0
                if(r < channelRadius)
                {
                    z -= eoDrift * channelRadius*channelRadius/(2*pow(channelRadius+Z, 2));
                } else {
                    u = (r-channelRadius);
                    norm = sqrt(u*u + Z*Z); // distance from edge of channel

                    v = eoDrift * channelRadius*channelRadius/(2*pow(channelRadius+norm, 2));

                    x -= x * u / (r * norm) * v;
                    y -= y * u / (r * norm) * v;
                    z -= Z / norm * v;
                }
            } else if(r <= vestibuleRadius + channelRadius) // if not above cap but inside pore
            {
                z -= eoDrift;
                if (z < vestibuleDepth && r> channelRadius)
                {
                    z = vestibuleDepth;
                }
            }
#endif
#if defined(CORE_2TAPE) || defined(CORE_1TAPE)
            /* Record to tape */
            tapeCheck--;
            if (tapeCheck == 0)
            {
#ifdef CORE_2TAPE
                tapeCheck = stepsPerTape;     // uncomment to allow for more than one snapshot
#endif
#ifdef CORE_1TAPE
                if (x*s->dX > -0.25 && x*s->dX < 0.25)  // only record if particle is within 0.5 nm slice
                {
#endif
                    tail->next = smds_AHL_coordinate_new(t,x,y,z);
                    tail = tail->next;
                    if(tail == NULL)
                    {
                        PyErr_SetString(PyExc_MemoryError, "Memory allocation error: coordinates.");
                        return NULL;
                    }
#ifdef CORE_1TAPE
                }
#endif
            }
#endif

            /* Check for hit */
            if (z <= hitPoint && sqrt(x*x+y*y) <= channelRadius)
            {
                hit = 1;
                break;
            }
        } // each timeStep

        /* Record Hit */
        if (hit)
        {
            hit = 0;

            /* Hit Count */
            s->hits++;

#ifdef CORE_2TAPE
            /* Hit Tape */
            if(s->tapeSelect <= 1)
            {
                tail->next = s->hitTape;
                s->hitTape = tempTape->next;
                tempTape->next = NULL;
                tail = tempTape;
            } else {smds_AHL_coordinate_free(tempTape->next);}
            tempTape->next = NULL;
            tail = tempTape;
#endif

            /* Time to Hit Distribution */
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
        } else {
            s->misses++;
#ifdef CORE_2TAPE
            /* Miss Tape */
            if(s->tapeSelect >= 1)
            {
                tail->next = s->missTape;
                s->missTape = tempTape->next;
            } else {smds_AHL_coordinate_free(tempTape->next);}
            tempTape->next = NULL;
            tail = tempTape;
#endif
        }
    } // each iteration

#ifdef CORE_1TAPE
    s->tapeTail = tail;
#endif

  Py_END_ALLOW_THREADS
  Py_RETURN_NONE;
}

PyObject * smds_core_getResults(PyObject *self)
{
  smds_core *c = (smds_core*)self;
  PyObject *r, *o, *ts = NULL, *counts = NULL;
  extent *e;
  int i, j, base, num;
#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
  PyObject *x = NULL, *y = NULL, *z = NULL, *t = NULL;
  coordinate *tape;
#endif

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

#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
  /* Hit Tape */
  num = 0;
  i = 0;
  t = PyList_New(0);
  if(!t) goto bail;
  x = PyList_New(0);
  if(!x) goto bail;
  y = PyList_New(0);
  if(!y) goto bail;
  z = PyList_New(0);
  if(!z) goto bail;
#ifdef CORE_2TAPE
  tape = c->hitTape;
#else
  tape = c->tape;
  tape = tape->next;    // skip empty node
#endif
  while (tape != NULL){
    o = PyFloat_FromDouble(tape->t*c->dT);
    if (!o) goto bail;
    if (PyList_Append(t, o)) goto bail;
    Py_DECREF(o);
    o = NULL;

    o = PyFloat_FromDouble(tape->x*c->dX);
    if (!o) goto bail;
    if (PyList_Append(x, o)) goto bail;
    Py_DECREF(o);
    o = NULL;

    o = PyFloat_FromDouble(tape->y*c->dX);
    if (!o) goto bail;
    if (PyList_Append(y, o)) goto bail;
    Py_DECREF(o);
    o = NULL;

    o = PyFloat_FromDouble(tape->z*c->dX);
    if (!o) goto bail;
    if (PyList_Append(z, o)) goto bail;
    Py_DECREF(o);
    o = NULL;
    tape = tape->next;
  }
#ifdef CORE_1TAPE
  if (PyObject_SetAttrString(r, "tape_t", t) == -1) goto bail;
  Py_DECREF(t);
  t = NULL;
  if (PyObject_SetAttrString(r, "tape_x", x) == -1) goto bail;
  Py_DECREF(x);
  x = NULL;
  if (PyObject_SetAttrString(r, "tape_y", y) == -1) goto bail;
  Py_DECREF(y);
  y = NULL;
  if (PyObject_SetAttrString(r, "tape_z", z) == -1) goto bail;
  Py_DECREF(z);
  z = NULL;
#else
  if (PyObject_SetAttrString(r, "hitTape_t", t) == -1) goto bail;
  Py_DECREF(t);
  t = NULL;
  if (PyObject_SetAttrString(r, "hitTape_x", x) == -1) goto bail;
  Py_DECREF(x);
  x = NULL;
  if (PyObject_SetAttrString(r, "hitTape_y", y) == -1) goto bail;
  Py_DECREF(y);
  y = NULL;
  if (PyObject_SetAttrString(r, "hitTape_z", z) == -1) goto bail;
  Py_DECREF(z);
  z = NULL;

  /* Miss Tape */
  num = 0;
  i = 0;
  t = PyList_New(0);
  if(!t) goto bail;
  x = PyList_New(0);
  if(!x) goto bail;
  y = PyList_New(0);
  if(!y) goto bail;
  z = PyList_New(0);
  if(!z) goto bail;
  tape = c->missTape;
  while (tape != NULL){
    o = PyFloat_FromDouble(tape->t*c->dT);
    if (!o) goto bail;
    if (PyList_Append(t, o)) goto bail;
    Py_DECREF(o);
    o = NULL;

    o = PyFloat_FromDouble(tape->x*c->dX);
    if (!o) goto bail;
    if (PyList_Append(x, o)) goto bail;
    Py_DECREF(o);
    o = NULL;

    o = PyFloat_FromDouble(tape->y*c->dX);
    if (!o) goto bail;
    if (PyList_Append(y, o)) goto bail;
    Py_DECREF(o);
    o = NULL;

    o = PyFloat_FromDouble(tape->z*c->dX);
    if (!o) goto bail;
    if (PyList_Append(z, o)) goto bail;
    Py_DECREF(o);
    o = NULL;
    tape = tape->next;
  }
  if (PyObject_SetAttrString(r, "missTape_t", t) == -1) goto bail;
  Py_DECREF(t);
  t = NULL;
  if (PyObject_SetAttrString(r, "missTape_x", x) == -1) goto bail;
  Py_DECREF(x);
  x = NULL;
  if (PyObject_SetAttrString(r, "missTape_y", y) == -1) goto bail;
  Py_DECREF(y);
  y = NULL;
  if (PyObject_SetAttrString(r, "missTape_z", z) == -1) goto bail;
  Py_DECREF(z);
  z = NULL;
#endif
#endif

  return r;
bail:
  Py_DECREF(r);
  if (o) { Py_DECREF(o); }
  if (ts) { Py_DECREF(ts); }
  if (counts) { Py_DECREF(counts); }
#if defined(CORE_1TAPE) || defined(CORE_2TAPE)
  if (t) { Py_DECREF(t); }
  if (x) { Py_DECREF(x); }
  if (y) { Py_DECREF(y); }
  if (z) { Py_DECREF(z); }
#endif
  return NULL;
}

void smds_core_free(smds_core *p)
{
  if (!p) return;
  smds_AHL_extent_free(p->timeToHitDist);
#ifdef CORE_POT
  if (p->potA) Py_DECREF((PyObject*)p->potB);
  if (p->potB) Py_DECREF((PyObject*)p->potA);
#endif
#ifdef CORE_2TAPE
  smds_AHL_coordinate_free(p->missTape);
  smds_AHL_coordinate_free(p->hitTape);
#endif
#ifdef CORE_1TAPE
  smds_AHL_coordinate_free(p->tape);
#endif
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
