/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Simulation Cores with Flow
 *    $Id: base_flow.c,v 1.4 2008/09/26 20:58:04 mculbert Exp $
 *
 *  You must define CORE_NAME as an identifier.
 *  Optionally define CORE_DOC to a double-quoted doc-string.
 *  Defaults to 2D core.  For 3D define CORE_3D.
 *  For intens profile, define CORE_INTENS.
 *  For photobleaching, define CORE_BLEACH.
 *  Note that these cores require that their corresponding Params
 *    and Results have the same name (CORE_NAME), but this isn't
 *    strictly necessary.  I chose this constraint for these cores
 *    to facilitate the generation of multiple cores from a single
 *    source file based on different defined macros.
 *  Flow is in the X direction, either pos or neg.
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

#include <math.h>
#include "random.h"

typedef int Int32;
typedef short Int16;

typedef struct _smds_molec smds_molec;
struct _smds_molec {
  double x, y;
#ifdef CORE_3D
  double z;
#endif
#ifdef CORE_BLEACH
  double photonTolerance;
#endif
  int species;
  smds_molec *next, *prev;
};

typedef struct {
  PyObject_HEAD
  int ready;			// Initialization succeeded.
  int numSpecies, numMolecs;
  double *X;			// fraction of each species
  double *a;			// a = Imax * dT
  double bkg, bw;		// bkg in phot/bin; bw in ms/bin
#ifdef CORE_BLEACH
  int *tol;
#endif
  int steps;			// Steps per bin = binWidth / dT
  Int32 *size;			// Space dimension in units of dX
  				// dX = sqrt(4*D*dT)
  				// size = sqrt(numMolecs/density) / dX
  Int32 *minX, *maxX;		// Locate the det vol in the flow axis (X)
  				// particle is on [-size, size), [minX, maxX)
#ifdef CORE_INTENS
  PyArrayObject *intens;
  int numBins_x;
  double *sc_x;          /* Size conversion factors for the intensity
                                   profile:
                                     sc_i[species] = dX / intens->binWidth_i
                                     I = (abs(x*sc_x) > intens->numBins_x ||
                                          abs(y*sc_x) > intens->numBins_x ||
                                          abs(z*sc_z) > intens->numBins_z ?
                                          0 : intens->data[
        (2*numBins_x+1) * (2*numBins_x+1) * (numBins_z + z * sc_z[spec]) +
        (2*numBins_x+1) * (numBins_x + y * sc_x[spec]) +
        (numBins_x + x * sc_x[spec]) ]
                                */
#ifdef CORE_3D
  double *sc_z;
  int numBins_z;
#endif
#else
  Int32 *threshold;		// Detection threshold in units of dX
  double *b;			// b = -2 / radius^2
#ifdef CORE_3D
  Int32 *threshold_z;
  double *c;			// c = -2 / radius_z^2
#endif
#endif
  double *flow_x;
  int dir_flow_x;
  double t_insert_avg;
  /*** WORKING SPACE ***/
  RandomStream rs;
  smds_molec *molecs, *avail;
  unsigned long int t_insert;
  /*** RESULTS ***/
  int dur;
  double I;
  int *counts;
#ifdef CORE_BLEACH
  int *numBleached;
#endif
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
#ifdef CORE_INTENS
static PyObject *IntensType = NULL;
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

#ifdef CORE_INTENS
  m = PyImport_ImportModule("smds.Task");
  if (m)
  {
    IntensType = PyObject_GetAttrString(m, "Intens");
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

// takes one arg: Params
int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds)
{
  static  char *kwlist[] = { "p", NULL };
  smds_core *c = (smds_core*)self;
  PyObject *p = NULL, *o = NULL, *oo = NULL, *s = NULL, *ss = NULL;
  double dX, dT, size, concentration;
  double preSize, postSize;
  int i;
#ifdef CORE_INTENS
  int binWidth_x;
#ifdef CORE_3D
  int binWidth_z;
#endif
#else
  int radius;
  int threshold;
#ifdef CORE_3D
  int Z;
  int threshold_Z;
#endif
#endif
  double flow_x;		// m/s

  if (!ParamsType)
  {
    PyErr_SetString(PyExc_SystemError, "Unable to find Params type.");
    return -1;
  }

#ifdef CORE_INTENS
  if (!IntensType)
  {
    PyErr_SetString(PyExc_SystemError, "Unable to find Intens type.");
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
  c->dur = 0;
  c->I = 0.0;

  o = PyObject_GetAttrString(p, "dT");
  if (!o) return -1;
  dT = PyFloat_AsDouble(o)*1e-9;		// s
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "binWidth");
  if (!o) return -1;
  c->bw = PyFloat_AsDouble(o);			// ms
  c->steps = (int)(c->bw * 1.0e-3 / dT + 0.5);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "bkg");
  if (!o) return -1;
  c->bkg = c->bw * PyFloat_AsDouble(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "numMolecs");
  if (!o) return -1;
  c->numMolecs = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "flowLength");
  if (!o) return -1;
  preSize = PyFloat_AsDouble(o);		// um
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "drainLength");
  if (!o) return -1;
  postSize = PyFloat_AsDouble(o);		// um
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "concentration");
  if (!o) return -1;
  concentration = PyFloat_AsDouble(o);
#ifdef CORE_3D
  size = sqrt((double)(c->numMolecs)/concentration/(preSize+postSize));
#else
  size = (double)(c->numMolecs)/concentration/(preSize+postSize);
#endif
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
#ifdef CORE_INTENS
  s = PyObject_GetAttrString(p, "intens");
  if (!s || !PyInstance_Check(s) ||
      !(o = PyObject_GetAttrString(s, "__class__")) ||
      PyObject_Compare(IntensType, o) != 0 || (i = (int)PyErr_Occurred()))
  {
    if (!i)  PyErr_SetString(PyExc_TypeError,
                      "Expected an instance of type smds.Task.Intens.");
    goto bail;
  }
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "data");
  if (!o)
  {
    PyErr_SetString(PyExc_ValueError, "Intensity Profile lacks data!");
    goto bail;
  }
  c->intens = (PyArrayObject*)PyArray_ContiguousFromObject(o, 
  				PyArray_DOUBLE, 1, 1);
			  	// c->intens is a new ref
  if (PyErr_Occurred())
  {
    PyErr_SetString(PyExc_ValueError, "Unintelligible Intensity profile.");
    goto bail;
  }
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "binWidth_x");
  if (!o) goto bail;
  binWidth_x = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_x");
  if (!o) goto bail;
  c->numBins_x = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);
#ifdef CORE_3D
  o = PyObject_GetAttrString(s, "binWidth_z");
  if (!o) goto bail;
  binWidth_z = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_z");
  if (!o) goto bail;
  c->numBins_z = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);
#endif
  Py_DECREF(s);  s = NULL;
#else
  o = PyObject_GetAttrString(p, "radius");
  if (!o) return -1;
  radius = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "threshold");
  if (!o) return -1;
  threshold = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
#ifdef CORE_3D
  o = PyObject_GetAttrString(p, "Z");
  if (!o) return -1;
  Z = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "threshold_Z");
  if (!o) return -1;
  threshold_Z = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
#endif
#endif

  o = PyObject_GetAttrString(p, "flow");		// m/s
  if (!o) return -1;
  flow_x = PyFloat_AsDouble(o) * 1e6;			// um/s
  if (flow_x == 0.0) c->dir_flow_x = 0;
  else c->dir_flow_x = (flow_x > 0.0 ? 1 : -1);
  flow_x = fabs(flow_x);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

#ifdef CORE_3D
  c->t_insert_avg = -1./(concentration*size*size*flow_x*dT);
#else
  c->t_insert_avg = -1./(concentration*size*flow_x*dT);
#endif
  do
  {
    c->t_insert = (int)(log(DRandom(&c->rs))*c->t_insert_avg);
  }
  while (!c->t_insert);

  ss = PyObject_GetAttrString(p, "SpeciesType");
  if (!ss) return -1;

  o = PyObject_GetAttrString(p, "sample");
  if (!o) return -1;
  if (PyList_Check(o))
  {
    PyObject *oo = PyList_AsTuple(o);
    Py_DECREF(o);
    o = oo;
  }
  if (!PyTuple_Check(o)) goto bail;
  c->numSpecies = PyTuple_Size(o);
  if (c->numSpecies <= 0)
  {
    PyErr_SetString(PyExc_ValueError, "You must have at least one species.");
    goto bail;
  }
  free(c->size); c->size = (Int32*)malloc(sizeof(Int32)*c->numSpecies);
  free(c->minX); c->minX = (Int32*)malloc(sizeof(Int32)*c->numSpecies);
  free(c->maxX); c->maxX = (Int32*)malloc(sizeof(Int32)*c->numSpecies);
  free(c->X);  c->X = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->a);  c->a = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->counts);  c->counts = (int*)malloc(sizeof(int)*c->numSpecies);
#ifdef CORE_INTENS
  free(c->sc_x);  c->sc_x = (double*)malloc(sizeof(double)*c->numSpecies);
#ifdef CORE_3D
  free(c->sc_z);  c->sc_z = (double*)malloc(sizeof(double)*c->numSpecies);
#endif
#else
  free(c->threshold);  
    c->threshold = (Int32*)malloc(sizeof(Int32)*c->numSpecies);
  free(c->b);  c->b = (double*)malloc(sizeof(double)*c->numSpecies);
#ifdef CORE_3D
  free(c->threshold_z);  
    c->threshold_z = (Int32*)malloc(sizeof(Int32)*c->numSpecies);
  free(c->c);  c->c = (double*)malloc(sizeof(double)*c->numSpecies);
#endif
#endif
#ifdef CORE_BLEACH
  free(c->tol);  c->tol = (int*)malloc(sizeof(int)*c->numSpecies);
  free(c->numBleached);
    c->numBleached = (int*)malloc(sizeof(int)*c->numSpecies);
#endif
  free(c->flow_x);
  c->flow_x = (double*)malloc(sizeof(double)*c->numSpecies);
  if (!(c->size && c->minX && c->maxX && c->X && c->a && c->counts
#ifdef CORE_INTENS
        && c->sc_x
#ifdef CORE_3D
        && c->sc_z
#endif
#else
        && c->threshold && c->b
#ifdef CORE_3D
        && c->threshold_z && c->c
#endif
#endif
#ifdef CORE_BLEACH
        && c->tol && c->numBleached
#endif
        && c->flow_x
     ))
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    goto bail;
  }
  for(i=0; i < c->numSpecies; i++)
  {
    int I = 0;
    s = PyTuple_GetItem(o, i);		// s is borrowed
    if (!s || !PyInstance_Check(s) || 
        !(oo = PyObject_GetAttrString(s, "__class__")) ||
        PyObject_Compare(ss, oo) != 0 || (I = (int)PyErr_Occurred()) )
    {
      if (!I)  PyErr_SetString(PyExc_TypeError,
      			"Wrong species type for this core.");
      goto bail;
    }
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "Fraction");
    if (!oo) goto bail;
    c->X[i] = PyFloat_AsDouble(oo);
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
#ifdef CORE_BLEACH
    c->numBleached[i] = 0;
    oo = PyObject_GetAttrString(s, "tolerance");
    if (!oo) goto bail;
    c->tol[i] = -PyFloat_AsDouble(oo);
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
#endif

                /* diffCoeff is in cm^2/s, 1.0e-4 converts to um^-1
                        multiply dX by a um distance to get 'dX'es  */
    oo = PyObject_GetAttrString(s, "D");
    if (!oo) goto bail;
#ifdef CORE_3D
    dX = 1.0 / sqrt(6.0 * PyFloat_AsDouble(oo) * dT) * 1e-4;
#else
    dX = 1.0 / sqrt(4.0 * PyFloat_AsDouble(oo) * dT) * 1e-4;
#endif
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "Imax");
    if (!oo) goto bail;
    c->a[i] = PyFloat_AsDouble(oo) * dT * 1.0e3;
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);  oo = NULL;

    c->size[i] = (int)((double)(size) * dX / 2.0);
    if (c->dir_flow_x < 0)
    {
      c->minX[i] = -(int)(postSize * dX);
      c->maxX[i] = (int)(preSize * dX);
    } else {
      c->minX[i] = -(int)(preSize * dX);
      c->maxX[i] = (int)(postSize * dX);
    }
    if (!( c->size[i] < sqrt((unsigned)(1 << 31)) &&
           c->maxX[i]-c->minX[i] < sqrt((unsigned)(1 << 31)) ))
    {
      PyErr_SetString(PyExc_ValueError, "Simulation space is too large.");
      goto bail;
    }

#ifdef CORE_INTENS
    c->sc_x[i] = 1.0 / binWidth_x / dX * 1.0e3;
#ifdef CORE_3D
    c->sc_z[i] = 1.0 / binWidth_z / dX * 1.0e3;
#endif
#else
    c->threshold[i] = (int)((double)(threshold * threshold) *
                    (dX * dX) * 1.0e-6);
    c->b[i] = -2.0 / (double)(radius * radius) / (dX * dX) * 1.0e6;
#ifdef CORE_3D
    c->threshold_z[i] = (int)( (double)(threshold_Z) * dX * 1.0e-3);
    c->c[i] = -2.0 / (double)(Z * Z) / (dX * dX) * 1.0e6;
#endif
#endif

    c->flow_x[i] = flow_x*dT*dX;

    c->counts[i] = 0;
    Py_DECREF(s);  s = NULL;
  }
  Py_DECREF(ss);  ss = NULL;

  if (c->molecs)
  {
    smds_molec *m = c->avail;
    if (!c->avail) c->avail = c->molecs;
    else 
    {
      while (m->next) m = m->next;
      m->next = c->molecs;
      c->molecs = NULL;
    }
  }
  for (i=0; i < c->numMolecs; i++)
  {
    smds_molec *m;
    if (c->avail)
    {
      m = c->avail;
      c->avail = c->avail->next;
    } else {
      m = (smds_molec*)malloc(sizeof(smds_molec));
    }
    m->prev = NULL;
    smds_core_create_molec(m, c, 0);
    if (c->molecs) c->molecs->prev = m;
    m->next = c->molecs;
    c->molecs = m;
  }

  c->ready = 1;
  return 0;
bail:
  if (o)  { Py_DECREF(o);  }
  if (oo) { Py_DECREF(oo); }
  if (s)  { Py_DECREF(s);  }
  if (ss) { Py_DECREF(ss); }
  return -1;
}

void smds_core_create_molec(smds_molec *m, smds_core *s, int edge)
{
  int r;
  double S;
   
  if (s->numSpecies > 1)
  {
    S = DRandom(&s->rs);
    for (r=0; r < s->numSpecies; r++)
      if (S < s->X[r])
        break;
      else
        S -= s->X[r];
    m->species = (r >= s->numSpecies ? s->numSpecies-1 : r);
  } else m->species = 0;
  s->counts[m->species]++;
                                                                                
#ifdef CORE_BLEACH
  /* Note that tolerance should be NEGATIVE */
  if (s->tol[m->species])
    m->photonTolerance = log(DRandom(&s->rs)) * s->tol[m->species];
#endif
                                                                                
  /* Place the new molecule at the edge */
  if (edge && (s->dir_flow_x || DRandom(&s->rs) < 0.5))
  {						// place on X extreme
    m->y = (int)((DRandom(&s->rs) - 0.5) * s->size[m->species] * 2.0);
    m->x = ( s->dir_flow_x < 0 ? s->maxX[m->species]-1
    			       : s->minX[m->species] );
#ifdef CORE_3D
    m->z = (int)((DRandom(&s->rs) - 0.5) * s->size[m->species] * 2.0);
#endif
  }
  else
  {						// place on Y extreme
    m->y = edge ? -s->size[m->species] :
           (int)((DRandom(&s->rs) - 0.5) * s->size[m->species] * 2.0);
    m->x = (int)(DRandom(&s->rs) * (s->maxX[m->species]-s->minX[m->species]) )
		+ s->minX[m->species];
#ifdef CORE_3D
    m->z = (int)((DRandom(&s->rs) - 0.5) * s->size[m->species] * 2.0);
#endif
  }
}

// takes three arguments: numbins, array[short], pos
//    PERFORMS NO CHECKING OF PARAMETERS!
PyObject * smds_core_run(PyObject *self, PyObject *args)
{
  smds_core *s = (smds_core*)self;
  PyArrayObject *o;
  Int16 *iData;
  double *data;
  int dur, pos, i;
  smds_molec *m;
  
  int size, minX, maxX;
  double a;
#ifdef CORE_INTENS
  double *intens = (double*)(s->intens->data);
  double sc_x;
  int numBins_x = s->numBins_x;
#ifdef CORE_3D
  double sc_z;
  int numBins_z = s->numBins_z;
#endif
#else
  int r;
  int threshold;
  double b;
#ifdef CORE_3D
  int threshold_z;
  double c;
#endif
#endif
  double flow_x;
  int dir_flow_x = s->dir_flow_x;
  double t_insert_avg = s->t_insert_avg;
  unsigned long int t_insert = s->t_insert;

  int t, dir, rndCounter = 0;
  unsigned int rnd = 0;
  double *binCount;
  double I;

  if (!PyArg_ParseTuple(args, "iO!i", &dur, &PyArray_Type, &o, &pos))
    return NULL;

  s->dur += dur;
  data = (double*)malloc(sizeof(double)*dur);
  if (!data)
  {
    PyErr_SetString(PyExc_MemoryError, "Unable to allocate working space.");
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS
  if (s->bkg > 0.0)
    for (i=0; i < dur; i++) data[i] = s->bkg;
  else
    for (i=0; i < dur; i++) data[i] = 0.0;

  binCount = data;
  for (i=dur; i; i--)
  {
    for (t=s->steps; t; t--, t_insert--)
    {
      /* Create new molecules */
      while (!t_insert)
      {
        if (s->avail)
        {
          m = s->avail;
          s->avail = s->avail->next;
        } else {
          m = (smds_molec*)malloc(sizeof(smds_molec));
        }
        m->prev = NULL;
        smds_core_create_molec(m, s, 1);
        if (s->molecs) s->molecs->prev = m;
        m->next = s->molecs;
        s->molecs = m;
        t_insert = (int)(log(DRandom(&s->rs))*t_insert_avg);
      }

      m = s->molecs;
      while (m)
      {
        size = s->size[m->species];
        minX = s->minX[m->species];
        maxX = s->maxX[m->species];
        a = s->a[m->species];
#ifdef CORE_INTENS
        sc_x = s->sc_x[m->species];
#ifdef CORE_3D
        sc_z = s->sc_z[m->species];
#endif
#else
        threshold = s->threshold[m->species];
        b = s->b[m->species];
#ifdef CORE_3D
        threshold_z = s->threshold_z[m->species];
        c = s->c[m->species];
#endif
#endif
        flow_x = s->flow_x[m->species];

        /* Flow */
        if ( ((dir_flow_x > 0) && ((m->x += flow_x) >= maxX)) ||
             ((dir_flow_x < 0) && ((m->x -= flow_x) <  minX)) )
        {
          smds_molec *p = m;
          m = m->next;
          if (p->prev) p->prev->next = p->next;
          else         s->molecs = p->next;		// head of list
          if (p->next) p->next->prev = p->prev;
          // p->prev = NULL;	s->avail is singly-linked.
          p->next = s->avail;
          s->avail = p;
          continue;
        }

        /* Reposition */
#ifdef CORE_3D
        do {
          if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
          dir = rnd & 0x7;  rndCounter--;
          rnd >>= 3;
        } while (dir == 6 || dir == 7);
#else
        if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 16; }
        dir = rnd & 0x3;  rndCounter--;
        rnd >>= 2;
#endif
        switch (dir)
        {
          case 0: m->x++;  if (m->x >= maxX) { m->x = minX;   }  break;
          case 1: m->x--;  if (m->x <  minX) { m->x = maxX-1; }  break;
          case 2: m->y++;  if (m->y >= size) { m->y = -size;  }  break;
          case 3: m->y--;  if (m->y < -size) { m->y = size-1; }  break;
#ifdef CORE_3D
          case 4: m->z++;  if (m->z >= size) { m->z = -size;  }  break;
          case 5: m->z--;  if (m->z < -size) { m->z = size-1; }  break;
#endif
        }

        /* Calculate Intensity */
#ifdef CORE_INTENS
#ifdef CORE_3D
        if (abs((int)(m->x*sc_x)) <= numBins_x &&
            abs((int)(m->y*sc_x)) <= numBins_x &&
            abs((int)(m->z*sc_z)) <= numBins_z)
        {
          I = a * intens[
            (2*numBins_x+1)*(2*numBins_x+1)*(int)(numBins_z + m->z * sc_z) +
            (2*numBins_x+1) * (int)(numBins_x + m->y * sc_x) +
            (int)(numBins_x + m->x * sc_x) ];
#else
        if (abs((int)(x*sc_x)) <= numBins_x &&
            abs((int)(y*sc_x)) <= numBins_x)
        {
          I = a * intens[
            (2*numBins_x+1) * (int)(numBins_x + m->y * sc_x) +
            (int)(numBins_x + m->x * sc_x)  ];
#endif
#else  // not INTENS
#ifdef CORE_3D
        if ( (r = m->x*m->x + m->y*m->y) < threshold && 
             abs(m->z) < threshold_z)
        {
          I = a * exp(b*(double)r + c*(double)(m->z*m->z));
#else
        if ( (r = m->x*m->x + m->y*m->y) < threshold)
        {
          I = a * exp(b * (double)r);
#endif
#endif
#ifdef CORE_BLEACH
          if ((m->photonTolerance -= I) < 0)
          {
            smds_molec *p = m;
            I += m->photonTolerance;
            s->numBleached[m->species]++;
            *binCount += I;

            m = m->next;
            if (p->prev) p->prev->next = p->next;
            else         s->molecs = p->next;		// head of list
            if (p->next) p->next->prev = p->prev;
            p->prev = NULL;
            p->next = s->avail;
            s->avail = p;
            continue;
          }
#endif
          *binCount += I;
        } // in threshold

        m = m->next;
      } // each molecule

    } // each step

    binCount++;
  } // each bin

  s->t_insert = t_insert;

  iData = (Int16*)o->data + pos;
  for (i=0; i < dur; i++)
    s->I += (double)(iData[i] = Poisson(data[i], &s->rs));

  Py_END_ALLOW_THREADS
  free(data);

  Py_RETURN_NONE;
}

PyObject * smds_core_getResults(PyObject *self)
{
  smds_core *c = (smds_core*)self;
  PyObject *r, *o, *l = NULL;
  int i;
  
  if (!ResultsType)
  {
    PyErr_SetString(PyExc_SystemError, "Couldn't find ResultsType.");
    return NULL;
  }
  r = PyInstance_New(ResultsType, NULL, NULL);
  if (!r) return NULL;

  o = PyFloat_FromDouble(c->I/(double)c->dur/c->bw);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "avg_I", o) == -1) goto bail;
  Py_DECREF(o);
  
  o = PyFloat_FromDouble((double)c->dur*c->bw*1e-3);
  if (!o) goto bail;
  if (PyObject_SetAttrString(r, "length", o) == -1) goto bail;
  Py_DECREF(o);
  
  l = PyList_New(c->numSpecies);
  if (!l) goto bail;
  for(i=0; i < c->numSpecies; i++)
  {
    o = PyInt_FromLong(c->counts[i]);
    if (!o) goto bail;
    if (PyList_SetItem(l, i, o)) goto bail;	// steals reference
  }
  o = NULL;
  if (PyObject_SetAttrString(r, "counts", l) == -1) goto bail;
  Py_DECREF(l);
  
#ifdef CORE_BLEACH
  l = PyList_New(c->numSpecies);
  if (!l) goto bail;
  for(i=0; i < c->numSpecies; i++)
  {
    o = PyInt_FromLong(c->numBleached[i]);
    if (!o) goto bail;
    if (PyList_SetItem(l, i, o)) goto bail;	// steals reference
  }
  o = NULL;
  if (PyObject_SetAttrString(r, "bleached", l) == -1) goto bail;
  Py_DECREF(l);
#endif
  
  return r;
bail:
  Py_DECREF(r);
  if (o) { Py_DECREF(o); }
  if (l) { Py_DECREF(l); }
  return NULL;    
}

void smds_core_free(smds_core *p)
{
  smds_molec *m;
  if (!p) return;
  free(p->X);
  free(p->a);
#ifdef CORE_BLEACH
  free(p->tol);
#endif
  free(p->size);
  free(p->minX);
  free(p->maxX);
#ifdef CORE_INTENS
  if (p->intens) {  Py_DECREF((PyObject*)p->intens);  }
  free(p->sc_x);
#ifdef CORE_3D
  free(p->sc_z);
#endif
#else
  free(p->threshold);
  free(p->b);
#ifdef CORE_3D
  free(p->c);
  free(p->threshold_z);
#endif
#endif
  free(p->counts);
#ifdef CORE_BLEACH
  free(p->numBleached);
#endif
  free(p->flow_x);
  while (p->molecs)
  {
    m = p->molecs;
    p->molecs = p->molecs->next;
    free(m);
  }
  while (p->avail)
  {
    m = p->avail;
    p->avail = p->avail->next;
    free(m);
  }
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
