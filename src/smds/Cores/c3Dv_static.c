/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Vesicle Simulation Cores
 *    $Id: base_vesicle.c,v 1.5 2008/09/26 20:58:04 mculbert Exp $
 *
 *  You must define CORE_NAME as an identifier.
 *  Optionally define CORE_DOC to a double-quoted doc-string.
 *  For intens profile, define CORE_INTENS.
 *  For photobleaching, define CORE_BLEACH.
 *  For blinking, define CORE_BLINK.
 *  For intersystem crossing, define CORE_TRIPLET.
 *  Note that these cores require that their corresponding Params
 *    and Results have the same name (CORE_NAME), but this isn't
 *    strictly necessary.  I chose this constraint for these cores
 *    to facilitate the generation of multiple cores from a single
 *    source file based on different defined macros.
 *
 *******************************************/

#if defined(CORE_TRIPLET) && !defined(CORE_INTENS)
#error "You must define CORE_INTENS to use CORE_TRIPLET."
#endif

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
#define SQR(x) ((x)*(x))

#include "random.h"
#include <math.h>

typedef int Int32;
typedef short Int16;

typedef struct {
  double x, y, z;
#ifdef CORE_BLEACH
  double photonTolerance;
#endif
#ifdef CORE_BLINK
  int state;
  int transitionTime;
#endif
#ifdef CORE_TRIPLET
  unsigned int tripletRelaxTime;
#endif
  int species;
} smds_molec;

typedef struct {
  double x, y, z;
  double r, r2;
  int inThresh;
  int numMolecs, maxMolecs;
  smds_molec *m;
  int type;
} smds_vesicle;

typedef struct {
  PyObject_HEAD
  int ready;			// Initialization succeeded.
  int numSpecies, numVesicleTypes;
  int numVesicles, numMolecs;	// free molecules
  double *vNumMolecs;		// avg num molecs/vesicle for each type
  double *vR, *vSigma;		// vesicle radius mean and S.D.
  double *X, *vX, *vmX;		// frac of free mol, ves types, mol in ves
  double *a;			// a = Imax * dT
  double bkg, bw;		// bkg in phot/bin; bw in ms/bin
#ifdef CORE_BLEACH		// tol is based on collected photons
  int *tol;			// unless using CORE_TRIPLET, in which
#endif				// case, tol is based on emitted photons
#ifdef CORE_BLINK		// For blinking:
  double *k_d, *k_b;		// dark/bright equilibrium constants
#endif
#ifdef CORE_TRIPLET
  double *k_T, *k_S0;		// Rate constants for S1->T and T->S0
  double phi;			// Collection efficiency function magnitude
  PyArrayObject *cef;		// Collection efficiency function profile
  int cef_numBins_x;		// when CORE_TRIPLET is used, intens is only
  double cef_sc_x;		// the laser intensity profile, not the
  int cef_numBins_z;		// entire molecule detection efficiency func
  double cef_sc_z;
#endif // CORE_TRIPLET
  int steps;			// Steps per bin = binWidth / dT
  double *dX;			// dX for each species
  double size;			// Space dimension for free mol
  				// size = sqrt(numMolecs/density)
  double vInsert, vFree;	// vesicle insertion and destruction points
  double *vdX;			// vesicle step size
#ifdef CORE_INTENS
  PyArrayObject *intens;
  int numBins_x;
  double sc_x;          /* Size conversion factors for the intensity
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
  double sc_z;
  int numBins_z;
#else
  double b, c;			// b = -2 / radius^2, c = -2 / radius_z^2
#endif
  double threshold, threshold_z;	// Detection threshold
  /*** WORKING SPACE ***/
  RandomStream rs;
  smds_molec *m;		// Free molecules
  smds_vesicle *v;		// Vesicles
  /*** RESULTS ***/
  int dur;
  double I;
  int *counts, *vCounts, *numFree;
#ifdef CORE_BLEACH
  int *numBleached;
#endif
} smds_core;

// takes one arg: Params
static int smds_core_init(PyObject *self, PyObject *args, PyObject *kwds);
// takes three arguments: numbins, array, pos
//    PERFORMS NO CHECKING OF PARAMETERS!
static void smds_core_create_molec(smds_molec *m, smds_core *s, 
					smds_vesicle *v, int edge);
static void smds_core_create_vesicle(smds_vesicle *v, smds_core *s, int edge);
static PyObject * smds_core_run(PyObject *self, PyObject *args);
static PyObject * smds_core_getResults(PyObject *self);
static void smds_core_free(smds_core *);
static PyObject * smds_core_getParamsType(PyObject *cls);
static PyObject * smds_core_getResultsType(PyObject *cls);
static PyObject * smds_core_getName(PyObject *cls);
static PyObject * smds_core_sRnd(PyObject *self, PyObject *args);
static double Gaussian(RandomStream *);

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
  smds_vesicle *v, *vv;
  PyObject *p = NULL, *o = NULL, *oo = NULL, *s = NULL, *ss = NULL;
  double dT, d, maxVdX, *x;
  int i, j;
  
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
  
  o = PyObject_GetAttrString(p, "concentration");
  if (!o) return -1;
  c->size = pow((double)(c->numMolecs)/PyFloat_AsDouble(o), 1.0/3.0) / 2.0;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "numVesicles");
  if (!o) return -1;
  c->numVesicles = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
  o = PyObject_GetAttrString(p, "vesicleConcentration");
  if (!o) return -1;
  c->vInsert = pow((double)(c->numVesicles)/PyFloat_AsDouble(o), 1.0/3) / 2.0;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
  
#ifdef CORE_TRIPLET
  o = PyObject_GetAttrString(p, "collectionEfficiency");
  if (!o) return -1;
  c->phi = PyFloat_AsDouble(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  s = PyObject_GetAttrString(p, "cef");
  if (!s || !PyInstance_Check(s) ||
      !(o = PyObject_GetAttrString(s, "__class__")) ||
      PyObject_Compare(IntensType, o) != 0 || (i = (int)PyErr_Occurred()))
  {
    if (!i)  PyErr_SetString(PyExc_TypeError,
             "Expected an instance of type smds.Task.Intens for CEF.");
    goto bail;
  }
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "data");
  if (!o)
  {
    PyErr_SetString(PyExc_ValueError, "CEF Profile lacks data!");
    goto bail;
  }
  c->cef = (PyArrayObject*)PyArray_ContiguousFromObject(o, 
  				PyArray_DOUBLE, 1, 1);
			  	// c->cef is a new ref
  if (PyErr_Occurred())
  {
    PyErr_SetString(PyExc_ValueError, "Unintelligible Intensity profile.");
    goto bail;
  }
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "binWidth_x");
  if (!o) goto bail;
  c->cef_sc_x = 1.0 / ( (double)PyInt_AsLong(o) * 1e-3 );
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_x");
  if (!o) goto bail;
  c->cef_numBins_x = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "binWidth_z");
  if (!o) goto bail;
  c->cef_sc_z = 1.0 / ( (double)PyInt_AsLong(o) * 1e-3 );
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_z");
  if (!o) goto bail;
  c->cef_numBins_z = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  Py_DECREF(s);  s = NULL;
#endif // CORE_TRIPLET

#ifdef CORE_INTENS
  s = PyObject_GetAttrString(p, "intens");
  if (!s || !PyInstance_Check(s) ||
      !(o = PyObject_GetAttrString(s, "__class__")) ||
      PyObject_Compare(IntensType, o) != 0 || (i = (int)PyErr_Occurred()))
  {
    if (!i)  PyErr_SetString(PyExc_TypeError,
             "Expected an instance of type smds.Task.Intens for intens.");
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

  o = PyObject_GetAttrString(s, "numBins_x");
  if (!o) goto bail;
  c->numBins_x = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "binWidth_x");
  if (!o) goto bail;
  d = (double)PyInt_AsLong(o) * 1e-3;			// um
  c->sc_x = 1.0 / d;
  d *= (double)c->numBins_x;
  c->threshold = d*d/2.0;
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_z");
  if (!o) goto bail;
  c->numBins_z = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "binWidth_z");
  if (!o) goto bail;
  d = (double)PyInt_AsLong(o) * 1e-3;
  c->sc_z = 1.0 / d;
  c->threshold_z = d * (double)c->numBins_z / 2.0;
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  Py_DECREF(s);  s = NULL;
#else
  o = PyObject_GetAttrString(p, "radius");
  if (!o) return -1;
  d = PyFloat_AsDouble(o) * 1e-3;
  c->b = -2.0 / (d*d);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "threshold");
  if (!o) return -1;
  d = PyFloat_AsDouble(o) * 1e-3;
  c->threshold = d*d;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "Z");
  if (!o) return -1;
  d = PyFloat_AsDouble(o) * 1e-3;
  c->c = -2.0 / (d*d);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "threshold_Z");
  if (!o) return -1;
  c->threshold_z = PyFloat_AsDouble(o) * 1e-3;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
#endif

  // Species
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
  free(c->X);  c->X = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->a);  c->a = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->dX); c->dX = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->counts);  c->counts = (int*)malloc(sizeof(int)*c->numSpecies);
  free(c->numFree); c->numFree = (int*)malloc(sizeof(int)*c->numSpecies);
#ifdef CORE_BLEACH
  free(c->tol);  c->tol = (int*)malloc(sizeof(int)*c->numSpecies);
  free(c->numBleached);
    c->numBleached = (int*)malloc(sizeof(int)*c->numSpecies);
#endif
#ifdef CORE_BLINK
  free(c->k_d); c->k_d = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->k_b); c->k_b = (double*)malloc(sizeof(double)*c->numSpecies);
#endif
#ifdef CORE_TRIPLET
  free(c->k_T);  c->k_T  = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->k_S0); c->k_S0 = (double*)malloc(sizeof(double)*c->numSpecies);
#endif
  free(c->m);  c->m = (smds_molec*)malloc(sizeof(smds_molec)*c->numMolecs);
  if (!(c->X && c->a && c->dX && c->counts && c->numFree && c->m
#ifdef CORE_BLEACH
        && c->tol && c->numBleached
#endif
#ifdef CORE_BLINK
        && c->k_d && c->k_b
#endif
#ifdef CORE_TRIPLET
        && c->k_T && c->k_S0
#endif
     ))
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    goto bail;
  }
  for(i=0; i < c->numSpecies; i++)
  {
    int I = 0;
    s = PyTuple_GetItem(o, i);		// s is borrowed
    Py_INCREF(s);
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

#ifdef CORE_BLINK
    oo = PyObject_GetAttrString(s, "k_dark");		// s^-1
    if (!oo) goto bail;
    c->k_d[i] = -1/(dT*PyFloat_AsDouble(oo));
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "k_bright");
    if (!oo) goto bail;
    c->k_b[i] = -1/(dT*PyFloat_AsDouble(oo));
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
#endif

#ifdef CORE_TRIPLET
    oo = PyObject_GetAttrString(s, "k_T");		// s^-1
    if (!oo) goto bail;
    c->k_T[i] = dT*PyFloat_AsDouble(oo);
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "k_S0");
    if (!oo) goto bail;
    c->k_S0[i] = dT*PyFloat_AsDouble(oo);
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
#endif

    oo = PyObject_GetAttrString(s, "D");
    if (!oo) goto bail;
    c->dX[i] = sqrt(6.0 * PyFloat_AsDouble(oo) * dT) * 1e4;	// um
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "Imax");
    if (!oo) goto bail;
    c->a[i] = PyFloat_AsDouble(oo) * dT * 1e3;
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);  oo = NULL;

    c->counts[i] = 0;
    c->numFree[i] = 0;
    Py_DECREF(s);  s = NULL;
  }
  Py_DECREF(o); o = NULL;
  Py_DECREF(ss);  ss = NULL;

  for (i=0; i < c->numMolecs; i++)
    smds_core_create_molec(&(c->m[i]), c, NULL, 0);

  // Vesicles
  ss = PyObject_GetAttrString(p, "VesicleType");
  if (!ss) return -1;

  o = PyObject_GetAttrString(p, "vesicles");
  if (!o) return -1;
  if (PyList_Check(o))
  {
    PyObject *oo = PyList_AsTuple(o);
    Py_DECREF(o);
    o = oo;
  }
  if (!PyTuple_Check(o)) goto bail;
  c->numVesicleTypes = PyTuple_Size(o);
  if (c->numVesicleTypes <= 0)
  {
    PyErr_SetString(PyExc_ValueError,
    			"You must have at least one vesicle type.");
    goto bail;
  }
  free(c->vNumMolecs);
    c->vNumMolecs = (double*)malloc(sizeof(double)*c->numVesicleTypes);
  free(c->vR);
    c->vR = (double*)malloc(sizeof(double)*c->numVesicleTypes);
  free(c->vSigma);
    c->vSigma = (double*)malloc(sizeof(double)*c->numVesicleTypes);
  free(c->vX);
    c->vX = (double*)malloc(sizeof(double)*c->numVesicleTypes);
  free(c->vmX);
    c->vmX = (double*)malloc(sizeof(double)*c->numVesicleTypes*c->numSpecies);
  free(c->vdX);
    c->vdX = (double*)malloc(sizeof(double)*c->numVesicleTypes);
  free(c->vCounts);
    c->vCounts = (int*)malloc(sizeof(int)*c->numVesicleTypes);
  free(c->v);
    c->v = (smds_vesicle*)malloc(sizeof(smds_vesicle)*c->numVesicles);
  if (!(c->vNumMolecs && c->vR && c->vSigma && c->vX && c->vmX &&
        c->vdX && c->vCounts && c->v
     ))
  {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
    goto bail;
  }
  for(i=0; i < c->numVesicleTypes; i++)
  {
    int I = 0;
    s = PyTuple_GetItem(o, i);		// s is borrowed
    Py_INCREF(s);
    if (!s || !PyInstance_Check(s) || 
        !(oo = PyObject_GetAttrString(s, "__class__")) ||
        PyObject_Compare(ss, oo) != 0 || (I = (int)PyErr_Occurred()) )
    {
      if (!I)  PyErr_SetString(PyExc_TypeError,
      			"Wrong vesicle type for this core.");
      goto bail;
    }
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "Fraction");
    if (!oo) goto bail;
    c->vX[i] = PyFloat_AsDouble(oo);
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "D");
    if (!oo) goto bail;
    c->vdX[i] = sqrt(6.0 * PyFloat_AsDouble(oo) * dT) * 1e4;	// um
    if (i == 0 || c->vdX[i] > maxVdX) maxVdX = c->vdX[i];
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);
    
    oo = PyObject_GetAttrString(s, "size");
    if (!oo) goto bail;
    c->vR[i] = PyFloat_AsDouble(oo);				// um
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);  oo = NULL;

    oo = PyObject_GetAttrString(s, "size_spread");
    if (!oo) goto bail;
    c->vSigma[i] = PyFloat_AsDouble(oo);			// um
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);  oo = NULL;

    oo = PyObject_GetAttrString(s, "numMolecs");
    if (!oo) goto bail;
    c->vNumMolecs[i] = PyFloat_AsDouble(oo);
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);  oo = NULL;

    oo = PyObject_GetAttrString(s, "composition");
    if (!oo) goto bail;
    if (PyList_Check(oo))
    {
      PyObject *tmp = PyList_AsTuple(oo);
      Py_DECREF(oo);
      oo = tmp;
    }
    if (!PyTuple_Check(oo)) goto bail;
    j = PyTuple_Size(oo);
    if (j != c->numVesicleTypes)
    {
      PyErr_SetString(PyExc_ValueError,
      			"Composition must match number of species.");
      goto bail;
    }
    x = c->vmX + i*c->numSpecies;
    for (j--; j >= 0; j--)
    {
      PyObject *n = PyTuple_GetItem(oo, j);		// n is borrowed
      x[j] = PyFloat_AsDouble(n);
      if (PyErr_Occurred()) goto bail;
    }
    Py_DECREF(oo); oo = NULL;

    c->vCounts[i] = 0;
    Py_DECREF(s);  s = NULL;
  }
  Py_DECREF(o); o = NULL;
  Py_DECREF(ss);  ss = NULL;

  c->vFree = c->vInsert + 5 * maxVdX;
  for (i=0, v=c->v; i < c->numVesicles; i++, v++)
  {
    int fail;
    v->maxMolecs = 0;
    v->m = NULL;
    do
    {
      fail = 0;
      smds_core_create_vesicle(v, c, 0);
      if (v->inThresh)
        for (j=0, vv=c->v; j < i && !fail; j++, vv++)
          if (vv->inThresh &&
              sqrt(SQR(v->x-vv->x) + SQR(v->y-vv->y) + SQR(v->z-vv->z))
                < v->r + vv->r)  fail = 1;
    } while (fail);
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

void smds_core_create_vesicle(smds_vesicle *v, smds_core *s, int edge)
{
  int r, num, type;
  double S;
  smds_molec *m;
  
  // Choose Type
  if (s->numVesicleTypes > 1)
  {
    S = DRandom(&s->rs);
    for (r=0; r < s->numVesicleTypes; r++)
      if (S < s->vX[r])
        break;
      else
        S -= s->vX[r];
    type = v->type = (r >= s->numVesicleTypes ? s->numVesicleTypes-1 : r);
  } else type = v->type = 0;
  s->vCounts[type]++;
  
  // Choose Size
  v->r = s->vR[type];
  v->r2 = v->r * v->r;

  // Position
  if (edge)
  {
    switch ((int)(DRandom(&s->rs)*6))
    {
      case 0:
        v->x = s->vInsert;
        v->y = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->z = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        break;
      case 1:
        v->x = -s->vInsert;
        v->y = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->z = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        break;
      case 2:
        v->x = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->y = s->vInsert;
        v->z = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        break;
      case 3:
        v->x = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->y = -s->vInsert;
        v->z = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        break;
      case 4:
        v->x = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->y = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->z = s->vInsert;
        break;
      case 5:
        v->x = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->y = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
        v->z = -s->vInsert;
        break;
    }
  }
  else
  {
    v->x = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
    v->y = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
    v->z = ((DRandom(&s->rs) - 0.5) * s->vInsert * 2.0);
  }
  
  v->inThresh = ((v->x*v->x + v->y*v->y - v->r2) < s->threshold &&
                 (fabs(v->z) - v->r) < s->threshold_z);
  
  // Create Molecules
  num = v->numMolecs = (int)(s->vNumMolecs[type]+0.5);
  if (num > v->maxMolecs)
  {
    free(v->m);
    v->maxMolecs = num;
    v->m = (smds_molec*)malloc(sizeof(smds_molec)*num);
  }
  for (r=0, m = v->m; r < num; r++, m++)
    smds_core_create_molec(m, s, v, 0);
}

void smds_core_create_molec(smds_molec *m, smds_core *s, 
				smds_vesicle *v, int edge)
{
  int r, sp;
  double S, *p;
   
  if (s->numSpecies > 1)
  {
    S = DRandom(&s->rs);
    if (v) p = s->vmX + v->type*s->numSpecies;
    else   p = s->X;
    for (r=0; r < s->numSpecies; r++, p++)
      if (S < *p)
        break;
      else
        S -= *p;
    sp = m->species = (r >= s->numSpecies ? s->numSpecies-1 : r);
  } else sp = m->species = 0;
  s->counts[sp]++;
  if (!v) s->numFree[sp]++;

#ifdef CORE_BLEACH
  /* Note that tolerance should be NEGATIVE */
  if (s->tol[sp])
    m->photonTolerance = log(DRandom(&s->rs)) * s->tol[sp];
#endif

#ifdef CORE_BLINK
  m->state = DRandom(&s->rs) > 
  	1/s->k_d[sp]/(1/s->k_d[sp]+1/s->k_b[sp]);
  m->transitionTime = (int)(log(DRandom(&s->rs))*
  		(m->state ? s->k_d[sp] : s->k_b[sp]) );
#endif

#ifdef CORE_TRIPLET
  m->tripletRelaxTime = 0;
#endif
                                                                                
  // Position
  if (v)
  {
    double theta, phi;

    theta = DRandom(&s->rs) * 2 * M_PI;
    phi = DRandom(&s->rs) * M_PI;
    S = DRandom(&s->rs) * v->r;
    
    m->x = S * cos(theta) * sin(phi);
    m->y = S * sin(theta) * sin(phi);
    m->z = S * cos(phi);
    
  } else {	// free molecule
    if (edge && DRandom(&s->rs) < 0.5)
    {
      m->x = ((DRandom(&s->rs) - 0.5) * s->size * 2.0);
      m->y = -s->size;
      m->z = ((DRandom(&s->rs) - 0.5) * s->size * 2.0);
    }
    else
    {
      m->x = edge ? -s->size :
             ((DRandom(&s->rs) - 0.5) * s->size * 2.0);
      m->y = ((DRandom(&s->rs) - 0.5) * s->size * 2.0);
      m->z = ((DRandom(&s->rs) - 0.5) * s->size * 2.0);
    }
  } // position molecule
  
} // create_molec

// takes three arguments: numbins, array[short], pos
//    PERFORMS NO CHECKING OF PARAMETERS!
PyObject * smds_core_run(PyObject *self, PyObject *args)
{
  smds_core *s = (smds_core*)self;
  PyArrayObject *o;
  Int16 *iData;
  double *data;
  int dur, pos, molec, vesicle, i, sp, tp, inThresh;
  int numMolecs = s->numMolecs, numVesicles = s->numVesicles;
  smds_molec *m;
  smds_vesicle *v;
  
#ifdef CORE_INTENS
  double *intens = (double*)(s->intens->data);
  double sc_x = s->sc_x, sc_z = s->sc_z;
  int numBins_x = s->numBins_x, numBins_z = s->numBins_z;
#else
  double r, b = s->b, c = s->c;
  double threshold = s->threshold, threshold_z = s->threshold_z;
#endif
#ifdef CORE_TRIPLET
  double cef_value, *cef = (double*)(s->cef->data);
  double cef_sc_x = s->cef_sc_x, cef_sc_z = s->cef_sc_z;
  int cef_numBins_x = s->cef_numBins_x, cef->numBins_z = s->cef_numBins_z;
#endif
  double size = s->size, vFree = s->vFree;

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
#if !defined(CORE_TRIPLET)
  if (s->bkg > 0.0)
    for (i=0; i < dur; i++) data[i] = s->bkg;
  else
#endif
    for (i=0; i < dur; i++) data[i] = 0.0;

  binCount = data;
  for (i=dur; i; i--)
  {
    for (t=s->steps; t; t--)
    {

      // Process free molecules
      for (molec=0, m=s->m; molec < numMolecs; molec++, m++)
      {
        sp = m->species;
      
        /* Reposition */
        do {
          if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
          dir = rnd & 0x7;  rndCounter--;
          rnd >>= 3;
        } while (dir == 6 || dir == 7);
        switch (dir)
        {
          case 0: m->x += s->dX[sp];
          		if (m->x >= size) { m->x = -size;  }  break;
          case 1: m->x -= s->dX[sp];
          		if (m->x < -size) { m->x = size-s->dX[sp]; }  break;
          case 2: m->y += s->dX[sp];
          		if (m->y >= size) { m->y = -size;  }  break;
          case 3: m->y -= s->dX[sp];
          		if (m->y < -size) { m->y = size-s->dX[sp]; }  break;
          case 4: m->z += s->dX[sp];
          		if (m->z >= size) { m->z = -size;  }  break;
          case 5: m->z -= s->dX[sp];
          		if (m->z < -size) { m->z = size-s->dX[sp]; }  break;
        }

        // Exclude motion that would bump into vesicles
        inThresh = 
#ifdef CORE_INTENS
            abs((int)(m->x*sc_x)) <= numBins_x &&
            abs((int)(m->y*sc_x)) <= numBins_x &&
            abs((int)(m->z*sc_z)) <= numBins_z);
#else  // not INTENS
            ( (r = m->x*m->x + m->y*m->y) < threshold && 
              fabs(m->z) < threshold_z );
#endif

        if (inThresh)
        {
          // Spatial exclusion for vesicles in threshold
          int i, fail = 0;
          smds_vesicle *vv = s->v;
          for (i=0; i < numVesicles; i++, vv++)
            if (vv->inThresh &&
                sqrt(SQR(m->x-vv->x) + SQR(m->y-vv->y) + SQR(m->z-vv->z))
                  < vv->r)
            {
              fail = 1;
              break;
            }
          if (fail)		// Refuse move
            switch (dir)	// Note, wraping shouldn't occur in the
            {			// threshold, so don't bother with it here.
              case 0:  m->x -= s->dX[sp];  break;
              case 1:  m->x += s->dX[sp];  break;
              case 2:  m->y -= s->dX[sp];  break;
              case 3:  m->y += s->dX[sp];  break;
              case 4:  m->z -= s->dX[sp];  break;
              case 5:  m->z += s->dX[sp];  break;
            }
        } // if inThresh

#ifdef CORE_TRIPLET
        if (m->tripletRelaxTime) m->tripletRelaxTime--;
#endif
#ifdef CORE_BLINK
        if (--(m->transitionTime) <= 0)
        {
          m->state = !m->state;
          m->transitionTime = (int)(log(DRandom(&s->rs)) *
  		(m->state ? s->k_d[sp] : s->k_b[sp]) );
        }
        if (!m->state) continue;		// next free molec
#endif
#ifdef CORE_TRIPLET
        if (m->tripletRelaxTime) continue;	// next free molec
#endif

        /* Calculate Intensity */
        if (inThresh)
        {
#ifdef CORE_INTENS
          I = s->a[sp] * intens[
            (2*numBins_x+1)*(2*numBins_x+1)*(int)(numBins_z + m->z * sc_z) +
            (2*numBins_x+1) * (int)(numBins_x + m->y * sc_x) +
            (int)(numBins_x + m->x * sc_x) ];
#else  // not INTENS
          I = s->a[sp] * exp(b*r + c*(m->z*m->z));
#endif

#ifdef CORE_TRIPLET
          if (abs((int)(m->x*cef_sc_x)) <= cef_numBins_x &&
              abs((int)(m->y*cef_sc_x)) <= cef_numBins_x &&
              abs((int)(m->z*cef_sc_z)) <= cef_numBins_z)
          {
            cef_value = cef[
              (2*cef_numBins_x+1) * (2*cef_numBins_x+1) * 
              	(int)(cef_numBins_z + m->z * cef_sc_z) +
              (2*cef_numBins_x+1) * (int)(cef_numBins_x + m->y * cef_sc_x) +
              (int)(cef_numBins_x + m->x * cef_sc_x) ];
          } else cef_value = 0.0;
#endif // CORE_TRIPLET

#ifdef CORE_BLEACH
          if ((m->photonTolerance -= I) < 0)
          {
            I += m->photonTolerance;
            s->numBleached[sp]++;
            smds_core_create_molec(m, s, NULL, 1);
#ifdef CORE_TRIPLET
          } else {
#endif
#endif // CORE_BLEACH
#ifdef CORE_TRIPLET
          {
            double p = I*s->k_T[sp]/(I*s->k_T[sp]+s->k_S0[sp]);
            if (p > 1e-10 && DRandom(&s->rs) < p)
              if (m->tripletRelaxTime =
                    (int)(-log(DRandom(&s->rs)) / s->k_S0[sp]))
                continue;
          }
#endif

#ifdef CORE_TRIPLET
          *binCount += I * cef_value;
#else
          *binCount += I;
#endif
        } // in threshold

      }  // each free molecule
  
      // Process vesicles
      for (vesicle=0, v=s->v; vesicle < numVesicles; vesicle++, v++)
      {
        tp = v->type;
      
        /* Reposition */
        do {
          if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
          dir = rnd & 0x7;  rndCounter--;
          rnd >>= 3;
        } while (dir == 6 || dir == 7);
        switch (dir)
        {
          case 0:
            v->x += s->vdX[tp];
            if (v->x >  vFree) { smds_core_create_vesicle(v, s, 1);  }
            break;
          case 1:
            v->x -= s->vdX[tp];
            if (v->x < -vFree) { smds_core_create_vesicle(v, s, 1);  }
            break;
          case 2:
            v->y += s->vdX[tp];
            if (v->y >  vFree) { smds_core_create_vesicle(v, s, 1);  }
            break;
          case 3:
            v->y -= s->vdX[tp];
            if (v->y < -vFree) { smds_core_create_vesicle(v, s, 1);  }
            break;
          case 4:
            v->z += s->vdX[tp];
            if (v->z >  vFree) { smds_core_create_vesicle(v, s, 1);  }
            break;
          case 5:
            v->z -= s->vdX[tp];
            if (v->z < -vFree) { smds_core_create_vesicle(v, s, 1);  }
            break;
        }
        
        v->inThresh = ((v->x*v->x + v->y*v->y - v->r2) < s->threshold &&
                       (fabs(v->z) - v->r) < s->threshold_z);

        if (v->inThresh)
        {
          // Spatial exclusion for vesicles in threshold
          int i, fail = 0;
          smds_vesicle *vv = s->v;
          for (i=0; i < numVesicles; i++, vv++)
            if (i != vesicle && vv->inThresh &&
                sqrt(SQR(v->x-vv->x) + SQR(v->y-vv->y) + SQR(v->z-vv->z))
                  < v->r + vv->r)
            {
              fail = 1;
              break;
            }
          if (fail)		// Refuse move
          {
            switch (dir)
            {
              case 0:  v->x -= s->vdX[tp];  break;
              case 1:  v->x += s->vdX[tp];  break;
              case 2:  v->y -= s->vdX[tp];  break;
              case 3:  v->y += s->vdX[tp];  break;
              case 4:  v->z -= s->vdX[tp];  break;
              case 5:  v->z += s->vdX[tp];  break;
            }
            v->inThresh = ((v->x*v->x + v->y*v->y - v->r2) < s->threshold &&
                       (fabs(v->z) - v->r) < s->threshold_z);
          }

          // Shove free molecules out of the way, if necessary
          for (i=0, m=s->m; i < numMolecs; i++, m++)
            if ((SQR(v->x-m->x) + SQR(v->y-m->y) + SQR(v->z-m->z)) < v->r2)
              switch (dir)
              {
                case 0:  m->x += s->vdX[tp];  break;
                case 1:  m->x -= s->vdX[tp];  break;
                case 2:  m->y += s->vdX[tp];  break;
                case 3:  m->y -= s->vdX[tp];  break;
                case 4:  m->z += s->vdX[tp];  break;
                case 5:  m->z -= s->vdX[tp];  break;
              }

        } // if v->inThresh

        // Process vesicle's molecules
        for (molec=0, m=v->m; molec < v->numMolecs; molec++, m++)
        {
          double x = v->x, y = v->y, z = v->z;
          sp = m->species;
        
          /* Reposition */
          do {
            if (!rndCounter) { rnd = Random(&s->rs);  rndCounter = 10; }
            dir = rnd & 0x7;  rndCounter--;
            rnd >>= 3;
          } while (dir == 6 || dir == 7);
          switch (dir)
          {
            case 0:
              m->x += s->dX[sp];
              if ((m->x*m->x + m->y*m->y + m->z*m->z) >= v->r2)
                m->x -= s->dX[sp];
              break;
            case 1:
              m->x -= s->dX[sp];
              if ((m->x*m->x + m->y*m->y + m->z*m->z) >= v->r2)
                m->x += s->dX[sp];
              break;
            case 2:
              m->y += s->dX[sp];
              if ((m->x*m->x + m->y*m->y + m->z*m->z) >= v->r2)
                m->y -= s->dX[sp];
              break;
            case 3:
              m->y -= s->dX[sp];
              if ((m->x*m->x + m->y*m->y + m->z*m->z) >= v->r2)
                m->y += s->dX[sp];
              break;
            case 4:
              m->z += s->dX[sp];
              if ((m->x*m->x + m->y*m->y + m->z*m->z) >= v->r2)
                m->z -= s->dX[sp];
              break;
            case 5:
              m->z -= s->dX[sp];
              if ((m->x*m->x + m->y*m->y + m->z*m->z) >= v->r2)
                m->z += s->dX[sp];
              break;
          }
        
  #ifdef CORE_TRIPLET
          if (m->tripletRelaxTime) m->tripletRelaxTime--;
  #endif
  #ifdef CORE_BLINK
          if (--(m->transitionTime) <= 0)
          {
            m->state = !m->state;
            m->transitionTime = (int)(log(DRandom(&s->rs)) *
    		(m->state ? s->k_d[sp] : s->k_b[sp]) );
          }
          if (!m->state) continue;		// next vesicle molec
  #endif
  #ifdef CORE_TRIPLET
          if (m->tripletRelaxTime) continue;	// next vesicle molec
  #endif

          /* Calculate Intensity */
          x += m->x;
          y += m->y;
          z += m->z;
  #ifdef CORE_INTENS
          if (abs((int)(x*sc_x)) <= numBins_x &&
              abs((int)(y*sc_x)) <= numBins_x &&
              abs((int)(z*sc_z)) <= numBins_z)
          {
            I = s->a[sp] * intens[
              (2*numBins_x+1)*(2*numBins_x+1)*(int)(numBins_z + z * sc_z) +
              (2*numBins_x+1) * (int)(numBins_x + y * sc_x) +
              (int)(numBins_x + x * sc_x) ];
  #else  // not INTENS
          if ( (r = x*x + y*y) < threshold && 
                fabs(z) < threshold_z)
          {
            I = s->a[sp] * exp(b*r + c*(z*z));
  #endif

  #ifdef CORE_TRIPLET
            if (abs((int)(x*cef_sc_x)) <= cef_numBins_x &&
                abs((int)(y*cef_sc_x)) <= cef_numBins_x &&
                abs((int)(z*cef_sc_z)) <= cef_numBins_z)
            {
              cef_value = cef[
                (2*cef_numBins_x+1) * (2*cef_numBins_x+1) * 
                	(int)(cef_numBins_z + z * cef_sc_z) +
                (2*cef_numBins_x+1) * (int)(cef_numBins_x + y * cef_sc_x) +
                (int)(cef_numBins_x + x * cef_sc_x) ];
            } else cef_value = 0.0;
  #endif // CORE_TRIPLET

  #ifdef CORE_BLEACH
            if ((m->photonTolerance -= I) < 0)
            {
              int k;
              I += m->photonTolerance;
              s->numBleached[sp]++;
	      // Delete molecule
              for (k=molec; k < v->numMolecs-1; k++)
                v->m[k] = v->m[k+1];
              v->numMolecs--;
              molec--;
              m--;
  #ifdef CORE_TRIPLET
            } else {
  #endif
  #endif // CORE_BLEACH
  #ifdef CORE_TRIPLET
            {
              double p = I*s->k_T[sp]/(I*s->k_T[sp]+s->k_S0[sp]);
              if (p > 1e-10 && DRandom(&s->rs) < p)
                if (m->tripletRelaxTime =
                      (int)(-log(DRandom(&s->rs)) / s->k_S0[sp]))
                  continue;
            }
  #endif

  #ifdef CORE_TRIPLET
            *binCount += I * cef_value;
  #else
            *binCount += I;
  #endif
          } // in threshold

        }  // each vesicle molecule

      } // each vesicle

    } // each step

    binCount++;
  } // each bin

  iData = (Int16*)o->data + pos;
  for (i=0; i < dur; i++)
#ifdef CORE_TRIPLET
    s->I += (double)(iData[i] = Poisson(s->phi*data[i]+s->bkg, &s->rs));
#else
    s->I += (double)(iData[i] = Poisson(data[i], &s->rs));
#endif

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
  
  l = PyList_New(c->numSpecies);
  if (!l) goto bail;
  for(i=0; i < c->numSpecies; i++)
  {
    o = PyInt_FromLong(c->numFree[i]);
    if (!o) goto bail;
    if (PyList_SetItem(l, i, o)) goto bail;	// steals reference
  }
  o = NULL;
  if (PyObject_SetAttrString(r, "freeCounts", l) == -1) goto bail;
  Py_DECREF(l);
  
  l = PyList_New(c->numVesicleTypes);
  if (!l) goto bail;
  for(i=0; i < c->numVesicleTypes; i++)
  {
    o = PyInt_FromLong(c->vCounts[i]);
    if (!o) goto bail;
    if (PyList_SetItem(l, i, o)) goto bail;	// steals reference
  }
  o = NULL;
  if (PyObject_SetAttrString(r, "vesicleCounts", l) == -1) goto bail;
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
  int i;
  if (!p) return;
  free(p->vNumMolecs);
  free(p->vR);
  free(p->vSigma);
  free(p->vX);
  free(p->vmX);
  free(p->dX);
  free(p->vdX);
  free(p->X);
  free(p->a);
#ifdef CORE_BLEACH
  free(p->tol);
#endif
#ifdef CORE_BLINK
  free(p->k_d);
  free(p->k_b);
#endif
#ifdef CORE_TRIPLET
  free(p->k_T);
  free(p->k_S0);
  if (p->cef) { Py_DECREF((PyObject*)p->cef); }
#endif
#ifdef CORE_INTENS
  if (p->intens) {  Py_DECREF((PyObject*)p->intens);  }
#endif
  free(p->m);
  for (i=0; i < p->numVesicles; i++)
    free(p->v[i].m);
  free(p->v);
  free(p->counts);
  free(p->vCounts);
  free(p->numFree);
#ifdef CORE_BLEACH
  free(p->numBleached);
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

// Box-Muller method for converting uniform to normal dist rand numbers
// Returns a normal random number with unit variance and zero mean
double Gaussian(RandomStream *rs)
{
  static int cached = 0;
  static double y1;
  double y2, r;

  if (cached)
  {
    cached = 0;
    return y1;
  }
  
  do
  {
    y1 = 2*DRandom(rs) - 1;
    y2 = 2*DRandom(rs) - 1;
    r = y1*y1 + y2*y2;
  } while (r >= 1 || r == 0);	// y1 and y2 must be in the unit circle

  r = sqrt(-2.0*log(r)/r);
  y1 *= r;
  cached = 1;
  return r*y2;
}
