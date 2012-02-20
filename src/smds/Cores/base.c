/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Basic Simulation Cores
 *    $Id: base.c,v 1.12 2009/09/08 21:27:28 mculbert Exp $
 *
 *  You must define CORE_NAME as an identifier.
 *  Optionally define CORE_DOC to a double-quoted doc-string.
 *  Defaults to 2D core.  For 3D define CORE_3D.
 *  For intens profile, define CORE_INTENS.
 *  For photobleaching, define CORE_BLEACH.
 *  For blinking, define CORE_BLINK.
 *  For alternate diffusion states, define CORE_STATES.
 *  For CCD detection, define CORE_CCD.
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
#if defined(CORE_CCD) && defined(CORE_3D)
#error "CCD detection may only be used with a 2D core."
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

#include "random.h"
#include <math.h>

typedef int Int32;
typedef short Int16;

typedef struct {
  Int32 x, y;
#ifdef CORE_3D
  Int32 z;
#endif
#ifdef CORE_STATES
  int species_state;
  unsigned int speciesTransitionTime;
#endif
#ifdef CORE_BLEACH
  double photonTolerance;
#endif
#ifdef CORE_BLINK
  int blink_state;
  int transitionTime;
#endif
#ifdef CORE_TRIPLET
  unsigned int tripletRelaxTime;
#endif
  int species;
} smds_molec;

typedef struct {
  PyObject_HEAD
  int ready;			// Initialization succeeded.
  int numSpecies, numMolecs;
  double *X;			// fraction of each species
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
  double *cef_sc_x;		// the laser intensity profile, not the
#ifdef CORE_3D			// entire molecule detection efficiency func
  int cef_numBins_z;
  double *cef_sc_z;
#endif
#endif // CORE_TRIPLET
  int steps;			// Steps per bin = binWidth / dT
  Int32 *size;			// Space dimension in units of dX
  				// dX = sqrt(4*D*dT)
  				// size = sqrt(numMolecs/density) / dX
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
#ifdef CORE_STATES
  double *k_sp;
  int *altSpecies;
#endif
#ifdef CORE_CCD
  int video_Nx, video_Ny;
  double *video_sc_x, *video_sc_y;
  int video_stepsPerFrame;
#endif
  /*** WORKING SPACE ***/
  RandomStream rs;
  smds_molec *m;
#ifdef CORE_CCD
  int frame, frame_stepsLeft;
  PyArrayObject *video;
#endif
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
static void smds_core_create_molec(smds_molec *m, smds_core *s, int edge, int sp);
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
  double dX, dT, size;
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
#ifdef CORE_TRIPLET
  int cef_binWidth_x;
#ifdef CORE_3D
  int cef_binWidth_z;
#endif
#endif
#ifdef CORE_CCD
  int video_sizeX, video_sizeY;
#endif
  
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
#ifdef CORE_CCD
  c->frame = 0;
#endif

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
  size = pow((double)(c->numMolecs)/PyFloat_AsDouble(o),
#ifdef CORE_3D
             1.0/3.0);
#else
             1.0/2.0);
#endif
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
  cef_binWidth_x = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_x");
  if (!o) goto bail;
  c->cef_numBins_x = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

#ifdef CORE_3D
  o = PyObject_GetAttrString(s, "binWidth_z");
  if (!o) goto bail;
  cef_binWidth_z = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);

  o = PyObject_GetAttrString(s, "numBins_z");
  if (!o) goto bail;
  c->cef_numBins_z = PyInt_AsLong(o);
  if (PyErr_Occurred()) goto bail;
  Py_DECREF(o);
#endif // CORE_3D
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

#ifdef CORE_CCD
  o = PyObject_GetAttrString(p, "video_framerate");
  if (!o) return -1;
  c->video_stepsPerFrame = (int)(PyFloat_AsDouble(o) * 1e-3 / dT);
  c->frame_stepsLeft = c->video_stepsPerFrame;
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "video_Nx");
  if (!o) return -1;
  c->video_Nx = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "video_Ny");
  if (!o) return -1;
  c->video_Ny = PyInt_AsLong(o);
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "video_sizeX");
  if (!o) return -1;
  video_sizeX = PyInt_AsLong(o);		// nm
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;

  o = PyObject_GetAttrString(p, "video_sizeY");
  if (!o) return -1;
  video_sizeY = PyInt_AsLong(o);		// nm
  Py_DECREF(o);
  if (PyErr_Occurred()) return -1;
#endif

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
#ifdef CORE_BLINK
  free(c->k_d); c->k_d = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->k_b); c->k_b = (double*)malloc(sizeof(double)*c->numSpecies);
#endif
#ifdef CORE_TRIPLET
  free(c->k_T);  c->k_T  = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->k_S0); c->k_S0 = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->cef_sc_x);  
  c->cef_sc_x = (double*)malloc(sizeof(double)*c->numSpecies);
#ifdef CORE_3D
  free(c->cef_sc_z);
  c->cef_sc_z = (double*)malloc(sizeof(double)*c->numSpecies);
#endif
#endif
#ifdef CORE_STATES
  free(c->k_sp);  c->k_sp = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->altSpecies); c->altSpecies = (int*)malloc(sizeof(int)*c->numSpecies);
#endif
#ifdef CORE_CCD
  free(c->video_sc_x);
  c->video_sc_x = (double*)malloc(sizeof(double)*c->numSpecies);
  free(c->video_sc_y);
  c->video_sc_y = (double*)malloc(sizeof(double)*c->numSpecies);
  if (c->video) { Py_DECREF(c->video); }
  {
    int dims[3] = { 1, c->video_Nx, c->video_Ny };

    PyObject *o_dur = PyObject_GetAttrString(p, "dur");
    if (!o_dur) return -1;
    dims[0] = ceil(PyFloat_AsDouble(o_dur) / (dT * c->video_stepsPerFrame));
    Py_DECREF(o_dur);
    if (PyErr_Occurred()) goto bail;

    c->video = (PyArrayObject*)PyArray_FromDims(3, dims, PyArray_DOUBLE);
  }
#endif
  free(c->m);  c->m = (smds_molec*)malloc(sizeof(smds_molec)*c->numMolecs);
  if (!(c->size && c->X && c->a && c->counts && c->m
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
#ifdef CORE_BLINK
        && c->k_d && c->k_b
#endif
#ifdef CORE_TRIPLET
        && c->k_T && c->k_S0 && c->cef_sc_x
#ifdef CORE_3D
        && c->cef_sc_z
#endif
#endif
#ifdef CORE_STATES
        && c->k_sp && c->altSpecies
#endif
#ifdef CORE_CCD
        && c->video_sc_x && c->video_sc_y && c->video
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

#ifdef CORE_STATES
    oo = PyObject_GetAttrString(s, "k_species");	// s^-1
    if (!oo) goto bail;
    c->k_sp[i] = -1/(dT*PyFloat_AsDouble(oo));
    if (PyErr_Occurred()) goto bail;
    Py_DECREF(oo);

    oo = PyObject_GetAttrString(s, "altSpecies");
    if (!oo) goto bail;
    c->altSpecies[i] = PyInt_AsLong(oo);
    if (PyErr_Occurred()) goto bail;
    if (c->altSpecies[i] < 0 || c->altSpecies[i] >= c->numSpecies)
    {
      PyErr_SetString(PyExc_ValueError, "Alternate species out of range");
      goto bail;
    }
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
    if (!(c->size[i] < sqrt((unsigned)(1 << 31))))
    {
      PyErr_SetString(PyExc_ValueError, "Simulation space is too large.");
      goto bail;
    }

#ifdef CORE_TRIPLET
    c->cef_sc_x[i] = 1.0 / cef_binWidth_x / dX * 1.0e3;
#ifdef CORE_3D
    c->cef_sc_z[i] = 1.0 / cef_binWidth_z / dX * 1.0e3;
#endif
#endif
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
#ifdef CORE_CCD
    c->video_sc_x[i] = 1.0 / video_sizeX / dX * 1.0e3;
    c->video_sc_y[i] = 1.0 / video_sizeY / dX * 1.0e3;
#endif

    c->counts[i] = 0;
    Py_DECREF(s);  s = NULL;
  }
  Py_DECREF(o); o = NULL;
  Py_DECREF(ss);  ss = NULL;

  {
    int sp, num;
    sp = -1;
    num = 0;
    for (i=0; i < c->numMolecs; i++)
    {
      while (sp > -2 && !num && sp < c->numSpecies)
      {
        sp++;
        num = (int)(c->X[sp] * c->numMolecs);
      }
      if (sp >= c->numSpecies) sp = -10;
      smds_core_create_molec(&(c->m[i]), c, 0, sp);
      num--;
    }
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

void smds_core_create_molec(smds_molec *m, smds_core *s, int edge, int sp)
{
  int r;
  double S;
   
  if (sp >= 0) m->species = sp;
  else
  {
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
  }
  s->counts[m->species]++;

#ifdef CORE_STATES
  m->species_state = DRandom(&s->rs) > 
  	1/s->k_sp[m->species] /
  	(1/s->k_sp[m->species] + 1/s->k_sp[s->altSpecies[m->species]]);
  sp = (m->species_state ? s->altSpecies[m->species] : m->species);
  m->speciesTransitionTime = (int)(log(DRandom(&s->rs))*s->k_sp[sp]);
#else
  sp = m->species;
#endif
                                                                                
#ifdef CORE_BLEACH
  /* Note that tolerance should be NEGATIVE */
  m->photonTolerance = log(DRandom(&s->rs)) * s->tol[m->species];
#endif

#ifdef CORE_BLINK
  m->blink_state = DRandom(&s->rs) > 
  	1/s->k_d[m->species]/(1/s->k_d[m->species]+1/s->k_b[m->species]);
  m->transitionTime = (int)(log(DRandom(&s->rs))*
  		(m->blink_state ? s->k_d[m->species] : s->k_b[m->species]) );
#endif

#ifdef CORE_TRIPLET
  m->tripletRelaxTime = 0;
#endif

  /* Place the new molecule at the edge */
  if (edge && DRandom(&s->rs) < 0.5)
  {
    m->x = (int)((DRandom(&s->rs) - 0.5) * s->size[sp] * 2.0);
    m->y = -s->size[sp];
#ifdef CORE_3D
    m->z = (int)((DRandom(&s->rs) - 0.5) * s->size[sp] * 2.0);
#endif
  }
  else
  {
    m->x = edge ? -s->size[sp] :
           (int)((DRandom(&s->rs) - 0.5) * s->size[sp] * 2.0);
    m->y = (int)((DRandom(&s->rs) - 0.5) * s->size[sp] * 2.0);
#ifdef CORE_3D
    m->z = (int)((DRandom(&s->rs) - 0.5) * s->size[sp] * 2.0);
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
  int dur, pos, molec, i;
  smds_molec *m;
  
  int x, y, size;
#ifdef CORE_3D
  int z;
#endif
  double a;
#ifdef CORE_INTENS
  double *intens = (double*)(s->intens->data);
  double sc_x;
  int numBins_x = s->numBins_x, intensMax_x = 2*s->numBins_x+1;
#ifdef CORE_3D
  double sc_z;
  int numBins_z = s->numBins_z, intensMax_z = 2*s->numBins_z+1;
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
  int species;
#ifdef CORE_STATES
  unsigned int speciesTransitionTime;
#endif
#ifdef CORE_BLINK
  int blink_state, transitionTime;
#endif
#ifdef CORE_TRIPLET
  double k_T, k_S0;
  unsigned int tripletRelaxTime;
  double cef_value, *cef = (double*)(s->cef->data);
  double cef_sc_x;
  int cef_numBins_x = s->cef_numBins_x, cefMax_x = 2*s->cef_numBins_x+1;
#ifdef CORE_3D
  double cef_sc_z;
  int cef_numBins_z = s->cef_numBins_z, cefMax_z = 2*s->cef_numBins_z+1;
#endif
#endif
  int t, dir, rndCounter = 0;
  Int32 rnd = 0;
  double *binCount;
  double I;
#ifdef CORE_CCD
  int frame, frame_stepsLeft, X, Y;
  int fX = s->video->strides[1]/sizeof(double),
      fY = s->video->strides[2]/sizeof(double);
  int video_Nx = s->video_Nx, video_Ny = s->video_Ny;
  double *video, *video_loc, video_sc_x, video_sc_y;
#endif

  if (!PyArg_ParseTuple(args, "iO!i", &dur, &PyArray_Type, &o, &pos))
    return NULL;

  s->dur += dur;
  data = (double*)malloc(sizeof(double)*dur);
  if (!data)
  {
    PyErr_SetString(PyExc_MemoryError, "Unable to allocate working space.");
    return NULL;
  }

  void loadSpecies()
  {
#ifdef CORE_STATES
    species = (m->species_state ? s->altSpecies[m->species] : m->species);
#else
    species = m->species;
#endif
    size = s->size[species];
    a = s->a[species];
#ifdef CORE_INTENS
    sc_x = s->sc_x[species];
#ifdef CORE_3D
    sc_z = s->sc_z[species];
#endif
#else
    threshold = s->threshold[species];
    b = s->b[species];
#ifdef CORE_3D
    threshold_z = s->threshold_z[species];
    c = s->c[species];
#endif
#endif
#ifdef CORE_CCD
    video_sc_x = s->video_sc_x[species];
    video_sc_y = s->video_sc_y[species];
#endif
#ifdef CORE_TRIPLET
    k_T = s->k_T[species];
    k_S0 = s->k_S0[species];
    cef_sc_x = s->cef_sc_x[species];
#ifdef CORE_3D
    cef_sc_z = s->cef_sc_z[species];
#endif
#endif
  }

  void loadMolecule()
  {
    x = m->x;
    y = m->y;
#ifdef CORE_3D
    z = m->z;
#endif
#ifdef CORE_STATES
    speciesTransitionTime = m->speciesTransitionTime;
#endif
#ifdef CORE_BLINK
    blink_state = m->blink_state;
    transitionTime = m->transitionTime;
#endif
#ifdef CORE_TRIPLET
    tripletRelaxTime = m->tripletRelaxTime;
#endif
  }

  void saveMolecule()
  {
    m->x = x;
    m->y = y;
#ifdef CORE_3D
    m->z = z;
#endif
#ifdef CORE_STATES
    m->speciesTransitionTime = speciesTransitionTime;
#endif
#ifdef CORE_BLINK
    m->blink_state = blink_state;
    m->transitionTime = transitionTime;
#endif
#ifdef CORE_TRIPLET
    m->tripletRelaxTime = tripletRelaxTime;
#endif
  }

  Py_BEGIN_ALLOW_THREADS
#if !defined(CORE_TRIPLET)
  if (s->bkg > 0.0)
    for (i=0; i < dur; i++) data[i] = s->bkg;
  else
#endif
    for (i=0; i < dur; i++) data[i] = 0.0;

  for (molec=0; molec < s->numMolecs; molec++)
  {
    m = &(s->m[molec]);
    loadSpecies();
    loadMolecule();
#ifdef CORE_CCD
    frame = s->frame;
    frame_stepsLeft = s->frame_stepsLeft;
    video = (double*)s->video->data + s->video->strides[0] * frame
                                          / sizeof(double);
#endif

    binCount = data;
    for (i=dur; i; i--)
    {
      for (t=s->steps; t; t--)
      {
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
          case 0: x++;  if (x >= size) { x = -size;  }  break;
          case 1: x--;  if (x < -size) { x = size-1; }  break;
          case 2: y++;  if (y >= size) { y = -size;  }  break;
          case 3: y--;  if (y < -size) { y = size-1; }  break;
#ifdef CORE_3D
          case 4: z++;  if (z >= size) { z = -size;  }  break;
          case 5: z--;  if (z < -size) { z = size-1; }  break;
#endif
        }
      
#ifdef CORE_CCD
        X = x*video_sc_x + video_Nx/2;
        Y = y*video_sc_y + video_Ny/2;
        if (X >= 0 && X < video_Nx && Y >= 0 && Y < video_Ny)
          video_loc = video + fX*X + fY*Y;
        else video_loc = NULL;
        if (!(--frame_stepsLeft))
        {
          frame++;
          frame_stepsLeft = s->video_stepsPerFrame;
          video += s->video->strides[0]/sizeof(double);
        }
#endif

#ifdef CORE_TRIPLET
        if (tripletRelaxTime) tripletRelaxTime--;
#endif
#ifdef CORE_STATES
        if (--speciesTransitionTime <= 0)
        {
          double scale;
          int fromSpecies = species;
          m->species_state = !m->species_state;
          loadSpecies();
          speciesTransitionTime =
                (int)(log(DRandom(&s->rs)) * s->k_sp[species]);
          scale = (double)s->size[species] / (double)s->size[fromSpecies];
          x = (int)(x*scale);
          y = (int)(y*scale);
#ifdef CORE_3D
          z = (int)(z*scale);
#endif
        }
#endif
#ifdef CORE_BLINK
        if (--transitionTime <= 0)
        {
          blink_state = !blink_state;
          transitionTime = (int)(log(DRandom(&s->rs)) *
  		(blink_state ? s->k_d[m->species] : s->k_b[m->species]) );
        }
        if (!blink_state) continue;
#endif
#ifdef CORE_TRIPLET
        if (tripletRelaxTime) continue;
#endif

        /* Calculate Intensity */
#ifdef CORE_INTENS
#ifdef CORE_3D
        if ((int)(numBins_x+x*sc_x) < intensMax_x &&
            (int)(numBins_x+x*sc_x) >= 0 &&
            (int)(numBins_x+y*sc_x) < intensMax_x &&
            (int)(numBins_x+y*sc_x) >= 0 &&
            (int)(numBins_z+z*sc_z) < intensMax_z &&
            (int)(numBins_z+z*sc_z) >= 0)
        {
          I = a * intens[
            intensMax_x * intensMax_x * (int)(numBins_z + z * sc_z) +
            intensMax_x * (int)(numBins_x + y * sc_x) +
            (int)(numBins_x + x * sc_x) ];
#else
        if ((int)(numBins_x+x*sc_x) < intensMax_x &&
            (int)(numBins_x+x*sc_x) >= 0 &&
            (int)(numBins_x+y*sc_x) < intensMax_x &&
            (int)(numBins_x+y*sc_x) >= 0)
        {
          I = a * intens[
            intensMax_x * (int)(numBins_x + y * sc_x) +
            (int)(numBins_x + x * sc_x)  ];
#endif
#else  // not INTENS
#ifdef CORE_3D
        if ( (r = x*x + y*y) < threshold && abs(z) < threshold_z)
        {
          I = a * exp(b*(double)r + c*(double)(z*z));
#else
        if ( (r = x*x + y*y) < threshold)
        {
          I = a * exp(b * (double)r);
#endif
#endif

#ifdef CORE_TRIPLET
#ifdef CORE_3D
          if ((int)(cef_numBins_x+x*cef_sc_x) < cefMax_x &&
              (int)(cef_numBins_x+x*cef_sc_x) >= 0 &&
              (int)(cef_numBins_x+y*cef_sc_x) < cefMax_x &&
              (int)(cef_numBins_x+y*cef_sc_x) >= 0 &&
              (int)(cef_numBins_z+z*cef_sc_z) < cefMax_z &&
              (int)(cef_numBins_z+z*cef_sc_z) >= 0)
          {
            cef_value = cef[
              cefMax_x * cefMax_x * (int)(cef_numBins_z + z * cef_sc_z) +
                         cefMax_x * (int)(cef_numBins_x + y * cef_sc_x) +
                                    (int)(cef_numBins_x + x * cef_sc_x) ];
          } else cef_value = 0.0;
#else
          if ((int)(cef_numBins_x+x*cef_sc_x) < cefMax_x &&
              (int)(cef_numBins_x+x*cef_sc_x) >= 0 &&
              (int)(cef_numBins_x+y*cef_sc_x) < cefMax_x &&
              (int)(cef_numBins_x+y*cef_sc_x) >= 0)
          {
            cef_value = cef[
              cefMax_x * (int)(cef_numBins_x + y * cef_sc_x) +
              (int)(cef_numBins_x + x * cef_sc_x)  ];
          } else cef_value = 0.0;
#endif
#endif // CORE_TRIPLET

#ifdef CORE_BLEACH
          if ((m->photonTolerance -= I) < 0)
          {
            I += m->photonTolerance;
            s->numBleached[m->species]++;
            smds_core_create_molec(m, s, 1, -1);
            loadSpecies();
            loadMolecule();
#ifdef CORE_TRIPLET
          } else
#else // not CORE_TRIPLET
          }
#endif
#endif // CORE_BLEACH
#ifdef CORE_TRIPLET
          {
            double p = I*k_T/(I*k_T+k_S0);
            if (p > 1e-10 && DRandom(&s->rs) < p)
              if (tripletRelaxTime = (int)(-log(DRandom(&s->rs)) / k_S0))
                continue;
          }
#endif

#ifdef CORE_TRIPLET
          I *= cef_value;
#endif
          *binCount += I;
#ifdef CORE_CCD
          if (video_loc) *video_loc += I;
#endif            

        } // in threshold

      } // each step

      binCount++;
    } // each bin

    saveMolecule();
  }  // each molecule

#ifdef CORE_CCD
  s->frame = frame;
  s->frame_stepsLeft = frame_stepsLeft;
#endif
  
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
#ifdef CORE_CCD
  PyArrayObject *a;
  int x, y, sI, dI, sX, dX, sY, dY;
  double *data;
  short *dest;
#endif
  
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

#ifdef CORE_CCD
  l = PyArray_FromDims(3, c->video->dimensions, PyArray_SHORT);
  if (!l) goto bail;
  a = (PyArrayObject*)l;
  sI = c->video->strides[0]/sizeof(double);
  sX = c->video->strides[1]/sizeof(double);
  sY = c->video->strides[2]/sizeof(double);
  data = (double*)c->video->data;
  dI = a->strides[0]/sizeof(short);
  dX = a->strides[1]/sizeof(short);
  dY = a->strides[2]/sizeof(short);
  dest = (short*)a->data;
  for (i=0; i < c->video->dimensions[0]; i++)
    for (x=0; x < c->video_Nx; x++)
      for (y=0; y < c->video_Ny; y++)
        dest[dI*i + dX*x + dY*y] = Poisson(data[sI*i + sX*x + sY*y], &c->rs);
  if (PyObject_SetAttrString(r, "video", l) == -1) goto bail;
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
  if (!p) return;
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
  free(p->cef_sc_x);
#ifdef CORE_3D
  free(p->cef_sc_z);
#endif
#endif
  free(p->size);
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
  free(p->m);
  free(p->counts);
#ifdef CORE_BLEACH
  free(p->numBleached);
#endif
#ifdef CORE_CCD
  free(p->video_sc_x);
  free(p->video_sc_y);
  if (p->video) {  Py_DECREF((PyObject*)p->video);  }
#endif
#ifdef CORE_STATES
  free(p->k_sp);
  free(p->altSpecies);
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
    if (!PyArg_ParseTuple(args, "i", &seed))
      return NULL;
    seedRnd(seed, NULL);
  }

  Py_RETURN_NONE;
}
