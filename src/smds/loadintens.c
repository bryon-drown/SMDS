/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    C routines to load and write intensity profile maps
 *
 *    Intensity files must be of the format:
 *         binWidth_x binWidth_z numBins_x numBins_z data
 *    where number of datapoints must be exactly
 *         (2*numBins_z+1) * (2*numBins_x+1)^2
 *    thus, data[x,y,z] goes from [-numBins, numBins]."""
 *
 *    aHL Potential Profiles must be of the format:
 *         binWidth_A_z binWidth_A_r numBins_A_z numBins_A_r data_A
 *         binWidth_B_z binWidth_B_r numBins_B_z numBins_B_z data_B
 *    where A refers to everywhere outside the channel (from the membrane,
 *    z=0, up) and B refers to in the channel (from the top of the cap,
 *    z=capHeight, down).  Note that the two regions overlap, and that
 *    the direction of data_A is in *increasing* z and the direction of
 *    data_B is in *decreasing* z.  Data extend radially from r=0 out.
 *    The two data segments have two tables of numBins_?_z rows and
 *    numBins_?_r columns.  The first table indicates the electrophoretic 
 *    flow in z and the second in r.
 *
 *******************************************/

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL smds_array_symbol
#include <numpy/arrayobject.h>
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE { Py_INCREF(Py_None); return Py_None; }
#endif

static PyObject * load(PyObject *self, PyObject *args);
static PyObject * writeIntens(PyObject *self, PyObject *args);
static PyObject * loadPotProfile(PyObject *self, PyObject *args);
static PyObject * writePotProfile(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"_loadintens", load, METH_VARARGS,
   "_loadintens(filename) -> (bwx, bwz, Nx, Nz, data)" },
  {"_writeintens", writeIntens, METH_VARARGS,
   "_writeintens(filename, bwx, bwz, Nx, Nz, data)" },
  {"_loadpot", loadPotProfile, METH_VARARGS,
   "_loadpot(filename) -> (bwAz, bwAr, NAz, NAr, dataA, bwBz, bwBr, NBz, NBr, dataB)" },
  {"_writepot", writePotProfile, METH_VARARGS,
   "_writepot(filename, bwAz, bwAr, NAz, NAr, dataA, bwBz, bwBr, NBz, NBr, dataB)" },
  {NULL}	// Sentinel
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initloadintens(void)
{
  PyObject *m;
  m = Py_InitModule3("loadintens", module_methods,
  			"Read/write SMDS intensity profile files.");
  import_array();
}

// returns (bwx, bwz, Nx, Nz, data)
static PyObject * load(PyObject *self, PyObject *args)
{
  FILE *f;
  PyArrayObject *a;
  PyObject *o;
  char *filename;
  double *data, max = 0.0;
  int binWidth_x, binWidth_z, numBins_x, numBins_z;
  int i, num;
  npy_intp dim;
  
  if (!PyArg_Parse(args, "(s)", &filename)) return NULL;
  
  f = fopen(filename, "r");
  if (!f) return PyErr_SetFromErrno(PyExc_IOError);
                                                                                
  if (4 != fscanf(f, "%d %d %d %d", &binWidth_x,
                        &binWidth_z, &numBins_x,  &numBins_z))
        {
          fclose(f);
          return PyErr_Format(PyExc_ValueError, 
              "Error processing Intensity Profile file %s.  "
                        "Check header.", filename);
        }
                                                                                
  dim = num = (2*numBins_x+1)*(2*numBins_x+1) * (2*numBins_z+1);
  a = (PyArrayObject*)PyArray_SimpleNew(1, &dim, PyArray_DOUBLE);
  if (!a)
        {
          fclose(f);
          return NULL;
        }
  data = (double*)(a->data);
  for (i=0; i < num; i++)
  {
    if (1 != fscanf(f, "%lf", &data[i]))
        {
          fclose(f);
          Py_DECREF(a);
          return PyErr_Format(PyExc_ValueError,
          	"Error reading Intensity Profile %s at point %d.",
                        filename, i);
        }
    if (data[i] > max) max = data[i];
    if (data[i] < 0.0) data[i] = 0.0;
  }
  fclose(f);
                                                                                
  if (max) for (i=0; i < num; i++)
    data[i] /= max;

  o = Py_BuildValue("(llllO)", binWidth_x, binWidth_z, numBins_x, 
  				numBins_z, a);
  Py_DECREF(a);
  return o;
}

// writeintens(filename, bwx, bwz, Nx, Nz, data) -> None
static PyObject * writeIntens(PyObject *self, PyObject *args)
{
  FILE *f;
  PyArrayObject *a;
  char *filename;
  double *data;
  int binWidth_x, binWidth_z, numBins_x, numBins_z;
  int i, x, y, z, num;
  
  if (!PyArg_Parse(args, "(siiiiO!)", &filename, &binWidth_x, &binWidth_z,
                   &numBins_x, &numBins_z, &PyArray_Type, &a)) return NULL;
  num = (2*numBins_x+1)*(2*numBins_x+1) * (2*numBins_z+1);
  if (a->nd != 1 || num != a->dimensions[0])
  {
    return PyErr_Format(PyExc_ValueError,
                        "Data not the correct shape (%d).", num);
  }
            // a's ref count is increased
  a = (PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)a,
  			PyArray_DOUBLE, 1, 1);
  if (!a) return NULL;
  data = (double*)(a->data);
  
  f = fopen(filename, "w");
  if (!f)
  {
    Py_DECREF(a); 
    return PyErr_SetFromErrno(PyExc_IOError);
  }
                                                                                
  if (fprintf(f, "%d %d %d %d\n", binWidth_x, binWidth_z, numBins_x,
              numBins_z) < 0)
        {
          fclose(f);
          Py_DECREF(a);
          return PyErr_SetFromErrno(PyExc_IOError);
        }
                                                                                
  numBins_x = numBins_x*2+1;
  numBins_z = numBins_z*2+1;
  i = 0;
  for (z=0; z < numBins_z; z++)
  {
    for (x=0; x < numBins_x; x++)
    {
      for (y=0; y < numBins_x; y++)
      {
        if (fprintf(f, "%g ", data[i++]) < 0)
        {
          fclose(f);
          Py_DECREF(a);
          return PyErr_SetFromErrno(PyExc_IOError);
        }
      } // y
      fprintf(f, "\n");
    } // x
    fprintf(f, "\n");
  } // z
  fclose(f);
                                                                                
  Py_DECREF(a);
  Py_RETURN_NONE;
}

// returns (bwAz, bwAr, NAz, NAr, dataA, bwBz, bwBr, NBz, NBr, dataB)
static PyObject * loadPotProfile(PyObject *self, PyObject *args)
{
  FILE *f;
  PyArrayObject *a, *b;
  PyObject *o;
  char *filename;
  double *data;
  double bwAz, bwAr, bwBz, bwBr;
  int NAz, NAr, NBz, NBr;
  int i, num;
  npy_intp dim;
  
  if (!PyArg_Parse(args, "(s)", &filename)) return NULL;
  
  f = fopen(filename, "r");
  if (!f) return PyErr_SetFromErrno(PyExc_IOError);
                                                                                
  if (4 != fscanf(f, "%lf %lf %d %d", &bwAz, &bwAr, &NAz, &NAr))
        {
          fclose(f);
          return PyErr_Format(PyExc_ValueError, 
              "Error processing Potential Profile file %s.  "
                        "Check header A.", filename);
        }
                                                                                
  dim = num = 2 * NAz * NAr;
  a = (PyArrayObject*)PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
  if (!a)
        {
          fclose(f);
          return NULL;
        }
  data = (double*)(a->data);
  for (i=0; i < num; i++)
  {
    if (1 != fscanf(f, "%lf", &data[i]))
        {
          fclose(f);
          Py_DECREF(a);
          return PyErr_Format(PyExc_ValueError,
          	"Error reading Potential Profile %s at point %d, Zone A.",
                        filename, i);
        }
  }
                                                                                
  if (4 != fscanf(f, "%lf %lf %d %d", &bwBz, &bwBr, &NBz, &NBr))
        {
          fclose(f);
          Py_DECREF(a);
          return PyErr_Format(PyExc_ValueError, 
              "Error processing Potential Profile file %s.  "
                        "Check header B.", filename);
        }
                                                                                
  dim = num = 2 * NBz * NBr;
  b = (PyArrayObject*)PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
  if (!b)
        {
          fclose(f);
          Py_DECREF(a);
          return NULL;
        }
  data = (double*)(b->data);
  for (i=0; i < num; i++)
  {
    if (1 != fscanf(f, "%lf", &data[i]))
        {
          fclose(f);
          Py_DECREF(a);
          Py_DECREF(b);
          return PyErr_Format(PyExc_ValueError,
          	"Error reading Potential Profile %s at point %d, Zone B.",
                        filename, i);
        }
  }
  fclose(f);
                                                                                
  o = Py_BuildValue("(ddllOddllO)", bwAz, bwAr, NAz, NAr, a,
                                    bwBz, bwBr, NBz, NBr, b);
  Py_DECREF(a);
  Py_DECREF(b);
  return o;
}

// writeintens(filename, bwAz, bwAr, NAz, NAr, dataA, bwBz, bwBr, NBz, NBr, dataB) -> None
static PyObject * writePotProfile(PyObject *self, PyObject *args)
{
  FILE *f;
  PyArrayObject *a, *b;
  char *filename;
  double *data;
  double bwAz, bwAr, bwBz, bwBr;
  int NAz, NAr, NBz, NBr;
  int i, numA, numB;
  
  if (!PyArg_Parse(args, "(sddiiO!ddiiO!)", &filename,
                   &bwAz, &bwAr, &NAz, &NAr, &PyArray_Type, &a,
                   &bwBz, &bwBr, &NBz, &NBr, &PyArray_Type, &b)) return NULL;
  numA = 2 * NAz * NAr;
  if (a->nd != 1 || numA != a->dimensions[0])
  {
    return PyErr_Format(PyExc_ValueError,
                        "Zone A Data not the correct shape (%d).", numA);
  }
  numB = 2 * NBz * NBr;
  if (b->nd != 1 || numB != b->dimensions[0])
  {
    return PyErr_Format(PyExc_ValueError,
                        "Zone B Data not the correct shape (%d).", numB);
  }

            // ref counts increased
  a = (PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)a,
  			PyArray_DOUBLE, 1, 1);
  if (!a) return NULL;
  b = (PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)b,
  			PyArray_DOUBLE, 1, 1);
  if (!b) { Py_DECREF(a); return NULL; }
  
  f = fopen(filename, "w");
  if (!f)
  {
    Py_DECREF(a); 
    Py_DECREF(b); 
    return PyErr_SetFromErrno(PyExc_IOError);
  }

  // Zone A
  if (fprintf(f, "%g %g %d %d\n", bwAz, bwAr, NAz, NAr) < 0)
        {
          fclose(f);
          Py_DECREF(a);
          Py_DECREF(b);
          return PyErr_SetFromErrno(PyExc_IOError);
        }

  data = (double*)(a->data);
  for (i=0; i < numA; i++)
    if (fprintf(f, "%g ", data[i]) < 0)
        {
          fclose(f);
          Py_DECREF(a);
          Py_DECREF(b);
          return PyErr_SetFromErrno(PyExc_IOError);
        }
  fprintf(f, "\n");

  // Zone B
  if (fprintf(f, "%g %g %d %d\n", bwBz, bwBr, NBz, NBr) < 0)
        {
          fclose(f);
          Py_DECREF(a);
          Py_DECREF(b);
          return PyErr_SetFromErrno(PyExc_IOError);
        }

  data = (double*)(b->data);
  for (i=0; i < numB; i++)
    if (fprintf(f, "%g ", data[i]) < 0)
        {
          fclose(f);
          Py_DECREF(a);
          Py_DECREF(b);
          return PyErr_SetFromErrno(PyExc_IOError);
        }
  fprintf(f, "\n");

  fclose(f);

  Py_DECREF(a);
  Py_DECREF(b);
  Py_RETURN_NONE;
}
