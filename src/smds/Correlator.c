/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Distributable Multi-tau Correlator
 *    $Id: Correlator.c,v 1.6 2009/06/11 18:28:52 mculbert Exp $
 *
 *    This multi-tau correlator is based on the design presented in:
 *      Thorsten Wohland, Rudolf Rigler, and Horst Vogel. ``The Standard
 *      Deviation in Fluorescence Correlation Spectroscopy.''  _Biophysical
 *      Journal,_ 80 p 2987 (2001).
 *    The combination algorithm is detailed in:
 *      Michael J. Culbertson and Daniel L. Burden. ``A Distributed
 *      Algorithm for Multi-tau Autocorrelation.'' _Review of Scientific
 *      Instruments,_ (2007).
 *    Implemented by Michael Culbertson, Wheaton College, 2004.
 *      Python module, 2005.
 *      Updated, 2007.
 *
 *    It is important to note that if the data from two correlators that
 *    are to be merged are misaligned, a slight error is introduced (half a
 *    bin is dropped) in the delays affected by the misalignment.  This
 *    shouldn't be a problem as long as the total record to be correlator
 *    is much greater than the number of merged correlators.
 *
 *******************************************/

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL smds_array_symbol
#include <numpy/arrayobject.h>
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE { Py_INCREF(Py_None); return Py_None; }
#endif

typedef int Int32;
typedef unsigned int UInt32;
typedef short Int16;

	// Groups are numbered [1, 16)
	// Group 1 has 16 elements numberd [-8, 8)
	// Groups [2, 16) have 8 elements numbered [0, 8)
	// Dynamic rage is 2^(K+3)
#define K 21			/* maximum group number */
typedef struct {
  PyObject_HEAD
  double shift[8*(K+1)], start[8*(K+1)], G[8*(K+1)], M_del[8*(K+1)], 
         M_dir[(K+1)], zero[(K+1)], acc[(K+1)];
  double G0, sum, bw;
  char pos[(K+1)];		// Starting positions in the round-robin arrays
  char fill[(K+1)];	// Bins of start left to fill
  UInt32 dur;
  Int32 parity;		// k'th bit true if (k+1)[0] is half-filled
} smds_Correlator;

// takes one arg: double binWidth(ms)
static int smds_Correlator_new(PyObject *self, 
				PyObject *args, PyObject *kwds);
static PyObject * smds_Correlator_data(PyObject *self, PyObject *o);
static PyObject * smds_Correlator_append(PyObject *self, PyObject *args);
static PyObject * smds_Correlator_output(PyObject *self);
static PyObject * smds_Correlator_getstate(PyObject *self);
static PyObject * smds_Correlator_setstate(PyObject *self, PyObject *args);
static PyObject * smds_Correlator_reduce(PyObject *self);

static PyTypeObject smds_CorrelatorType = {
  PyObject_HEAD_INIT(NULL)
  0,				// ob_size
  "smds.Analysis.Correlator",	// tp_name
  sizeof(smds_Correlator),	// tp_basicsize
  0,				// tp_itemsize
  0,				// tp_dealloc
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

static PyObject *ACFType = NULL;
static PyObject *RTRType = NULL;

static PyMethodDef smds_Correlator_methods[] = {
  {"append", (PyCFunction)smds_Correlator_append, METH_VARARGS,
   "Append more data to the Correlator.\n  "
   "Argument is either of type Correlator or of type RTR."},
  {"getACF", (PyCFunction)smds_Correlator_output, METH_NOARGS,
   "Returns the smds.Analysis.ACF corresponding the the correlated data."},
  {"__getstate__", (PyCFunction)smds_Correlator_getstate, METH_NOARGS, ""},
  {"__setstate__", (PyCFunction)smds_Correlator_setstate, METH_VARARGS, ""},
  {"__reduce__", (PyCFunction)smds_Correlator_reduce, METH_NOARGS, ""},
  {NULL}	// Sentinel
};

static PyMethodDef module_methods[] = {
  {NULL}	// Sentinel
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initCorrelator(void)
{
  PyObject *m;
  smds_CorrelatorType.tp_new = PyType_GenericNew;
  smds_CorrelatorType.tp_init = smds_Correlator_new;
  smds_CorrelatorType.tp_methods = smds_Correlator_methods;
  if (PyType_Ready(&smds_CorrelatorType) < 0)
    return;
  m = Py_InitModule3("Correlator", module_methods,
  			"A Distributable Multi-Tau Auto-Correlator");
  Py_INCREF(&smds_CorrelatorType);
  PyModule_AddObject(m, "Correlator", (PyObject*)&smds_CorrelatorType);
  import_array();

  m = PyImport_ImportModule("smds.Analysis");
  if (m)
  {
    ACFType = PyObject_GetAttrString(m, "ACF");
    Py_DECREF(m);
  }
  m = PyImport_ImportModule("smds.Task");
  if (m)
  {
    RTRType = PyObject_GetAttrString(m, "RTR");
    Py_DECREF(m);
  }
}

// k is on [1, (K+1)); If (k==1), i is on [-8, 8)
#define POS(k, i) ((k)*8+(i))

int smds_Correlator_new(PyObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = { "binWidth", NULL };
  smds_Correlator *p = (smds_Correlator*)self;
  double b;
  int i;
  
  // Get binWidth argument (ms)
  if (! PyArg_ParseTupleAndKeywords(args, kwds, "d", kwlist, &b))
    return -1;

  p->bw = b;
  for (i=0; i < 8*(K+1); i++)
    p->shift[i] = p->start[i] = p->G[i] = p->M_del[i] = 0.0;
  for (i=0; i < (K+1); i++)
  {
    p->M_dir[i] = p->zero[i] = p->acc[i] = 0.0;
    p->pos[i] = 0;
    p->fill[i] = 8;
  }
  p->fill[0] = 0;
  p->fill[1] = 16;
  p->pos[1] = -8;
  p->dur = 0;
  p->parity = 0;
  p->G0 = p->sum = 0.0;

  return 0;
}

PyObject * smds_Correlator_output(PyObject *self)
{
  smds_Correlator *p = (smds_Correlator *)self;
  PyObject *r = NULL, *o = NULL;
  PyArrayObject *Gs = NULL, *ts = NULL;
  double norm, *G, *t, bin = 0.0;
  int num;
  int i, k, l = 0;			// set l = 1 if g(0) is included

  if (ACFType) r = PyInstance_New(ACFType, NULL, NULL);
  else PyErr_SetString(PyExc_SystemError, "Unable to create ACF object.");
  if (!r)  return NULL;
  
  num = 8*(K+1);				// 129 if g(0) is included
  for (k=1; k < (K+1); k++)
    num -= p->fill[k];
  if (num <= 0)
  {
    Py_DECREF(r);
    Py_RETURN_NONE;
  }
  
  Gs = (PyArrayObject*)PyArray_FromDims(1, &num, PyArray_DOUBLE);
  if (!Gs) goto bail;
  ts = (PyArrayObject*)PyArray_FromDims(1, &num, PyArray_DOUBLE);
  if (!ts) goto bail;
  if (PyObject_SetAttrString(r, "G", (PyObject*)Gs) == -1 ||
      PyObject_SetAttrString(r, "t", (PyObject*)ts) == -1)
    goto bail;
  t = (double*)(ts->data);
  G = (double*)(Gs->data);

/*
  t[0] = 0.0;
  if (!p->sum)
    G[0] = 0.0;
  else
    G[0] = p->G0 / (p->sum*p->sum) * (double)(p->dur);
*/
  
  for (k=1; k < (K+1); k++)
  {
    norm = (double)(p->dur/(1<<(k-1))) / p->M_dir[k];
    for (i=(k==1?-8:0); i < 8-p->fill[k]; i++)
    {
      bin = t[l] = bin + (double)(1 << (k-1)) * p->bw * 1e-3;
      G[l] = norm * p->G[POS(k,i)] / p->M_del[POS(k,i)];
      if (finite(G[l])) l++;
    }
  }
  
  if (l != num)
    Gs->dimensions[0] = ts->dimensions[0] = l;
  
  Py_DECREF(Gs);
  Py_DECREF(ts);
  return r;

bail:
  if (o)  { Py_DECREF(o); }
  if (r)  { Py_DECREF(r); }
  if (Gs) { Py_DECREF(Gs); }
  if (ts) { Py_DECREF(ts); }
  return NULL;
}

// argument is of type smds.Task.RTR
PyObject * smds_Correlator_data (PyObject *self, PyObject *o)
{
  smds_Correlator *c = (smds_Correlator*)self;
//  PyObject *o = NULL;
  PyArrayObject *a = NULL;
  short *data = NULL;
  int k, i, n, rr, num;

/*   *** This func is called only from append() when the arg is a RTR.
  if (!PyArg_ParseTuple(args, "O!", RTRType, &o))  return NULL;
  // o is borrowed
  if (!PyInstance_Check(o))
  {
    PyErr_SetString(PyExc_TypeError, "Argument must be of type smds.RTR.");
    return NULL;
  }
*/
  if (!(o = PyObject_GetAttrString(o, "data"))) return NULL;
  // o is ours
  a = (PyArrayObject*)PyArray_ContiguousFromObject(o, PyArray_SHORT, 1, 1);
  Py_DECREF(o);  o = NULL;
  if (!a) return NULL;
  
  data = (short*)a->data;
  num = a->dimensions[0];
  if (num <= 0) { Py_DECREF(a);  Py_RETURN_NONE; }

  for (n=0; n < num; n++)
  {
    c->dur++;
    c->zero[1] = (double)data[n];
    c->sum += c->zero[1];
    c->G0 += c->zero[1] * c->zero[1];
    for (k=1; k < (K+1); k++)
    {
      // New datum has entered zero, trickle up
      if (k != K)
        c->zero[k+1] += c->zero[k];
      if (!c->fill[k-1] && c->fill[k])
        c->start[POS(k, 8-(c->fill[k]--))] = c->zero[k];
      
      // Do multiplications
      c->M_dir[k] += c->zero[k];
      i = (k==1?0:k*8);
      rr = c->pos[k];
      do
      {
        c->M_del[i] += c->shift[POS(k,rr)];
        c->G[i] += c->zero[k] * c->shift[POS(k,rr)];
        i++;
        if (++rr == 8) rr = (k==1?-8:0);
      } while (rr != c->pos[k]);
      
      // Shift
      rr = c->pos[k]-1;		// oldest position becomes newest
      if (rr < (k==1?-8:0)) rr = 7;
      if (k != K)		// Shift oldest position to next group's acc
        c->acc[k+1] += c->shift[POS(k,rr)];
				// Shift this group's acc into newest pos
      c->shift[POS(k,rr)] = (k==1?c->zero[k]:c->acc[k]);
      c->pos[k] = rr;		// save newest pos
      c->parity ^= 1 << k;

      // Recurse
      c->acc[k] = c->zero[k] = 0.0;
      if (c->parity & (1 << k)) break;
    } // each k
  } // each datum

  Py_DECREF(a);
  Py_RETURN_NONE;
} // smds_Correlator_data()

// argument is of type smds_Correlator or smds.Task.RTR
PyObject * smds_Correlator_append(PyObject *self, PyObject *args)
{
  smds_Correlator *a = (smds_Correlator*)self, *b = NULL;
  PyObject *o = NULL, *c = NULL;
  double swap[(K+1)], *t1 = swap, *t2 = swap+8;
  int k, i, rr, phi, rem = 0, used = 0;	   // k, sigma, tau, phi,  rho, xsi
  
  if (!PyArg_ParseTuple(args, "O", &o))  return NULL;
  // o is borrowed
  if (PyInstance_Check(o))
  {
    c = PyObject_GetAttrString(o, "__class__");
    if (!c) return NULL;
    i = PyObject_Compare(c, RTRType);
    Py_DECREF(c);
    if (PyErr_Occurred()) return NULL;
    if (i == 0) return smds_Correlator_data(self, o);
  }
  if (!PyObject_TypeCheck(o, &smds_CorrelatorType))
  {
    PyErr_SetString(PyExc_TypeError, 
      "Argument must be of type smds.Task.RTR or smds.Analysis.Correlator!");
    return NULL;
  }
  b = (smds_Correlator*)o;

  if (a->bw != b->bw)
  {
    PyErr_SetString(PyExc_ValueError,
      "Binwidth mis-match in smds_Correlator_append().");
    return NULL;
  }

  // Reorder a->shift so that pos[k] == 0 for all k
  for (k=1; k < (K+1); k++)
  {
    if (a->pos[k] == (k==1?-8:0)) continue;
    rr = a->pos[k];
    i = 0;
    do
    {
      swap[i++] = a->shift[POS(k,rr)];
      if (++rr == 8) rr = (k==1?-8:0);
    } while (rr != a->pos[k]);
    if (k==1)
    {
      for (i=0; i < 16; i++)
        a->shift[i] = swap[i];
      a->pos[k] = -8;
    }
    else
    {
      for (i=0; i < 8; i++)
        a->shift[POS(k,i)] = swap[i];
      a->pos[k] = 0;
    }
  }

  // Compute overlap correlation contributions for each lag, Steps 1-3
  {
    phi = (b->fill[1] < 8 ? 8 : 16-(b->fill[1]));
    // first 8 of k==1 (requires no summation afterward)
    for (rr = 0; rr < 8; rr++)
      for (i=0; i <= (rr < phi ? rr : phi-1); i++)
      {
        a->M_del[rr] += a->shift[rr-i];
        a->G[rr] += b->start[i] * a->shift[rr-i];
      }
  }
  for (k=1; k < (K+1); k++)
  {
    if (k==1) phi = (b->fill[1] < 8 ? 8-(b->fill[1]) : 0);
    else      phi = 8-(b->fill[k]);

    for (rr = 0; rr < 8; rr++)				// Lag
      for (i=0; i <= (rr < phi ? rr+8 : phi+rem-1); i++)  // Overlap position
      {
        a->M_del[POS(k,rr)] += a->shift[POS(k,rr)-i];
        a->G[POS(k,rr)] += b->start[POS(k-1,i)] * a->shift[POS(k,rr)-i];
      }
    if (a->fill[k] && !b->fill[k-1])
      for (rr=used; a->fill[k] && rr < phi+8; a->fill[k]--, rr++)
      {
        a->start[POS(k,8-(a->fill[k]))] = b->start[POS(k-1,rr)];
        used++;
      }
    
    // Sum shift/start for the next k
    rem = (rem+phi)/2;
    if (b->fill[k])
    {
      // Save a->shift for later use
      for (i=0; i < 8; i++)
        t2[i] = a->shift[POS(k,i)];
    }
    if (k != K)
      for (rr = POS(k,7), i=0; i < 8; i++,rr--)
      {
        a->shift[rr] = a->shift[rr-i] + a->shift[rr-i-1];
        b->start[rr] = b->start[rr-i] + b->start[rr-i-1];
      }
    if (b->fill[k-1])
    {
      // Restore a->shift from last iteration
      for (i=0; i < 8; i++)
        a->shift[POS(k-1,i)] = t1[i];
    }
    if (b->fill[k])
      for (i=0; i < 8; i++)
        t1[i] = t2[i];

    used >>= 1;
  }
  
  // Step 4
  for (i=0; i < 8*(K+1); i++)
  {
    a->G[i] += b->G[i];
    a->M_del[i] += b->M_del[i];
  }
  for (i=0; i < (K+1); i++)
  {
    a->M_dir[i] += b->M_dir[i];
  }
  a->G0 += b->G0;
  a->sum += b->sum;
  a->dur += b->dur;
  a->parity = b->parity;

  // Step 5, copy B's shift register to A
  for (k=1; k < (K+1) && !(b->fill[k]); k++)
  {
    if (k == 1)
      for (i=0, rr=0; i < 8; i++, rr += 2)
        t1[i] = a->shift[rr] + a->shift[rr+1];
    else
    {
      for (i=0, rr=0; i < 4; i++, rr += 2)
        t1[i] = t1[rr] + t1[rr+1];
      for (i=4, rr=0; i < 8; i++, rr += 2)
        t1[i] = a->shift[POS(k,rr)] + a->shift[POS(k,rr+1)];
    }
    a->pos[k] = b->pos[k];
    a->acc[k] = b->acc[k];
    for (rr=(k==1?-8:0); rr < 8; rr++)
      a->shift[POS(k,rr)] = b->shift[POS(k,rr)];
  }

  if (k == 1)		// k is the first group of b->start that isn't full
  {				// Step 6
    phi = 16-(b->fill[k]);
    rr = 16-phi;
    if (phi % 2)
    {
      a->acc[k+1] = a->shift[rr];
      a->parity |= (1 << k);
      rr++;
    } else {
      a->acc[k+1] = 0;
      a->parity &= !(1 << k);
    }

    for (i=0; rr < 16; rr += 2, i++)
      t1[i] = a->shift[rr] + a->shift[rr+1];

    a->pos[k] = 16-phi;
    for (i=(a->pos[k]), rr=(b->pos[k]); i < 16; i++)
    {
      a->shift[i] = b->shift[POS(k,rr)];
      if (++rr == 8) rr = -8;
    }
    rem = phi/2;

  } else if (k < (K+1)) {		// Step 7
    phi = 8-(b->fill[k]);
    for (rr=0, i=0; rr < 8; rr += 2, i++)
      t2[i] = a->shift[POS(k,rr)] + a->shift[POS(k,rr+1)];
    
    for (i=0, rr=(b->pos[k]); i < phi; i++)
    {
      a->shift[POS(k,i)] = b->shift[POS(k,rr)];
      if (++rr == 8) rr = 0;
    }
    for (rr=0; i < 8; i++, rr++)
      a->shift[POS(k,i)] = t1[rr];
    a->acc[k] = b->acc[k];
    
    if (phi % 2)
    {
      if (k != K) a->acc[k+1] = t1[rr];
      a->parity |= (1 << k);
      rr++;
    } else {
      if (k != K) a->acc[k+1] = 0;
      a->parity &= !(1 << k);
    }

    for(i=0; rr < 8; rr += 2, i++)
      t1[i] = t1[rr] + t1[rr+1];
    for(rr=0; rr < 4; rr++, i++)
      t1[i] = t2[rr];
    
    rem = phi/2+4;
  }  

  // Step 8
  for (k++; k < (K+1); k++)
  {
    if (!rem)
    {
      if (k != K) a->acc[k+1] = 0;
      a->parity &= !(1 << k);
      continue;
    }
    
    a->pos[k] = 8-rem;
    for (i=0, rr=(a->pos[k]); rr < 8; i++, rr++)
    {
      t2[i] = a->shift[POS(k,rr)];
      a->shift[POS(k,rr)] = t1[i];
    }
    
    if (rem % 2)
    {
      if (k != K) a->acc[k+1] = t2[0];
      a->parity |= (1 << k);
      rr = 1;
    } else {
      if (k != K) a->acc[k+1] = 0;
      a->parity &= !(1 << k);
      rr = 0;
    }
    for (i=0; rr < rem; rr+=2, i++)
      t1[i] = t2[rr] + t2[rr+1];

    rem >>= 1;
  }

  // Step 10, fill the zero-delay register
  for (k=2, i=1, rr=-8; k < (K+1); k++)
  {
    // k is the current zero we're filling
    // i/rr is the current position in the shift register (xsi)
    // bin width of level k is 2**(k-1)
    // if parity[k], fill up zero with half a bin width, i.e. 2**(k-2)
    // when rr == 8, add in acc[i+1] if parity[i]
    a->zero[k] = 0.0;
    if (a->parity & (1 << (k-1)))
    {
      for (rem = 1 << (k-2); rem > 0; rem -= (1 << (i-1)))
      {
        if (rr == 8)
        {
          if (a->parity & (1 << i))
          {
            a->zero[k] += a->acc[i+1];
            rem -= 1 << (i-1);
            if (rem < 0) break;	// Error
          }
          i++;
          rr = 0;
          if (rem == 0) break;	// Normal termination
        }
        a->zero[k] += a->shift[POS(i,rr++)];
      }
      if (rem < 0) 
      {
        PyErr_SetString(PyExc_Exception, 
          "smds_Correlator_append(): Error calculating zero[]!");
        return NULL;
//        g_error("smds_Correlator_append(): Error calculating zero[%d], i=%d, "
//                "rr=%d, rho=%d, par = %d", k, i, rr, rem, a->parity);
      }
    }
  }

  // Done.
  Py_RETURN_NONE;
}

PyObject * smds_Correlator_getstate(PyObject *self)
{
  smds_Correlator *c = (smds_Correlator*)self;
  PyObject *r = NULL, *o = NULL;
  PyArrayObject *a = NULL;
  double *d;
  int i, dim;

  r = PyTuple_New(4);
  if (!r) return NULL;

  dim = 4*8*(K+1);
  o = PyArray_FromDims(1, &dim, PyArray_DOUBLE);
  if (!o) goto bail;
  a = (PyArrayObject*)o;
  d = (double*)a->data;
  for (i=0; i < 8*(K+1); i++)
  {
    d[0*8*(K+1)+i] = c->shift[i];
    d[1*8*(K+1)+i] = c->start[i];
    d[2*8*(K+1)+i] = c->G[i];
    d[3*8*(K+1)+i] = c->M_del[i];
  }
  if (PyTuple_SetItem(r, 0, o)) goto bail; // true on error; steals reference

  dim = 3*(K+1);
  o = PyArray_FromDims(1, &dim, PyArray_DOUBLE);
  if (!o) goto bail;
  a = (PyArrayObject*)o;
  d = (double*)a->data;
  for (i=0; i < (K+1); i++)
  {
    d[0*(K+1)+i] = c->M_dir[i];
    d[1*(K+1)+i] = c->zero[i];
    d[2*(K+1)+i] = c->acc[i];
  }
  if (PyTuple_SetItem(r, 1, o)) goto bail;
  
  dim = 2*(K+1);
  o = PyArray_FromDims(1, &dim, PyArray_CHAR);
  if (!o) goto bail;
  a = (PyArrayObject*)o;
  for (i=0; i < (K+1); i++)
  {
    a->data[0*(K+1)+i] = c->pos[i];
    a->data[1*(K+1)+i] = c->fill[i];
  }
  if (PyTuple_SetItem(r, 2, o)) goto bail;
  
  o = Py_BuildValue("(dddii)", c->G0, c->sum, c->bw, (int)c->dur,
  					(int)c->parity);
  if (!o) goto bail;
  if (PyTuple_SetItem(r, 3, o)) goto bail;
  o = NULL;

  return r;
bail:
  if (r) { Py_DECREF(r); }
  if (o) { Py_DECREF(o); }
  return NULL;
}

PyObject * smds_Correlator_setstate(PyObject *self, PyObject *args)
{
  smds_Correlator *c = (smds_Correlator*)self;
  PyObject *state = NULL, *o = NULL;
  PyArrayObject *a = NULL;
  double *d;
  int i;

#define STR(x) #x
#define DIE(line) { PyErr_SetString(PyExc_TypeError, \
		"Invalid Correlator state object. (" STR(line) ")"); return NULL; }

/*
  double shift[8*(K+1)], start[8*(K+1)], G[8*(K+1)], M_del[8*(K+1)], 
         M_dir[(K+1)], zero[(K+1)], acc[(K+1)];
  double G0, sum, bw;
  char pos[(K+1)];		// Starting positions in the round-robin arrays
  char fill[(K+1)];	// Bins of start left to fill
  Int32 dur;
  Int16 parity;		// k'th bit true if (k+1)[0] is half-filled
*/

  if (!PyArg_ParseTuple(args, "O!", &PyTuple_Type, &state)) return NULL;
		// state is borrowed
  if (PyTuple_Size(state) != 4) DIE(__LINE__);
  
  o = PyTuple_GetItem(state, 0);	// o is borrowed
  if (!o) return NULL;
  a = (PyArrayObject*)PyArray_ContiguousFromObject(o, PyArray_DOUBLE, 1, 1);
  if (!a || a->dimensions[0] != 4*8*(K+1)) DIE(__LINE__);
  d = (double*)a->data;
  for (i=0; i < 8*(K+1); i++)
  {
    c->shift[i] = d[0*8*(K+1)+i];
    c->start[i] = d[1*8*(K+1)+i];
    c->G[i]     = d[2*8*(K+1)+i];
    c->M_del[i] = d[3*8*(K+1)+i];
  }
  Py_DECREF(a);
  
  o = PyTuple_GetItem(state, 1);
  if (!o) return NULL;
  a = (PyArrayObject*)PyArray_ContiguousFromObject(o, PyArray_DOUBLE, 1, 1);
  if (!a) return NULL;
  if (a->dimensions[0] != 3*(K+1)) { Py_DECREF(a); DIE(__LINE__); }
  d = (double*)a->data;
  for (i=0; i < (K+1); i++)
  {
    c->M_dir[i] = d[0*(K+1)+i];
    c->zero[i]  = d[1*(K+1)+i];
    c->acc[i]   = d[2*(K+1)+i];
  }
  Py_DECREF(a);
  
  o = PyTuple_GetItem(state, 2);
  if (!o) return NULL;
  a = (PyArrayObject*)PyArray_ContiguousFromObject(o, PyArray_CHAR, 1, 1);
  if (!a) return NULL;
  if (a->dimensions[0] != 2*(K+1)) { Py_DECREF(a); DIE(__LINE__); }
  for (i=0; i < (K+1); i++)
  {
    c->pos[i]  = a->data[0*(K+1)+i];
    c->fill[i] = a->data[1*(K+1)+i];
  }
  Py_DECREF(a);

  o = PyTuple_GetItem(state, 3);
  if (!o) return NULL;
  if (!PyTuple_Check(o) || PyTuple_Size(o) != 5) DIE(__LINE__);
  if (!PyArg_ParseTuple(o, "dddii", &c->G0, &c->sum, &c->bw, &c->dur, 
  				&c->parity)) DIE(__LINE__);

  Py_RETURN_NONE;
}

PyObject * smds_Correlator_reduce(PyObject *self)
{
  PyObject *o, *r;
  o = smds_Correlator_getstate(self);
  if (!o) return NULL;
  r = Py_BuildValue("(O(d)O)", &smds_CorrelatorType, 0.0, o);
  Py_DECREF(o);
  return r;
}
