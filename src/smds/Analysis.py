# smds.Analysis module
#  Structures and routines for data analysis
#  $Id: Analysis.py,v 1.15 2009/06/11 18:31:15 mculbert Exp $

# online analyses must provide four methods:
#    prepare(task)
#    analyze(numBins, data[array(Short)], posInDataToProcess)
#    finalize()
#    combine(analysis)
#       finalize() is called after all calls to analyze() and again after
#       all calls to combine()
# offline analyses must provide one method:
#    analyze(task)

from warnings import warn
from sys import stdout
from Scientific.Functions.LeastSquares import leastSquaresFit
from types import FunctionType
import xml.dom, re, copy, smds
import smds.Calc

analyses = {}		# Name -> class

class Analysis :
  """Generic class for data analysis results."""
  online = False
  def fromXML(n) :
    """fromXML(n) --> analysis

    Determines the appropriate analysis type and returns an analysis object 
    from the given node."""
    if (n == None or n.tagName != 'Analysis') : return None
    return analyses[n.getAttribute('Name')].fromXML(n)
  fromXML = staticmethod(fromXML)
  def copy(self) :
    r = self.__class__()
    r.__dict__.update(copy.deepcopy(self.__dict__))
    return r
  def strip(self) :
    return self

class SEDH (Analysis) :
  """Single-Event Duration Histogram"""
  online = True
  def __init__(self, binWidth = 0.1, threshold = 0) :
    """SEDH([binWidth], [threshold]) --> sedh

    Takes optional binWidth (ms) and threshold (phot/ms)."""
    self.binWidth = binWidth
    self.threshold = threshold	
    self.data = []
    self.tabulator = None
  def toXML(self, doc) :
    """sedh.toXML(doc)

    Returns this SEDH as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("SEDH")
    n.setAttribute("binWidth", str(self.binWidth))
    n.setAttribute("threshold", str(self.threshold))
    if len(self.data) :
      d = doc.createTextNode(
        '\n'.join([ "%g,%d" % (self.binWidth*i, self.data[i]) for
        			i in range(len(self.data))
        			if self.data[i] > 0 ])
      )
      n.appendChild(d)
    return n
  def fromXML(n) :
    """SEDH.fromXML(n) --> sedh

    Returns a new smds.Analysis.SEDH created from the xml.dom node n."""
    if (n == None or n.tagName != 'SEDH') : return None
    r = SEDH(float(n.getAttribute('binWidth')), 
                  float(n.getAttribute('threshold')) )
    d = smds.strFromTxtNodes(n.childNodes)
    d = re.sub("\r", "\n", d)
    d = re.sub("[ \t]", '', d)
    d = re.sub("\n{2,}", "\n", d)
    if d == '' : data = []
    else : data = d.split("\n")
    if len(data) and data[0] == '' : data = data[1:]
    if len(data) and data[-1] == '' : data = data[:-1]
    if len(data) :
      for (t, x) in [ p.split(",") for p in data ] :
        i = int(float(t)/r.binWidth+0.5)
        if i >= len(r.data) : r.data += [0]*(i-len(r.data)+1)
        r.data[i] = int(x)
    return r
  fromXML = staticmethod(fromXML)
  def prepare(self, task) :
    self.tabulator = Tabulator(task.p.binWidth, self.threshold,	self.binWidth)
  def finalize(self) :
    if self.tabulator :
      self.data = self.tabulator.fetch().data
      self.tabulator = None
  def analyze(self, numBins, data, pos) :
    rtr = smds.Task.RTR()
    rtr.data = data[pos:pos+numBins]
    self.tabulator.tabulate(rtr)
  def combine(self, a) :
    if (len(a.data) > len(self.data)) :
      data = a.data.copy()
      for i in range(len(self.data)) :
        data[i] += self.data[i]
      self.data = data
    else :
      for i in range(len(a.data)) :
        self.data[i] += a.data[i]
  def strip(self) :
    r = SEDH()
    r.__dict__ = self.__dict__.copy()
    r.tabulator = None
    r.data = []
    return r
  def diff(a, b) :
    """sedh_a.diff(sedh_b) --> diff
Returns a figure of comparison between SEDH instances a and b.
    
    Assumes that a is a "true" curve.  The figure of comparison is defined as:
      Sum(t){ [(b(t)-a(t))/a(t)]^2 if a(t) > 0 and b(t) > 0,
              1 if a(t) > 0 xor b(t) > 0,
              0 otherwise }"""
    if (len(a.data) < len(b.data)) :
      swap = a
      a = b
      b = swap
    sum = 0
    for i in range(len(b.data)) :
      if (a.data[i] > 0) :
        if (b.data[i] > 0) :
          sum += ((float(b.data[i])-float(a.data[i]))/float(a.data[i]))**2
        else :
          sum += 1.0
      else :
        if (b.data[i] > 0) :
          sum += 1.0
    for i in range(len(b.data), len(a.data)) :
      if (a.data[i] > 0) :
        sum += 1.0
    return sum
  def ssr(a, b) :
    """sedh_a.ssr(sedh_b) --> ssr

Returns the sum of squared residuals between smds.SEDH instances a and b.
    
    Assumes that a is a "true" curve."""
    if (len(a.data) < len(b.data)) :
      swap = a
      a = b
      b = swap
    sum = 0
    for i in range(len(b.data)) :
      sum += (a.data[i] - b.data[i])**2
    for i in range(len(b.data), len(a.data)) :
      sum += a.data[i]**2
    return sum

class PCH (Analysis) :
  """Photon Counting Histogram"""
  online = True
  def __init__(self, binWidth = None) :
    """PCH([binWidth (ms)]) --> pch"""
    self.binWidth = binWidth
    self.data = []
    self.tabulator = None
  def toXML(self, doc) :
    """pch.toXML(doc)

    Returns this PCH as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("Analysis")
    n.setAttribute("Name", "PCH")
    n.setAttribute("binWidth", str(self.binWidth))
    if len(self.data) :
      d = doc.createTextNode(
        '\n'.join([ "%g,%d" % (i, self.data[i]) for
        			i in range(len(self.data))
        			if self.data[i] > 0 ])
      )
      n.appendChild(d)
    return n
  def fromXML(n) :
    """PCH.fromXML(n) --> pch

    Returns a new smds.Analysis.PCH created from the xml.dom node n."""
    r = PCH(float(n.getAttribute('binWidth')))
    d = smds.strFromTxtNodes(n.childNodes)
    d = re.sub("\r", "\n", d)
    d = re.sub("[ \t]", '', d)
    d = re.sub("\n{2,}", "\n", d)
    if d == '' : data = []
    else : data = d.split("\n")
    if len(data) and data[0] == '' : data = data[1:]
    if len(data) and data[-1] == '' : data = data[:-1]
    if len(data) :
      for (t, x) in [ p.split(",") for p in data ] :
        i = int(t)
        if i >= len(r.data) : r.data += [0]*(i-len(r.data)+1)
        r.data[i] = int(x)
    return r
  fromXML = staticmethod(fromXML)
  def prepare(self, task) :
    if not self.binWidth : self.binWidth = task.p.binWidth
    self.tabulator = Tabulator_PCH(task.p.binWidth, self.binWidth)
  def finalize(self) :
    if self.tabulator : 
      self.data = self.tabulator.fetch().data
      self.tabulator = None
  def analyze(self, numBins, data, pos) :
    rtr = smds.Task.RTR()
    rtr.data = data[pos:pos+numBins]
    self.tabulator.tabulate(rtr)
  def combine(self, a) :
    if (len(a.data) > len(self.data)) :
      data = a.data.copy()
      data[:len(self.data)] += self.data
      self.data = data
    else :
      self.data[:len(a.data)] += a.data
  def strip(self) :
    r = PCH()
    r.__dict__ = self.__dict__.copy()
    r.tabulator = None
    r.data = []
    return r

analyses['PCH'] = PCH

class IPCH (Analysis) :
  """Integrated Pulse Counting Histogram"""
  online = True
  def __init__(self, threshold = 0, binWidth = None) :
    """IPCH([threshold (phot/ms)], [binWidth (ms)]) --> ipch"""
    self.threshold = threshold
    self.binWidth = binWidth
    self.data = []
    self.tabulator = None
  def toXML(self, doc) :
    """ipch.toXML(doc)

    Returns this IPCH as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("Analysis")
    n.setAttribute("Name", "IPCH")
    n.setAttribute("threshold", str(self.threshold))
    n.setAttribute("binWidth", str(self.binWidth))
    if len(self.data) :
      d = doc.createTextNode(
        '\n'.join([ "%g,%d" % (i, self.data[i]) for
        			i in range(len(self.data))
        			if self.data[i] > 0 ])
      )
      n.appendChild(d)
    return n
  def fromXML(n) :
    """IPCH.fromXML(n) --> ipch

    Returns a new smds.Analysis.PCH created from the xml.dom node n."""
    r = IPCH(float(n.getAttribute('threshold')),
            float(n.getAttribute('binWidth')))
    d = smds.strFromTxtNodes(n.childNodes)
    d = re.sub("\r", "\n", d)
    d = re.sub("[ \t]", '', d)
    d = re.sub("\n{2,}", "\n", d)
    if d == '' : data = []
    else : data = d.split("\n")
    if len(data) and data[0] == '' : data = data[1:]
    if len(data) and data[-1] == '' : data = data[:-1]
    if len(data) :
      for (t, x) in [ p.split(",") for p in data ] :
        i = int(t)
        if i >= len(r.data) : r.data += [0]*(i-len(r.data)+1)
        r.data[i] = int(x)
    return r
  fromXML = staticmethod(fromXML)
  def prepare(self, task) :
    if not self.binWidth : self.binWidth = task.p.binWidth
    self.tabulator = Tabulator_IPCH(task.p.binWidth, self.threshold,
                                    self.binWidth)
  def finalize(self) :
    if self.tabulator : 
      self.data = self.tabulator.fetch().data
      self.tabulator = None
  def analyze(self, numBins, data, pos) :
    rtr = smds.Task.RTR()
    rtr.data = data[pos:pos+numBins]
    self.tabulator.tabulate(rtr)
  def combine(self, a) :
    if (len(a.data) > len(self.data)) :
      data = a.data.copy()
      data[:len(self.data)] += self.data
      self.data = data
    else :
      self.data[:len(a.data)] += a.data
  def strip(self) :
    r = IPCH()
    r.__dict__ = self.__dict__.copy()
    r.tabulator = None
    r.data = []
    return r

analyses['IPCH'] = IPCH

class ACF (Analysis) :
  """Autocorrelation Function"""
  online = True
  def __init__(self) :
    self.t = []
    self.G = []
  def toXML(self, doc, n = None) :
    """acf.toXML(doc) --> node

    Returns this ACF as a node for insertion into xml.dom.Document doc."""
    if (n == None) : n = doc.createElement("ACF")
    if len(self.t) :
      d = doc.createTextNode(
        '\n'.join([ "%g,%g" % (self.t[i], self.G[i])
        			for i in range(len(self.t)) ]) )
      n.appendChild(d)
    return n
  def fromXML(n, r = None) :
    """ACF.fromXML(n) --> acf

    Returns a new smds.Analysis.ACF created from the xml.dom node n."""
    if (n == None or (n.tagName != 'ACF' and r == None)) : return None
    from numpy import asarray
    if (r == None) : r = ACF()
    d = smds.strFromTxtNodes(n.childNodes)
    d = re.sub("\r", "\n", d)
    d = re.sub("[ \t]", '', d)
    d = re.sub("\n{2,}", "\n", d)
    if d == '' : data = []
    else :       data = d.split("\n")
    if len(data) and data[0] == '' : data = data[1:]
    if len(data) and data[-1] == '' : data = data[:-1]
    if len(data) :
      for (t, G) in [ p.split(",") for p in data ] :
        r.t.append(float(t))
        r.G.append(float(G))
    r.t = asarray(r.t)
    r.G = asarray(r.G)
    return r
  fromXML = staticmethod(fromXML)
  def prepare(self, task) :
    self.correlator = Correlator(task.p.binWidth)
  def finalize(self) :
    r = self.correlator.getACF()
    self.t = r.t
    self.G = r.G
  def analyze(self, numBins, data, pos) :
    rtr = smds.Task.RTR()
    rtr.data = data[pos:pos+numBins]
    self.correlator.append(rtr)
  def combine(self, a) :
    self.correlator.append(a.correlator)
  def strip(self) :
    return ACF()
  def write(self, f) :
    """acf.write(file)

    Writes this ACF to the given open file, one t,G pair per line."""
    for i in range(len(self.t)) :
      f.write("%G,%G\n" % (self.t[i], self.G[i]))
  def read(f) :
    """ACF.read(file) --> acf

    Creates a new ACF instance based on the data read in from the given open
    file.  The file should have one t,G pair per line."""
    from numpy import asarray
    r = ACF()
    for (t,G) in [ line.split(',') for line in f.readlines() ] :
      r.t.append(float(t))
      r.G.append(float(G))
    r.t = asarray(r.t)
    r.G = asarray(r.G)
    return r
  read = staticmethod(read)
  def ssr(self, b) :
    """acf_a.ssr(acf_b) --> ssr

    Calculates the difference between two autocorrelation curve.

    Difference is reported as the sum of squared residuals using a as the
    model curve.  If the abscissas from b do not match those from a, linear
    interpolations are constructed between each point of a for comparison."""
    def getValue(t, G, x) :
      low = -1
      high = len(t)
      while (high > low+1) :
        i = int((high+low)/2)
        if (t[i] == x) : return G[i]
        if (x < t[i])  : high = i
        else           : low = i
      i = low
      if (t[i] == x) : return G[i]
      if i == len(t)-1 : i = i-1
      return (lambda (x1, y1), (x2, y2), x : (y2-y1)/(x2-x1)*(x-x1)+y1) \
      		((t[i], G[i]), (t[i+1], G[i+1]), x)
    sum = 0
    for i in range(len(b.t)) :
      sum += (getValue(self.t, self.G, b.t[i])-b.G[i])**2
    return sum

class ACF_func (ACF) :
  """Autocorrelation of a Function of the Data """
  def __init__(self, func = lambda (rtr,working) : rtr) :
    ACF.__init__(self)
    self.func = func
    self.name = func.__name__
    self.working = {}
    self.freeze()
  def freeze(self) :
    """Private."""
    from base64 import b64encode
    if (type(self.func) == FunctionType) :
      self.func = b64encode(smds.Calc.mDump(self.func.func_code, 1))
  def thaw(self) :
    """Private."""
    from base64 import b64decode
    if (type(self.func) != FunctionType) :
      self.func = FunctionType(smds.Calc.mLoad(b64decode(self.func)),
                               smds.Calc.env)
      self.func.__name__ = self.name
  def analyze(self, numBins, data, pos) :
    rtr = smds.Task.RTR()
    rtr.data = self.func(data[pos:pos+numBins].copy(),self.working).astype('s')
    self.correlator.append(rtr)
  def prepare(self, task) :
    self.working['task'] = task
    self.thaw()
    ACF.prepare(self, task)
  def finalize(self) :
    self.working = {}		# Free working memory
    self.freeze()
    ACF.finalize(self)
  def strip(self) :
    ret = ACF_func()
    self.freeze()
    ret.func = self.func
    ret.name = self.name
    return ret
  def toXML(self, doc) :
    n = doc.createElement("Analysis")
    n.setAttribute("Name", "ACF_func")
    n.setAttribute("func_name", self.name)
    if (type(self.func) == FunctionType) :
      from base64 import b64encode
      n.setAttribute("func_data", b64encode(smds.Calc.mDump(
                                    self.func.func_code, 1)) )
    else :
      n.setAttribute("func_data", self.func)
    return ACF.toXML(self, doc, n)
  def fromXML(n) :
    r = ACF_func()
    r = ACF.fromXML(n, r)
    r.name = n.getAttribute("func_name")
    r.func = n.getAttribute("func_data")
    return r
  fromXML = staticmethod(fromXML)

analyses['ACF_func'] = ACF_func

#####################
#### Fitting Routines
#####################

# This must be subclassed for the actual fits.
#  Supply initialGuess, chooseXY, func, and varCoeffs[]
#  Add the new fit to FitTypes[]
#  There could be a spiffy-cool way to map the coeffNames to the coeff array
#   dynamically, but I don't feel like putting the effort into the extra code.
class Fit (Analysis) :
  """A generic class for data fitting results."""
  name = "No Fit"
  def __init__(self) :
    self.coeff = {}
    self.chiSq = 0
  def analyze(self, task) :
    self.initialGuess(task)
    p = [ self.coeff[c] for c in self.varCoeffs ]
    (p, self.chisq) = leastSquaresFit(self.func, p, self.chooseXY(task))
    for i in range(len(self.varCoeffs)) : 
      self.coeff[self.varCoeffs[i]] = p[i]
    del self.func
  def toXML(self, doc) :
    """fit.toXML(doc) --> node

    Returns this Fit as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("Fit")
    n.setAttribute("type", self.name)
    if self.chiSq != 0.0 : n.setAttribute("chiSq", self.chisq)
    for i in self.coeff.keys() :
      c = doc.createElement("Coeff")
      c.setAttribute("name", i)
      c.setAttribute("value", "%g" % (self.coeff[i],))
      n.appendChild(c)
    return n
  def fromXML(n) :
    """Fit.fromXML(n) --> fit

    Returns a new smds.Analysis.Fit created from the xml.dom node n."""
    if (n == None or n.tagName != 'Fit') : return None
    try :
      r = FitTypes[n.getAttribute("type")]()
    except KeyError, k :
      r = Fit()
      warn("Unknown Fit type: " + k[0], stacklevel = 2)
    if n.hasAttribute('chiSq') : r.chiSq = float(n.getAttribute('chiSq'))
    for coeff in n.childNodes :
      if coeff.nodeType == xml.dom.Node.ELEMENT_NODE and 	\
         coeff.tagName == 'Coeff' :
        (name, value) = coeff.getAttribute('name'), 		\
			float(coeff.getAttribute('value'))
#        		float(smds.strFromTxtNodes(coeff.childNodes))
        r.coeff[name] = value
    return r
  fromXML = staticmethod(fromXML)

# A generic ACF fitting class
#   supply initialGuess, func, and coeffNames[]
class FitACF (Fit) :
  def chooseXY(self, task) :
    for a in task.anal :
      if a.__class__ == ACF :
        return [ (a.t[i], a.G[i]) for i in range(len(a.t)) ]
    raise AttributeError, "ACF not found for ACF fit."

class FCS_2D (FitACF) :
  name = "FCS 2D"
  varCoeffs = ['N', 'D']
  def initialGuess(self, task) :
    self.coeff = { 'N' : 0.1, 'D' : 1.0e-7, 'R' : task.p.radius }
    self.func = lambda p, x :						\
	(1./p[0])*(1.+4.*p[1]*x/(self.coeff['R']*1e-7)**2)**(-1)+1.0

class FCS_2D_bkg (FCS_2D) :
  name = "FCS 2D bkg"
  def initialGuess(self, task) :
    FCS_2D.initialGuess(self, task)
    self.coeff['I'] = task.results.avg_I
    self.coeff['B'] = task.p.bkg
    self.func = lambda p, x : ((self.coeff['I']-self.coeff['B']) /
    			self.coeff['I'])**2 * (self.func(p,x)-1.0) + 1.0

class FCS_2D_bkg_gamma (FCS_2D_bkg) :
  name = "FCS 2D bkg gamma"
  def initialGuess(self, task) :
    FCS_2D_bkg.initialGuess(self, task)
    self.func = lambda p, x : 0.5 * (self.func(p,x)-1.0) + 1.0

class FCS_3D (FitACF) :
  name = "FCS 3D"
  varCoeffs = ['N', 'D']
  def initialGuess(self, task) :
    self.coeff = { 'N' : 0.1, 'D' : 1.0e-7, 
    			'R' : task.p.radius, 'Z' : task.p.Z }
    self.func = lambda p, x :						\
	(1./p[0]) * (1.+4.*p[1]*x/(self.coeff['R']*1e-7)**2)**(-1) *	\
        (1.+4.*p[1]*x/(self.coeff['Z']*1e-7)**2)**(-0.5) + 1.0

class FCS_3D_bkg (FCS_3D) :
  name = "FCS 3D bkg"
  def initialGuess(self, task) :
    FCS_3D.initialGuess(self, task)
    self.coeff['I'] = task.results.avg_I
    self.coeff['B'] = task.p.bkg
    self.func = lambda p, x : ((self.coeff['I']-self.coeff['B']) /
                     self.coeff['I'])**2 * (self.func(p,x)-1.0) + 1.0

class FCS_3D_bkg_gamma (FCS_3D_bkg) :
  name = "FCS 3D bkg gamma"
  def initialGuess(self, task) :
    FCS_3D_bkg.initialGuess(self, task)
    self.func = lambda p, x : 0.35 * (self.func(p,x)-1.0) + 1.0

FitTypes = { None : Fit, '' : Fit, 'No Fit' : Fit,
             'FCS 2D' : FCS_2D,
             'FCS 2D bkg' : FCS_2D_bkg,
             'FCS 2D bkg gamma' : FCS_2D_bkg_gamma,
             'FCS 3D' : FCS_3D,
             'FCS 3D bkg' : FCS_3D_bkg,
             'FCS 3D bkg gamma' : FCS_3D_bkg_gamma,
           }

## Utility fits.

from numpy.oldnumeric.linear_algebra import solve_linear_equations
from numpy import dot, transpose

def quadFit(X, Y) :
  """quadFit(X, Y) --> [a1, a2, ... , b1, b2, ... , c]

  Fits (X,Y) data to a generalized quadratic:
    X, Y -- arrays of the (x,y) data, where x is either a scalar
					or a tuple (x1, x2, x3, ...)
    Returns [a1, a2, ..., b1, b2, ... , c]
    such that y = a1*x1^2 + a2*x2^2 + .... b1*x1 + b2*x2 + ... + c"""
  from numpy import asarray, reshape, ones, concatenate
  X = asarray(X)
  if len(X.shape) == 1 : X = reshape(X, (X.shape[0], 1))
  A = concatenate((X**2, X, ones((X.shape[0], 1))), 1)
  alpha = dot(transpose(A), A)
  beta  = dot(transpose(A), Y)
  return solve_linear_equations(alpha, beta)

def quadMin(A) :
  """quadMin(A) -> (x,y)

  Returns the (x,y) minimum of the generalized quadratic defined by
    A = [a1, a2, ... , b1, b2, ... , c]
    such that y = a1*x1^2 + a2*x2^2 + ... + b1*x1 + b2*x2 + ... + c"""
  from numpy import asarray, concatenate
  A = asarray(A)
  num = (len(A)-1)/2
  (a, b, c) = (A[:num], A[num:-1], A[-1])
  minX = -b/(2.*a)
  return concatenate((minX, (sum(a*minX**2 + b*minX) + c,)))

#######################
##### Matching Routines
#######################

from smds.opt import fmin

def compareSEDH(sedh, t) :
  """compareSEDH(sedh, t) --> ssr

  Calculate the SSR between an SEDH with the first SEDH in a given task."""
  return sedh.ssr([ a for a in t.anal if isinstance(a, SEDH) ][0])

def compareACF(ac, t) :
  """compareACF(ac, t) --> ssr

  Calculates the SSR between an ACF and the ACF of a given task."""
  return ac.ssr([ a for a in t.anal if isinstance(a, ACF) ][0])

def match(data, guess, mapper, compare, xtol=0.01, ftol=0.01) :
  """match(data, guess, mapper, compare, [xtol], [ftol])
Search for the parameters that match the given data.

  data    : the data to match
  guess   : a list comprising the values of the initial guess
  mapper  : a function taking a single list of values for the parameters
            that will vary during the match and map them into (i.e. return)
            an smds.Task or None if out of bounds
  compare : a function taking data and the task returned during a given
            function evaluation and return a figure of merit
  Returns : (x_min, task_matched, merit)
    x_min        : list of matched parameters
    task_matched : the Task of the match
    merit        : the minimum figure of merit for the match
  """
  store = {}
  def func(v) :
    try:
      t = store[tuple(v)]
    except KeyError:
      t = mapper(v)
      if not t : return 1e100
      batch = smds.Batch([t])
      batch.addToQueue()
      smds.clearQueue()
      t = batch.data[0]
      store[tuple(v)] = t
    ret = compare(data, t)
    print "func eval ["+', '.join(["%g" % (x,) for x in v])+"] -> "+str(ret)
    stdout.flush()
    return ret
  x_min = fmin(func, list(guess), disp=0, xtol=xtol, ftol=ftol)
  t_min = store[tuple(x_min)]
  return (x_min, t_min, compare(data, t_min))

