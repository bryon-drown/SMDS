# smds.Task module
#  Structures for simulation tasks
#  $Id: Task.py,v 1.9 2008/09/15 17:19:42 mculbert Exp $

import socket, time, re, copy, numpy, xml.dom, smds
from smds.loadintens import _loadintens, _writeintens, _loadpot, _writepot

Status_Ready = 1
Status_Waiting = 2
Status_Running = 3
Status_Processing = 4
Status_Done = 5
Status_Error = 0

intensCache = {}

class Intens :
  """Intensity Profiles."""
  def __init__(self, fn = '') :
    self.name = fn
    try :
      self.load()
      self.okay = True
    except Exception as e:
      self.okay = False
      if fn != '' :
        smds.msg("Error loading intensity profile <%s>!\n%s: %s" %
                    (fn, type(e).__name__, str(e)), smds.MSG_WARN)
  def load(self) :
    """Load the intensity profile.
   File must be of the format:
        binWidth_x binWidth_z numBins_x numBins_z data
   where number of datapoints must be exactly
        (2*numBins_z+1) * (2*numBins_x+1)^2
   thus, data[x,y,z] goes from [-numBins, numBins]."""
    if (self.name is '') :
      raise ValueError, "No filename supplied for Intensity profile load."
    if self.name in intensCache :
      if not isinstance(intensCache[self.name], Intens) :
        raise ValueError, name + " exists in memory, but is not an smds.Task.Intens"
      self.__dict__ = intensCache[self.name].__dict__.copy()
      return
    (self.binWidth_x, self.binWidth_z, self.numBins_x, self.numBins_z,
	self.data) = _loadintens(self.name)
    intensCache[self.name] = self
  def write(self, filename) :
    """intens.write(filename)
    Write the intensity profile to a file in standard SMDS format."""
    if not self.okay :
      raise ValueError, "Intensiy profile is not valid."
    _writeintens(filename, self.binWidth_x, self.binWidth_z, self.numBins_x,
                 self.numBins_z, self.data)
  def __deepcopy__(self, memo) :
    return copy.copy(self)
  def strip(self) :
    r = Intens()
    r.name = self.name
    return r
  def getData(self) :
    """intens.getData() -> array[x,y,z]
    Returns the intensity profile as a numpy array with three dimensions."""
    if not self.okay :
      raise ValueError, "Intensity profile is not valid."
    numX = numY = self.numBins_x*2+1
    numZ = self.numBins_z*2+1
    return numpy.transpose(numpy.reshape(self.data, 
                                    (numZ, numX, numY)), (1, 2, 0))
  def fromData(name, bwx, bwz, data) :
    """smds.Task.Intens.fromData(name, binWidth_xy, binWidth_z, array[x,y,z])
    Create an smds.Task.Intens from a three-dimensional array.
    The length in X and Y must be equal, and all dimensions must be odd.
    Bin widths are in nm and must be integers.  The name must be unique.
    Profile will be normalized."""
    if name in intensCache :
      raise ValueError, "Intens <%s> already exists." % name
    data = numpy.array(data, 'd', 0)
    (Nx, Ny, Nz) = data.shape
    if not (Nx == Ny and Nx % 2 and Nz % 2) :
      raise ValueError, "Data is not the correct shape."
    r = Intens()
    r.data = numpy.reshape(numpy.transpose(data, (2, 0, 1)), (Nx*Ny*Nz,))
    r.data /= numpy.maximum.reduce(r.data)
    r.binWidth_x = bwx
    r.binWidth_z = bwz
    r.numBins_x = (Nx-1)/2
    r.numBins_z = (Nz-1)/2
    r.name = name
    r.okay = True
    smds.Task.intensCache[r.name] = r
    return r    
  fromData = staticmethod(fromData)
  def fromFunc(name, Nx, Nz, bwx, bwz, func) :
    """smds.Task.Intens.fromFunc( name, Nx, Nz, bwx, bwz, func(X,Y,Z) )
    Create an smds.Task.Intens with the given size from a function.
    Binwidths are in nm and must be integers.  The name must be unique.
    Dimensions will be incremented if they are even.
    The function takes three parameters (X, Y, Z) in nm and returns the
    value of the intensity profile.  Profile will be normalized."""
    if name in intensCache :
      raise ValueError, "Intens <%s> already exists." % name
    if not(Nx % 2) : Nx += 1
    if not(Nz % 2) : Nz += 1
    Ny = Nx
    from numpy import transpose, repeat
    X = (transpose(repeat([repeat([range(-(Nx/2), Nx/2+1)], Ny, 0)], Nz, 0),
                        (2, 1, 0))*float(bwx))
    Y = (transpose(repeat([repeat([range(-(Ny/2), Ny/2+1)], Nx, 0)], Nz, 0),
                        (1, 2, 0))*float(bwx))
    Z = (repeat([repeat([range(-(Nz/2), Nz/2+1)], Ny, 0)], Nx, 0)*float(bwz))
    return Intens.fromData(name, bwx, bwz, func(X,Y,Z))
  fromFunc = staticmethod(fromFunc)

class PotProfile :
  """aHL Potential Profiles."""
  def __init__(self, fn = '') :
    self.name = fn
    try :
      self.load()
      self.okay = True
    except Exception as e:
      self.okay = False
      if fn != '' :
        smds.msg("Error loading potential profile <%s>!\n%s: %s" %
                    (fn, type(e).__name__, str(e)), smds.MSG_WARN)
  def load(self) :
    """Load the potential profile.
     aHL Potential Profiles must be of the format:
          binWidth_A_z binWidth_A_r numBins_A_z numBins_A_r data_A
          binWidth_B_z binWidth_B_r numBins_B_z numBins_B_z data_B
where A refers to everywhere outside the channel (from the membrane, z=0,
up) and B refers to in the channel (from the top of the cap, z=capHeight,
down).  Note that the two regions overlap, and that the direction of data_A
is in *increasing* z and the direction of data_B is in *decreasing* z.  Data
extend radially from r=0 out.  The two data segments have two tables of
numBins_?_z rows and numBins_?_r columns.  The first table indicates the
electrophoretic flow in z and the second in r."""
    if (self.name is '') :
      raise ValueError, "No filename supplied for Potential profile load."
    if self.name in intensCache :
      if not isinstance(intensCache[self.name], PotProfile) :
        raise ValueError, name + " exists in memory, but is not an smds.Task.PotProfile"
      self.__dict__ = intensCache[self.name].__dict__.copy()
      return
    (self.bwAz, self.bwAr, self.NAz, self.NAr, self.dataA,
     self.bwBz, self.bwBr, self.NBz, self.NBr, self.dataB) = _loadpot(self.name)
    intensCache[self.name] = self
  def write(self, filename) :
    """potProfile.write(filename)
    Write the potential profile to a file in standard SMDS format."""
    if not self.okay :
      raise ValueError, "Potential profile is not valid."
    _writepot(filename, self.bwAz, self.bwAr, self.NAz, self.NAr, self.dataA,
                        self.bwBz, self.bwBr, self.NBz, self.NBr, self.dataB)
  def __deepcopy__(self, memo) :
    return copy.copy(self)
  def strip(self) :
    r = PotProfile()
    r.name = self.name
    return r
  def getData(self) :
    """intens.getData() -> (arrayA[z,r,flow dimension], arrayB[z,r,flow dimension])
    Returns the potential profile as two numpy arrays with three dimensions.
    Last dimension indicates flow direction, 0 == z, 1 == r."""
    if not self.okay :
      raise ValueError, "Potential profile is not valid."
    return (numpy.transpose(numpy.reshape(self.dataA, 
                                    (2, self.NAz, self.NAr)), (1, 2, 0)),
            numpy.transpose(numpy.reshape(self.dataB,
                                    (2, self.NBz, self.NBr)), (1, 2, 0)) )
  def fromData(name, bwAz, bwAr, dataA, bwBz, bwBr, dataB) :
    """smds.Task.PotProfile.fromData(name, bwAz, bwAr, arrayA[z,r,flow dim], bwBz, bwBr, arrayB[z,r,flow dim])
    Create an smds.Task.PotProfile from two three-dimensional arrays.
    Bin widths are in nm.  The name must be unique."""
    if name in intensCache :
      raise ValueError, "Intens <%s> already exists." % name
    dataA = numpy.array(dataA, 'd', copy=False)
    dataB = numpy.array(dataB, 'd', copy=False)
    (NAz, NAr, NAf) = dataA.shape
    (NBz, NBr, NBf) = dataB.shape
    if (NAf != 2) :
      raise ValueError, "Data for Zone A is not the correct shape."
    if (NBf != 2) :
      raise ValueError, "Data for Zone B is not the correct shape."
    r = PotProfile()
    r.dataA = numpy.reshape(numpy.transpose(dataA, (2, 0, 1)), (NAz*NAr*2,))
    r.bwAz = bwAz
    r.bwAr = bwAr
    r.NAz = NAz
    r.NAr = NAr
    r.dataB = numpy.reshape(numpy.transpose(dataB, (2, 0, 1)), (NBz*NBr*2,))
    r.bwBz = bwBz
    r.bwBr = bwBr
    r.NBz = NBz
    r.NBr = NBr
    r.name = name
    r.okay = True
    smds.Task.intensCache[r.name] = r
    return r    
  fromData = staticmethod(fromData)
  def fromFunc(name, NAz, NAr, bwAz, bwAr, funcA, NBz, NBr, bwBz, bwBr, funcB) :
    """smds.Task.PotProfile.fromFunc( name, NAz, NAr, bwAz, bwAr, funcA(z,r), NBz, NBr, bwBz, bwBr, funcB(z,r) )
    Create an smds.Task.PotProfile with the given size from functions.
    Binwidths are in nm.  The name must be unique.
    The function takes two arrays (Z, R) in nm and returns a pair of arrays
    for the flow in z and r.  Note that Z will be increasing from 0,
    but that in Zone A, this represents moving *up* from the membrane
    and that in Zone B, this represents moving *down* from the top of cap."""
    if name in intensCache :
      raise ValueError, "Intens <%s> already exists." % name
    from numpy import transpose, repeat, reshape
    Za = reshape(transpose(repeat([range(NAz)], NAr, 0), (1,0)),
                 (NAz*NAr,))*float(bwAz)
    Ra = reshape(repeat([range(NAr)], NAz, 0), (NAz*NAr,))*float(bwAr)
    Zb = reshape(transpose(repeat([range(NBz)], NBr, 0), (1,0)),
                 (NBz*NBr,))*float(bwBz)
    Rb = reshape(repeat([range(NBr)], NBz, 0), (NBz*NBr,))*float(bwBr)
    dataA = transpose(reshape(funcA(Za,Ra), (2, NAz, NAr)), (1, 2, 0))
    dataB = transpose(reshape(funcB(Zb,Rb), (2, NBz, NBr)), (1, 2, 0))
    return PotProfile.fromData(name, bwAz, bwAr, dataA, bwBz, bwBr, dataB)
  fromFunc = staticmethod(fromFunc)

class Run :
  """Run-time information."""
  def __init__(self, numSegments = 1) :
    self.host = socket.gethostname()
    self.date = time.strftime("%c")
    self.numSegments = numSegments
  def setTime(self) :
    self.date = time.strftime("%c")
  def incSegments(self, num = 1) :
    self.numSegments += num
  def toXML(self, doc) :
    """Returns this Run as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("Run")
    n.setAttribute("host", self.host)
    n.setAttribute("date", self.date)
    n.setAttribute("numSegments", str(self.numSegments))
    return n
  def fromXML(n) :
    """Returns a new smds.Task.Run created from the xml.dom node n."""
    if (n == None or n.tagName != 'Run') : return None
    r = Run()
    r.host = n.getAttribute('host')
    r.date = n.getAttribute('date')
    r.numSegments = int(n.getAttribute('numSegments'))
    return r
  fromXML = staticmethod(fromXML)

class RTR :
  def __init__(self) :
    self.data = None
  def toXML(self, doc) :
    """Returns this RTR as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("RTR")
    if self.data != None :
      t = ''
      for i in range(len(self.data)) :
        t += "\n" + str(self.data[i])
      n.appendChild(doc.createTextNode(t))
    return n
  def fromXML(n) :
    """Returns a new smds.Task.RTR created from the xml.dom node n."""
    if (n == None or n.tagName != 'RTR') : return None
    r = RTR()
    d = re.findall("\d+", smds.strFromTxtNodes(n.childNodes))
    r.data = numpy.zeros(len(d), numpy.int16)
    for i in range(len(d)) :
      r.data[i] = int(d[i])
    return r
  fromXML = staticmethod(fromXML)

class Task :
  """A simulation task."""
  def __init__(self, ID = None, p = None, rtr = None, anal = None) :
    if ID == None : self.ID = time.strftime('%g%m%d-%H%M%S')
    else          : self.ID = ID
    if p == None  : self.p = smds.Params.c3D()
    else          : self.p = p.copy()
    if p : self.rem = p.getNumBins()
    else : self.rem = 0
    self.completed = 0
    self.solo = True
    self.run = None
    if rtr : self.rtr = RTR()
    else   : self.rtr = None
    self.results = None
    self.anal = []
    if anal :
      for a in anal :
        if isinstance(a, smds.Analysis.Analysis) : self.anal.append(a.copy())
        else :
          smds.msg("Rejecting unknown object %s from analysis list of <%s>"
          		% (str(a), self.ID), smds.MSG_WARN)
    if self.p.okay() : self.status = Status_Ready
    else             : self.status = Status_Error
  def copy(self) :
    return copy.deepcopy(self)
  def toXML(self, doc) :
    """Returns this Task as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("Task")
    n.setAttribute("Core", self.p.getCore().getName())
    n.setAttribute("ID", self.ID)
    n.appendChild(self.p.toXML(doc))
    if self.run != None : n.appendChild(self.run.toXML(doc))
    if self.results != None: n.appendChild(self.results.toXML(doc))
    if self.rtr != None : n.appendChild(self.rtr.toXML(doc))
    for a in self.anal :
      n.appendChild(a.toXML(doc))
    return n
  def fromXML(n) :
    """Returns a new smds.Task.Task created from the xml.dom node n."""
    if (n == None or n.tagName != 'Task') : return None
    core = smds.Cores.CoreTypes[n.getAttribute('Core')]
    r = Task(n.getAttribute('ID'))
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE :
        if c.tagName == 'Params' : r.p = core.getParamsType().fromXML(c)
        if c.tagName == 'Results' :
          r.results = core.getResultsType().fromXML(c)
        if c.tagName == 'Run'  : r.run = Run.fromXML(c)
        if c.tagName == 'RTR'  : r.rtr = RTR.fromXML(c)
        if c.tagName == 'ACF'  : r.anal.append(smds.Analysis.ACF.fromXML(c))
        if c.tagName == 'SEDH' : r.anal.append(smds.Analysis.SEDH.fromXML(c))
        if c.tagName == 'Fit'  : r.anal.append(smds.Analysis.Fit.fromXML(c))
        if c.tagName == 'Analysis' :
          r.anal.append(smds.Analysis.Analysis.fromXML(c))
    if r.run :
      r.rem = 0
      r.status = Status_Done
    else :
      r.rem = r.p.getNumBins()
      if r.p.okay() : r.status = Status_Ready
      else          : r.status = Status_Error
    r.completed = 0
    return r
  fromXML = staticmethod(fromXML)
  def strip(self, minimal) :
    """task.strip(minimal) -> smds.Task
    Returns a stripped version of the task for transmitting over the
    network.  If minimal is true, only paramaters are stripped, otherwise
    all data is stripped.
"""
    # In principle, this is for getting rid of the Intensity profile (or any
    # future large set of parameters that need only be transmited across
    # the network once.  Set minimal to False when transmiting to a remote 
    # processor (to get rid of accumulating results from other processors)
    # and to True when returning data from a remote processor.
    r = Task()
    r.__dict__ = self.__dict__.copy()
    r.p = r.p.strip()
    if not minimal :
      if r.rtr : r.rtr = smds.Task.RTR()
      r.anal = [ a.strip() for a in r.anal ]
      r.run = None
      r.results = None
    return r
