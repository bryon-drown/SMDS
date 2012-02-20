# smds.Results module
#  This module contains the class definitions for simulation results.
#  That is, data generated directly from the simulation cores.
#  $Id: Results.py,v 1.9 2009/09/08 21:27:22 mculbert Exp $

import smds, numpy
from cPickle import dumps, loads

class Base :
  def __init__ (self) :
    self.counts = []		# actual number of each species
    self.avg_I = 0		# average Intensity (phot/ms)
    self.length = 0		# actual time completed (s)
  def cat(self, r) :
    if self.__class__ != r.__class__ :
      raise ValueError, "Cannot concatenate results of different types."
    if len(self.counts) != len(r.counts) :
      raise ValueError, 	\
      	"Cannot concatenate results that have different numbers of species."
    for i in range(len(self.counts)) :
      self.counts[i] += r.counts[i]
    self.avg_I = (self.avg_I*self.length + r.avg_I*r.length) 	\
    		/ (self.length + r.length)
    self.length += r.length
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Results")
    n.appendChild(smds.createXMLvalue(doc, "Length", str(self.length)))
    n.appendChild(smds.createXMLvalue(doc, "Counts", 
    		', '.join([str(x) for x in self.counts]) ))
    n.appendChild(smds.createXMLvalue(doc, "Avg_I", str(self.avg_I)))    
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Results created from the xml.dom node n."""
    if (n == None or n.tagName != 'Results') : return None
    r = cls()
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'Avg_I' :
      self.avg_I = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'Length' :
      self.length = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'Counts' :
      d = smds.strFromTxtNodes(c.childNodes)
      self.counts = [ int(x) for x in d.split(', ') ]
      return True
    return False		# value not found

class Base_phot :
  def __init__ (self) :
    self.bleached = []		# num of each species that bleached
  def cat(self, r) :
    for i in range(len(self.bleached)) :
      self.bleached[i] += r.bleached[i]
  def toXML(self, doc, n) :
    "Appends these Results to node n for insertion into xml.dom.Document doc."
    n.appendChild(smds.createXMLvalue(doc, "Bleached",
    		', '.join([ str(x) for x in self.bleached ]) ))
  def interpretXMLvalue(self, c) :
    if c.getAttribute('Name') == 'Bleached' :
      d = smds.strFromTxtNodes(c.childNodes)
      self.counts = [ int(x) for x in d.split(', ') ]
      return True
    return False

class Base_vesicle :
  def __init__ (self) :
    self.vesicleCounts = []	# number of each vesicle type
    self.freeCounts = []	# number of free molecs of each species
  def cat(self, r) :
    for i in range(len(self.vesicleCounts)) :
      self.vesicleCounts[i] += r.vesicleCounts[i]
    for i in range(len(self.freeCounts)) :
      self.freeCounts[i] += r.freeCounts[i]
  def toXML(self, doc, n) :
    "Apends these Results to node n for insertion into xml.dom.Document doc."
    n.appendChild(smds.createXMLvalue(doc, "vesicleCounts",
    		', '.join([ str(x) for x in self.vesicleCounts ]) ))
    n.appendChild(smds.createXMLvalue(doc, "freeCounts",
    		', '.join([ str(x) for x in self.freeCounts ]) ))
  def interpretXMLvalue(self, c) :
    if c.getAttribute('Name') == 'vesicleCounts' :
      d = smds.strFromTxtNodes(c.childNodes)
      self.vesicleCounts = [ int(x) for x in d.split(', ') ]
      return True
    if c.getAttribute('Name') == 'freeCounts' :
      d = smds.strFromTxtNodes(c.childNodes)
      self.freeCounts = [ int(x) for x in d.split(', ') ]
      return True
    return False

class Base_ccd:
  def __init__ (self) :
    self.video = numpy.zeros(0)
  def toXML(self, doc, n) :
    "Apends these Results to node n for insertion into xml.dom.Document doc."
    n.appendChild(smds.createXMLvalue(doc, "video",
    		dumps(self.video), True ))
  def interpretXMLvalue(self, c) :
    if c.getAttribute('Name') == 'video' :
      d = smds.strFromTxtNodes(c.childNodes, True)
      self.video = loads(str(d))
      return True

class c2D (Base) :
  pass

class c3D (Base) :
  pass

class c2Di (c2D) :
  pass

class c3Di (c3D) :
  pass

class c2Dph (c2D, Base_phot) :
  def __init__ (self) :
    c2D.__init__(self)
    Base_phot.__init__(self)
  def cat(self, r) :
    c2D.cat(self, r)
    Base_phot.cat(self, r)
  def toXML(self, doc) :
    n = c2D.toXML(self, doc)
    Base_phot.toXML(self, doc, n)
    return n
  def interpretXMLvalue(self, c) :
    return ( c2D.interpretXMLvalue(self, c) or
             Base_phot.interpretXMLvalue(self, c) )

class c3Dph (c3D, Base_phot) :
  def __init__ (self) :
    c3D.__init__(self)
    Base_phot.__init__(self)
  def cat(self, r) :
    c3D.cat(self, r)
    Base_phot.cat(self, r)
  def toXML(self, doc) :
    n = c3D.toXML(self, doc)
    Base_phot.toXML(self, doc, n)
    return n
  def interpretXMLvalue(self, c) :
    return ( c3D.interpretXMLvalue(self, c) or
             Base_phot.interpretXMLvalue(self, c) )

class c2Dphi (c2Dph) :
  pass

class c3Dphi (c3Dph) :
  pass

class c3Df (c3D) :
  pass

class c3Db (c3D) :
  pass

class c3Dit (c3Di) :
  pass

class c3Dv (c3D, Base_vesicle) :
  def __init__ (self) :
    c3D.__init__(self)
    Base_vesicle.__init__(self)
  def cat(self, r) :
    c3D.cat(self, r)
    Base_vesicle.cat(self, r)
  def toXML(self, doc) :
    n = c3D.toXML(self, doc)
    Base_vesicle.toXML(self, doc, n)
    return n
  def interpretXMLvalue(self, c) :
    return ( c3D.interpretXMLvalue(self, c) or
             Base_vesicle.interpretXMLvalue(self, c) )

class c3Dv_static (c3Dv) :
  pass

class c3Dvph (c3Dv, Base_phot) :
  def __init__ (self) :
    c3Dv.__init__(self)
    Base_phot.__init__(self)
  def cat(self, r) :
    c3Dv.cat(self, r)
    Base_phot.cat(self, r)
  def toXML(self, doc) :
    n = c3Dv.toXML(self, doc)
    Base_phot.toXML(self, doc, n)
    return n
  def interpretXMLvalue(self, c) :
    return ( c3Dv.interpretXMLvalue(self, c) or
             Base_phot.interpretXMLvalue(self, c) )

class c2D_ccd (c2D, Base_ccd) :
  def __init__ (self) :
    c2D.__init__(self)
    Base_ccd.__init__(self)
  def cat(self, r) :
    raise "c3D_ccd cannot be split."
  def toXML(self, doc) :
    n = c2D.toXML(self, doc)
    Base_ccd.toXML(self, doc, n)
    return n
  def interpretXMLvalue(self, c) :
    return ( c2D.interpretXMLvalue(self, c) or
             Base_ccd.interpretXMLvalue(self, c) )

class c2D_sp (c2D) :
  pass

class cAHL :
  def __init__ (self) :
    self.num = 0		# actual number of simulations completed
    self.hits = 0
    self.misses = 0
    self.hitTimeDist_t = numpy.array([]) # times for time to hist dist (s)
    self.hitTimeDist_count = numpy.array([]) # counts for time to hist dist
  def cat(self, r) :
    if self.__class__ != r.__class__ :
      raise ValueError, "Cannot concatenate results of different types."
    self.num += r.num
    self.hits += r.hits
    self.misses += r.misses
    # Merge time to hist distribution
    bins = len(self.hitTimeDist_t)+len(r.hitTimeDist_t)
    t = [0]*bins
    counts = [0]*bins
    bins = i = j = 0
    while i < len(self.hitTimeDist_t) and j < len(r.hitTimeDist_t) :
      if self.hitTimeDist_t[i] == r.hitTimeDist_t[j] :
        t[bins] = self.hitTimeDist_t[i]
        counts[bins] = self.hitTimeDist_count[i] + r.hitTimeDist_count[j]
        i += 1
        j += 1
        bins += 1
      elif self.hitTimeDist_t[i] < r.hitTimeDist_t[j] :  # select "self"
        t[bins] = self.hitTimeDist_t[i]
        counts[bins] = self.hitTimeDist_count[i]
        i += 1
        bins += 1
      else :  # select "r"
        t[bins] = r.hitTimeDist_t[j]
        counts[bins] = r.hitTimeDist_count[j]
        j += 1
        bins += 1
    if i < len(self.hitTimeDist_t) :
      t[bins:(bins+len(self.hitTimeDist_t)-i)] = self.hitTimeDist_t[i:]
      counts[bins:(bins+len(self.hitTimeDist_t)-i)] = self.hitTimeDist_count[i:]
      bins += len(self.hitTimeDist_t)-i
    if j < len(r.hitTimeDist_t) :
      t[bins:(bins+len(r.hitTimeDist_t)-j)] = r.hitTimeDist_t[j:]
      counts[bins:(bins+len(r.hitTimeDist_t)-j)] = r.hitTimeDist_count[j:]
      bins += len(r.hitTimeDist_t)-j
    self.hitTimeDist_t = t[:bins]
    self.hitTimeDist_count = counts[:bins]
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Results")
    n.appendChild(smds.createXMLvalue(doc, "Num", str(self.num)))
    n.appendChild(smds.createXMLvalue(doc, "Hits", str(self.hits)))
    n.appendChild(smds.createXMLvalue(doc, "Misses", str(self.misses)))
    n.appendChild(smds.createXMLvalue(doc, "TimeToHitDist", 
    		'\n'.join(["%g,%d" % (self.hitTimeDist_t[i],
    		                      self.hitTimeDist_count[i])
    		           for i in range(len(self.hitTimeDist_t))]) ))
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Results created from the xml.dom node n."""
    if (n == None or n.tagName != 'Results') : return None
    r = cls()
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'Num' :
      self.num = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'Hits' :
      self.hits = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'Misses' :
      self.misses = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'TimeToHitDist' :
      import re
      d = smds.strFromTxtNodes(c.childNodes)
      d = re.sub("\r", "\n", d)
      d = re.sub("[ \t]", '', d)
      d = re.sub("\n{2,}", "\n", d)
      if d == '' : data = []
      else :       data = d.split("\n")
      if len(data) and data[0] == '' : data = data[1:]
      if len(data) and data[-1] == '' : data = data[:-1]
      self.hitTimeDist_t = []
      self.hitTimeDist_count = []
      if len(data) :
        for (t, count) in [ p.split(",") for p in data ] :
          self.hitTimeDist_t.append(float(t))
          self.hitTimeDist_count.append(int(count))
      return True
    return False		# value not found

class cAHL_fp (cAHL) :
  pass

class cAHL_lef :
  def __init__ (self) :
    self.length = 0		# actual number of simulations completed
    self.hits = 0
  def cat(self, r) :
    if self.__class__ != r.__class__ :
      raise ValueError, "Cannot concatenate results of different types."
    self.length += r.length
    self.hits += r.hits
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Results")
    n.appendChild(smds.createXMLvalue(doc, "Length", str(self.length)))
    n.appendChild(smds.createXMLvalue(doc, "Hits", str(self.hits)))
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Results created from the xml.dom node n."""
    if (n == None or n.tagName != 'Results') : return None
    r = cls()
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'Length' :
      self.length = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'Hits' :
      self.hits = int(smds.strFromTxtNodes(c.childNodes))
      return True
    return False		# value not found

class cAHL_elef (cAHL_lef) :
  pass

class cAHL_chelef (cAHL_elef) :
  pass

class Base_AHL_1Tape :
  def __init__ (self) :
    self.tape_t = numpy.array([])
    self.tape_x = numpy.array([])
    self.tape_y = numpy.array([])
    self.tape_z = numpy.array([])
  def cat(self, r) :
    if self.__class__ != r.__class__ :
      raise ValueError, "Cannot concatenate results of different types."
    # Merge tape
    self.tape_t = numpy.concatenate((self.tape_t,r.tape_t))
    self.tape_x = numpy.concatenate((self.tape_x,r.tape_x))
    self.tape_y = numpy.concatenate((self.tape_y,r.tape_y))
    self.tape_z = numpy.concatenate((self.tape_z,r.tape_z))
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Results")
    n.appendChild(smds.createXMLvalue(doc, "Tape",
    		'\n'.join(["%g,%g,%g,%g" % (self.tape_t[i],self.tape_x[i],self.tape_y[i],self.tape_z[i])
    		           for i in range(len(self.tape_x))]) ))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'Tape' :
      import re
      d = smds.strFromTxtNodes(c.childNodes)
      d = re.sub("\r", "\n", d) # replace linefeed with newline
      d = re.sub("[ \t]", '', d)    # remove all tabs
      d = re.sub("\n{2,}", "\n", d) # replace all repeating newlines with one
      if d == '' : data = []
      else :       data = d.split("\n")
      if len(data) and data[0] == '' : data = data[1:]
      if len(data) and data[-1] == '' : data = data[:-1]
      self.tape_t = []
      self.tape_x = []
      self.tape_y = []
      self.tape_z = []
      if len(data) :
        for (t, x, y, z) in [ p.split(",") for p in data ] :
          self.tape_t.append(float(t))
          self.tape_x.append(float(x))
          self.tape_y.append(float(y))
          self.tape_z.append(float(z))
      return True
    return False		# value not found

class Base_AHL_2Tape :
  def __init__ (self) :
    self.hitTape_t = numpy.array([])
    self.hitTape_x = numpy.array([])
    self.hitTape_y = numpy.array([])
    self.hitTape_z = numpy.array([])
    self.missTape_t = numpy.array([])
    self.missTape_x = numpy.array([])
    self.missTape_y = numpy.array([])
    self.missTape_z = numpy.array([])
  def cat(self, r) :
    if self.__class__ != r.__class__ :
      raise ValueError, "Cannot concatenate results of different types."
    # Merge tapes
    self.hitTape_t = numpy.concatenate((self.hitTape_t,r.hitTape_t))
    self.hitTape_x = numpy.concatenate((self.hitTape_x,r.hitTape_x))
    self.hitTape_y = numpy.concatenate((self.hitTape_y,r.hitTape_y))
    self.hitTape_z = numpy.concatenate((self.hitTape_z,r.hitTape_z))
    self.missTape_t = numpy.concatenate((self.missTape_t,r.missTape_t))
    self.missTape_x = numpy.concatenate((self.missTape_x,r.missTape_x))
    self.missTape_y = numpy.concatenate((self.missTape_y,r.missTape_y))
    self.missTape_z = numpy.concatenate((self.missTape_z,r.missTape_z))
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n.appendChild(smds.createXMLvalue(doc, "Hit Tape",
    		'\n'.join(["%g,%g,%g,%g" % (self.hitTape_t[i],self.hitTape_x[i],self.hitTape_y[i],self.hitTape_z[i])
    		           for i in range(len(self.hitTape_x))]) ))
    n.appendChild(smds.createXMLvalue(doc, "Miss Tape",
    		'\n'.join(["%g,%g,%g,%g" % (self.missTape_t[i],self.missTape_x[i],self.missTape_y[i],self.missTape_z[i])
    		           for i in range(len(self.missTape_x))]) ))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'Hit Tape' :
      import re
      d = smds.strFromTxtNodes(c.childNodes)
      d = re.sub("\r", "\n", d) # replace linefeed with newline
      d = re.sub("[ \t]", '', d)    # remove all tabs
      d = re.sub("\n{2,}", "\n", d) # replace all repeating newlines with one
      if d == '' : data = []
      else :       data = d.split("\n")
      if len(data) and data[0] == '' : data = data[1:]
      if len(data) and data[-1] == '' : data = data[:-1]
      self.hitTape_t = []
      self.hitTape_x = []
      self.hitTape_y = []
      self.hitTape_z = []
      if len(data) :
        for (t, x, y, z) in [ p.split(",") for p in data ] :
          self.hitTape_t.append(float(t))
          self.hitTape_x.append(float(x))
          self.hitTape_y.append(float(y))
          self.hitTape_z.append(float(z))
      return True
    if name == 'Miss Tape' :
      import re
      d = smds.strFromTxtNodes(c.childNodes)
      d = re.sub("\r", "\n", d) # replace linefeed with newline
      d = re.sub("[ \t]", '', d)    # remove all tabs
      d = re.sub("\n{2,}", "\n", d) # replace all repeating newlines with one
      if d == '' : data = []
      else :       data = d.split("\n")
      if len(data) and data[0] == '' : data = data[1:]
      if len(data) and data[-1] == '' : data = data[:-1]
      self.missTape_t = []
      self.missTape_x = []
      self.missTape_y = []
      self.missTape_z = []
      if len(data) :
        for (t, x, y, z) in [ p.split(",") for p in data ] :
          self.missTape_t.append(float(t))
          self.missTape_x.append(float(x))
          self.missTape_y.append(float(y))
          self.missTape_z.append(float(z))
      return True
    return False		# value not found


class cAHL_sdx (cAHL) :
    pass

class cAHL_sdxp (cAHL_sdx) :
    pass

class cAHL_sphdx (cAHL_sdx) :
    pass

class cAHL_sdxeo (cAHL_sdx) :
    pass

class cAHL_sphdxeo (cAHL_sdxeo):
    pass

class cAHL_sdxt (cAHL, Base_AHL_1Tape) :
  def __init__ (self) :
    cAHL.__init__(self)
    Base_AHL_1Tape.__init__(self)
  def cat(self, r) :
    cAHL.cat(self, r)
    Base_AHL_1Tape.cat(self, r)
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n = cAHL.toXML(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "Tape",
    		'\n'.join(["%g,%g,%g,%g" % (self.tape_t[i],self.tape_x[i],self.tape_y[i],self.tape_z[i])
    		           for i in range(len(self.tape_x))]) ))
    return n
  def interpretXMLvalue(self, c) :
    if cAHL.interpretXMLvalue(self, c) : return True
    if Base_AHL_1Tape.interpretXMLvalue(self, c) : return True
    return False		# value not found

class cAHL_sdxeot (cAHL_sdxt) :
    """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Tape Record"""
    pass

class cAHL_sphdxeot (cAHL_sdxeot) :
    """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Spherical Particle, Tape Record"""
    pass

class cAHL_sdx2t (cAHL, Base_AHL_2Tape) :
  def __init__ (self) :
    cAHL.__init__(self)
    Base_AHL_2Tape.__init__(self)
  def cat(self, r) :
    cAHL.cat(self, r)
    Base_AHL_2Tape.cat(self, r)
  def toXML(self, doc) :
    "Returns these Results as a node for insertion into xml.dom.Document doc."
    n = cAHL.toXML(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "Hit Tape",
    		'\n'.join(["%g,%g,%g,%g" % (self.hitTape_t[i],self.hitTape_x[i],self.hitTape_y[i],self.hitTape_z[i])
    		           for i in range(len(self.hitTape_x))]) ))
    n.appendChild(smds.createXMLvalue(doc, "Miss Tape",
    		'\n'.join(["%g,%g,%g,%g" % (self.missTape_t[i],self.missTape_x[i],self.missTape_y[i],self.missTape_z[i])
    		           for i in range(len(self.missTape_x))]) ))
    return n
  def interpretXMLvalue(self, c) :
    if cAHL.interpretXMLvalue(self, c) : return True
    if Base_AHL_2Tape.interpretXMLvalue(self, c) : return True
    return False		# value not found

class cAHL_sphdx2t (cAHL_sdx2t) :
    """Simulation results, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distribution, Hit & Miss Tape Records, Spherical Particle"""
    pass

class cAHL_sdxeo2t (cAHL_sdx2t) :
  """Simulation results, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Hit & Miss Tape Records"""
  pass

class cAHL_sphdxeo2t (cAHL_sdx2t) :
  """Simulation results, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Hit & Miss Tape Records, Spherical Particle"""
  pass