# smds.Params module
#  This module contains the class definitions for the Parameter structures
#  $Id: Params.py,v 1.11 2009/06/16 22:36:11 mculbert Exp $

from warnings import warn
import copy, math, xml.dom, smds

class Species :
  """Species parameters"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    self.D = 1e-7				# cm^2/s
    self.Imax = 1600			# phot/ms
    self.Fraction = 1.0
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def toXML(self, doc) :
    "Returns this Species as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Species")
    n.setAttribute("D", '%g'%(self.D,))
    n.setAttribute("Fraction", str(self.Fraction))
    n.setAttribute("Imax", str(self.Imax))
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Params.Species created from the xml.dom node n."""
    if (n == None or n.tagName != 'Species') : return None
    r = cls()
    r.interpretXMLhead(n)
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def interpretXMLhead(self, n) :
    self.D = float(n.getAttribute("D"))
    self.Imax = float(n.getAttribute("Imax"))
    self.Fraction = float(n.getAttribute("Fraction"))
  def interpretXMLvalue(self, c) :
    return False
  def strip(self) :
    return self

class Species_phot (Species) :
  """Species with photobleaching"""
  def __init__(self, p = None) :
    self.tolerance = 1000				# phot
    Species.__init__(self, p)
  def toXML(self, doc) :
    "Returns this Species as a node for insertion into xml.dom.Document doc."
    n = Species.toXML(self, doc)
    n.setAttribute("Tolerance", str(self.tolerance))
    return n
  def interpretXMLhead(self, n) :
    self.tolerance = float(n.getAttribute("Tolerance"))
    Species.interpretXMLhead(self, n)

class Species_blink (Species) :
  """Species with blinking"""
  def __init__(self, p = None) :
    self.k_dark = 10				# s^-1
    self.k_bright = 100				# s^-1
    Species.__init__(self, p)
  def toXML(self, doc) :
    "Returns this Species as a node for insertion into xml.dom.Document doc."
    n = Species.toXML(self, doc)
    n.setAttribute("k_dark", str(self.k_dark))
    n.setAttribute("k_brigt", str(self.k_bright))
    return n
  def interpretXMLhead(self, n) :
    self.k_dark = float(n.getAttribute("k_dark"))
    self.k_bright = float(n.getAttribute("k_bright"))
    Species.interpretXMLhead(self, n)

class Species_triplet (Species) :
  """Species with intersystem crossing"""
  def __init__(self, p = None) :
    self.k_T = 10				# s^-1
    self.k_S0 = 100				# s^-1
    Species.__init__(self, p)
  def toXML(self, doc) :
    "Returns this Species as a node for insertion into xml.dom.Document doc."
    n = Species.toXML(self, doc)
    n.setAttribute("k_T", str(self.k_T))
    n.setAttribute("k_S0", str(self.k_S0))
    return n
  def interpretXMLhead(self, n) :
    self.k_T = float(n.getAttribute("k_T"))
    self.k_S0 = float(n.getAttribute("k_S0"))
    Species.interpretXMLhead(self, n)

class Species_states (Species) :
  """Species with alternating characteristics"""
  def __init__(self, p = None) :
    self.k_species = 10				# s^-1
    self.altSpecies = 0
    Species.__init__(self, p)
  def toXML(self, doc) :
    "Returns this Species as a node for insertion into xml.dom.Document doc."
    n = Species.toXML(self, doc)
    n.setAttribute("k_species", str(self.k_species))
    n.setAttribute("altSpecies", str(self.altSpecies))
    return n
  def interpretXMLhead(self, n) :
    self.k_species = float(n.getAttribute("k_species"))
    self.altSpecies = int(float(n.getAttribute("altSpecies"))+0.5)
    Species.interpretXMLhead(self, n)

class Vesicle :
  """Vesicle parameters"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    self.D = 1e-7				# cm^2/s
    self.Fraction = 1.0
    self.numMolecs = 10.0			# avg num of molecules
    self.size = 1.0				# avg radius, um
    self.size_spread = 0.2			# std dev of radius, um
    self.composition = [ 1.0 ]			# fraction of each species
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def toXML(self, doc) :
    "Returns this Vesicle as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Vesicle")
    n.setAttribute("D", '%g'%(self.D,))
    n.setAttribute("Fraction", str(self.Fraction))
    n.setAttribute("numMolecs", str(self.numMolecs))
    n.setAttribute("size", str(self.size))
    n.setAttribute("size_spread", str(self.size_spread))
    n.setAttribute("composition",
    		', '.join([ str(x) for x in self.composition ]) )
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Params.Vesicle created from the xml.dom node n."""
    if (n == None or n.tagName != 'Vesicle') : return None
    r = cls()
    r.interpretXMLhead(n)
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def interpretXMLhead(self, n) :
    self.D = float(n.getAttribute("D"))
    self.Fraction = float(n.getAttribute("Fraction"))
    self.numMolecs = float(n.getAttribute("numMolecs"))
    self.size = float(n.getAttribute("size"))
    self.size_spread = float(n.getAttribute("size_spread"))
    self.composition = [ float(x) for x in
    			n.getAttribute("composition").split(', ') ]
  def interpretXMLvalue(self, c) :
    return False
  def strip(self) :
    return self

class Base:
  """Simulation parameters, Base functions"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    self.dT = 500				# nm
    self.binWidth = 0.1				# ms
    self.numMolecs = 48
    self.dur = 100				# s
    self.bkg = 7.15				# phot/ms
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = Species
  Unity = False
  def copy(self) :
    return copy.deepcopy(self)
  def getCore(cls) :
    return smds.Cores.paramsToCore[cls]
  getCore = classmethod(getCore)
  def initCore(self) :
    return self.getCore()(self)
  def getNumBins(self) :
    return int(self.dur * 1e3 / self.binWidth + 0.5)
  def toXML(self, doc) :
    "Returns these Params as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Params")
    n.appendChild(self.createXMLlaser(doc))
    n.appendChild(self.createXMLsample(doc))
    n.appendChild(self.createXMLsim(doc))
    ## add extra values with smds.createXMLvalue() here ##
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Params created from the xml.dom node n."""
    if (n == None or n.tagName != 'Params') : return None
    r = cls()
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE :
        if c.tagName == 'Laser' : r.interpretXMLlaser(c)
        if c.tagName == 'Sample' : r.interpretXMLsample(c)
        if c.tagName == 'Sim' : r.interpretXMLsim(c)
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def createXMLlaser(self, doc) :
    n = doc.createElement("Laser")
    n.setAttribute("bkg", str(self.bkg))
    return n
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.setAttribute("dT", str(self.dT))
    n.setAttribute("binWidth", str(self.binWidth))
    n.setAttribute("numMolecs", str(self.numMolecs))
    n.setAttribute("dur", str(self.dur))
    return n
  def interpretXMLlaser(self, n) :
    self.bkg = float(n.getAttribute("bkg"))
  def interpretXMLsim(self, n) :
    self.dT = int(n.getAttribute("dT"))
    self.binWidth = float(n.getAttribute("binWidth"))
    self.numMolecs = int(n.getAttribute("numMolecs"))
    self.dur = int(n.getAttribute("dur"))
  def interpretXMLvalue(self, c) :
    return False
  def okay(self) :
    return True
  def strip(self) :
    return self

class Base_intens :
  """Simulation parameters with Intensity profile"""
  def __init__(self, p = None) :
    self.intens = smds.Task.Intens("")				# filename
    if (p and p.has_key('intens')) : self.intens = p['intens']
  def copy(self) :
    intens = self.intens
    self.intens = None
    r = copy.deepcopy(self)
    self.intens = intens
    r.intens = intens
    return r
  def createXMLlaser(self, doc, n) :
    n.setAttribute("IntensityProfile", self.intens.name)
  def interpretXMLlaser(self, n) :
    self.intens = smds.Task.Intens(n.getAttribute("IntensityProfile"))
  def okay(self) :
    return self.intens.okay
  def strip(self) :
    r = self.__class__()
    r.__dict__ = self.__dict__.copy()
    r.intens = r.intens.strip()
    return r

class Base_pot :
  """Simulation parameters with Potential profile"""
  def __init__(self, p = None) :
    self.pot = smds.Task.PotProfile("")				# filename
    if (p and p.has_key('pot')) : self.pot = p['pot']
  def createXMLchannel(self, doc) :
    n = doc.createElement("Channel")
    n.appendChild(smds.createXMLvalue(doc, "PotProfile", self.pot.name))
    n.appendChild(smds.createXMLvalue(doc, "potA",
                                      "%g" % self.potA))
    n.appendChild(smds.createXMLvalue(doc, "potB",
                                      "%g" % self.potB))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'PotProfile' :
      self.potExtent = smds.Task.PotProfile(smds.strFromTxtNodes(c.childNodes).strip())
      return True
    if name == 'potA' :
      self.potA = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potB' :
      self.potB = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False
  def copy(self) :
    pot = self.pot
    self.pot = None
    r = copy.deepcopy(self)
    self.pot = pot
    r.pot = pot
    return r
  def okay(self) :
    return self.pot.okay
  def strip(self) :
    r = self.__class__()
    r.__dict__ = self.__dict__.copy()
    r.pot = r.pot.strip()
    return r

class Base_triplet :
  """Simulation parameters with CEF profile"""
  def __init__(self, p = None) :
    self.collectionEfficiency = 1.0
    self.cef = smds.Task.Intens("")				# filename
    if (p and p.has_key('cef')) : self.cef = p['cef']
  def copy(self) :
    intens = self.intens
    cef = self.cef
    self.intens = None
    self.cef = None
    r = copy.deepcopy(self)
    self.intens = intens
    self.cef = cef
    r.intens = intens
    r.cef = cef
    return r
  def createXMLlaser(self, doc, n) :
    n.setAttribute("CollectionEfficiencyProfile", self.cef.name)
    n.setAttribute("collectionEfficiency", 
			"%g" % self.collectionEfficiency)
  def interpretXMLlaser(self, n) :
    self.cef = smds.Task.Intens(n.getAttribute("CollectionEfficienyProfile"))
    self.collectionEfficiency = float(n.getAttribute("collectionEfficiency"))
  def okay(self) :
    return self.cef.okay
  def strip(self) :
    r = self.__class__()
    r.__dict__ = self.__dict__.copy()
    r.cef = r.cef.strip()
    r.intens = r.intens.strip()
    return r

class Base_ccd :
  """Simulation parameters for CCD detection"""
  def __init__(self, p = None) :
    self.video_framerate = 1.25		# ms/frame
    self.video_Nx = 40			# pixels
    self.video_Ny = 40			# pixels
    self.video_sizeX = 125		# nm
    self.video_sizeY = 125		# nm
  def toXML(self, doc, n) :
    n.appendChild(smds.createXMLvalue(doc, "video_framerate",
                                      "%g"%self.video_framerate))
    n.appendChild(smds.createXMLvalue(doc, "video_Nx", 
						"%d"%self.video_Nx))
    n.appendChild(smds.createXMLvalue(doc, "video_Ny", 
						"%d"%self.video_Ny))
    n.appendChild(smds.createXMLvalue(doc, "video_sizeX", 
						"%d"%self.video_sizeX))
    n.appendChild(smds.createXMLvalue(doc, "video_sizeY", 
						"%d"%self.video_sizeY))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'video_framerate' :
      self.video_framerate = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'video_Nx' :
      self.video_Nx = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'video_Ny' :
      self.video_Ny = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'video_sizeX' :
      self.video_sizeX = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'video_sizeY' :
      self.video_sizeY = int(smds.strFromTxtNodes(c.childNodes))
      return True
    return False		# value not found

class c2D (Base) :
  """Basic 2D simulation"""
  def __init__(self, p = None) :
    self.sample = []				# Params.Species()
    self.concentration = 1.0			# molec/um^2
    self.radius = 264				# nm
    self.threshold = 750			# nm
    Base.__init__(self, p)
    self.sample = list(self.sample)
    for i in range(len(self.sample)) :
      if type(self.sample[i]) == type({}) :
        self.sample[i] = self.SpeciesType(self.sample[i])
#    self.size = math.sqrt(self.numMolecs/self.concentration)
  def createXMLlaser(self, doc) :
    r = Base.createXMLlaser(self, doc)
    r.setAttribute("Radius", str(self.radius))
    r.setAttribute("Threshold", str(self.threshold))
    return r
  def createXMLsample(self, doc) :
    r = doc.createElement("Sample")
    r.setAttribute("Concentration", str(self.concentration))
    for s in self.sample :
      r.appendChild(s.toXML(doc))
    return r
  def interpretXMLlaser(self, n) :
    self.radius = float(n.getAttribute("Radius"))
    self.threshold = int(n.getAttribute("Threshold"))
    Base.interpretXMLlaser(self, n)
  def interpretXMLsample(self, n) :
    self.concentration = float(n.getAttribute("Concentration"))
    self.sample = []
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE and    \
         c.tagName == 'Species' :
        self.sample.append(self.SpeciesType.fromXML(c))

class c3D (c2D) :
  """Basic 3D simulation"""
  def __init__(self, p = None) :
    self.Z = 1000				# nm
    self.threshold_Z = 2000			# nm
    c2D.__init__(self, p)	# concentration is now in molec/um^3
#    self.size = (self.numMolecs/self.concentration)**(1.0/3.0)
  def createXMLlaser(self, doc) :
    r = c2D.createXMLlaser(self, doc)
    r.setAttribute("Z", str(self.Z))
    r.setAttribute("Threshold_Z", str(self.threshold_Z))
    return r
  def interpretXMLlaser(self, n) :
    self.Z = float(n.getAttribute("Z"))
    self.threshold_Z = int(n.getAttribute("Threshold_Z"))
    c2D.interpretXMLlaser(self, n)

class c2Dph (c2D) :
  """2D simulation with photobleaching."""
  SpeciesType = Species_phot

class c3Dph (c3D) :
  """3D simulation with photobleaching."""
  SpeciesType = Species_phot

class c2Di (c2D, Base_intens) :
  """2D simulation with intensity map."""
  def __init__(self, p = None) :
    c2D.__init__(self, p)
    del self.threshold
    Base_intens.__init__(self, p)
  def copy(self) :
    return Base_intens.copy(self)
  def createXMLlaser(self, doc) :
    r = Base.createXMLlaser(self, doc)
    r.setAttribute("Radius", str(self.radius))
    Base_intens.createXMLlaser(self, doc, r)
    return r
  def interpretXMLlaser(self, n) :
    Base_intens.interpretXMLlaser(self, n)
    if n.hasAttribute("Radius") :
      self.radius = float(n.getAttribute("Radius"))
    Base.interpretXMLlaser(self, n)
  def okay(self) :
    return c2D.okay(self) and Base_intens.okay(self)
  def strip(self) :
    return Base_intens.strip(self)

class c3Di (c3D, Base_intens) :
  """3D simulation with intensity map."""
  def __init__(self, p = None) :
    c3D.__init__(self, p)
    del self.threshold
    del self.threshold_Z
    Base_intens.__init__(self, p)
  def copy(self) :
    return Base_intens.copy(self)
  def createXMLlaser(self, doc) :
    r = Base.createXMLlaser(self, doc)
    r.setAttribute("Radius", str(self.radius))
    r.setAttribute("Z", str(self.Z))
    Base_intens.createXMLlaser(self, doc, r)
    return r
  def interpretXMLlaser(self, n) :
    Base_intens.interpretXMLlaser(self, n)
    if n.hasAttribute("Radius") :
      self.radius = float(n.getAttribute("Radius"))
    if n.hasAttribute("Z") :
      self.Z = float(n.getAttribute("Z"))
    Base.interpretXMLlaser(self, n)
  def okay(self) :
    return c3D.okay(self) and Base_intens.okay(self)
  def strip(self) :
    return Base_intens.strip(self)

class c2Dphi (c2Di) :
  """2D simulation with photobleaching and intensity map."""
  SpeciesType = Species_phot

class c3Dphi (c3Di) :
  """3D simulation with photobleaching and intensity map."""
  SpeciesType = Species_phot

class c3Df (c3D) :
  """3D simulation with flow."""
  def __init__(self, p = None) :
    self.flow = 0.0			# m/s
    self.flowLength = 1.0;		# um before center of det vol
    self.drainLength = 1.0;		# um after center of det vol
    c3D.__init__(self, p)
  def toXML(self, doc) :
    n = c3D.toXML(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "flow", "%g"%self.flow))
    n.appendChild(smds.createXMLvalue(doc, "flowLength", 
						"%g"%self.flowLength))
    n.appendChild(smds.createXMLvalue(doc, "drainLength", 
						"%g"%self.drainLength))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'flow' :
      self.flow = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'flowLength' :
      self.flowLength = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'drainLength' :
      self.drainLength = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False		# value not found

class c3Db (c3D) :
  """3D simulation with blinking."""
  SpeciesType = Species_blink

class c3Dit (c3Di, Base_triplet) :
  """3D simulation with intersystem crossing."""
  def __init__(self, p = None) :
    c3Di.__init__(self, p)
    Base_triplet.__init__(self, p)
  def copy(self) :
    return Base_triplet.copy(self)
  def createXMLlaser(self, doc) :
    r = c3Di.createXMLlaser(self, doc)
    Base_triplet.createXMLlaser(self, doc, r)
    return r
  def interpretXMLlaser(self, n) :
    Base_triplet.interpretXMLlaser(self, n)
    c3Di.interpretXMLlaser(self, n)
  def okay(self) :
    return c3Di.okay(self) and Base_triplet.okay(self)
  def strip(self) :
    return Base_triplet.strip(self)
  SpeciesType = Species_triplet

class c3Dv (c3D) :
  """3D simulation with vesicles"""
  def __init__(self, p = None) :
    self.vesicleConcentration = 0.1		# um^-3
    self.numVesicles = 10
    self.vesicles = []
    c3D.__init__(self, p)
    self.vesicles = list(self.vesicles)
    for i in range(len(self.vesicles)) :
      if type(self.vesicles[i]) == type({}) :
        self.vesicles[i] = self.VesicleType(self.vesicles[i])
  VesicleType = Vesicle
  def createXMLsample(self, doc) :
    r = c3D.createXMLsample(self, doc)
    r.setAttribute("vesicleConcentration", str(self.vesicleConcentration))
    for s in self.vesicles :
      r.appendChild(s.toXML(doc))
    return r
  def createXMLsim(self, doc) :
    r = c3D.createXMLsim(self, doc)
    r.setAttribute("numVesicles", str(self.numVesicles))
    return r
  def interpretXMLsample(self, n) :
    self.vesicleConcentration = float(n.getAttribute("vesicleConcentration"))
    self.vesicles = []
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE and    \
         c.tagName == 'Vesicle' :
        self.vesicles.append(self.VesicleType.fromXML(c))
    c3D.interpretXMLsample(self, n)
  def interpretXMLsim(self, n) :
    self.numVesicles = int(n.getAttribute("numVesicles"))
    c3D.interpretXMLsim(self, n)

class c3Dv_static (c3Dv) :
  pass

class c3Dvph (c3Dv) :
  """3D simulation with vesicles and photobleaching."""
  SpeciesType = Species_phot

class c2D_ccd (c2D, Base_ccd) :
  """2D simulation with CCD detection."""
  def __init__(self, p = None) :
    Base_ccd.__init__(self, p)
    c2D.__init__(self, p)
  Unity = True
  def toXML(self, doc) :
    n = c2D.toXML(self, doc)
    Base_ccd.toXML(self, doc, n)
    return n
  def interpretXMLvalue(self, c) :
    return ( c2D.interpretXMLvalue(self, c) or
             Base_ccd.interpretXMLvalue(self, c) )

class c2D_sp (c2D) :
  """2D simulation with Alternating species"""
  SpeciesType = Species_states

class cAHL:
  """Simulation parameters, Alpha-Hemolycin Core"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    self.num = 100				# iterations
    self.dT = 0.001				# ns
    self.binWidth = 1				# ns
    self.timeToMiss = 10			# ms
    self.initHeight = 0				# nm above cap height
    self.D = 3e-6				# cm^2/s
    self.capHeight = 7.				# nm above membrane
    self.capWidth = 3.3				# nm
    self.channelRadius = 1.			# nm
    self.hitPoint = 7.5				# nm below cap height
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = None
  Unity = False
  def copy(self) :
    return copy.deepcopy(self)
  def getCore(cls) :
    return smds.Cores.paramsToCore[cls]
  getCore = classmethod(getCore)
  def initCore(self) :
    return self.getCore()(self)
  def getNumBins(self) :
    return self.num
  def toXML(self, doc) :
    "Returns these Params as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Params")
    n.appendChild(self.createXMLchannel(doc))
    n.appendChild(self.createXMLsim(doc))
    ## add extra values with smds.createXMLvalue() here ##
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Params created from the xml.dom node n."""
    if (n == None or n.tagName != 'Params') : return None
    r = cls()
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE :
        if c.tagName == 'Channel' : r.interpretXMLchannel(c)
        if c.tagName == 'Sim' : r.interpretXMLsim(c)
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def createXMLchannel(self, doc) :
    n = doc.createElement("Channel")
    n.appendChild(smds.createXMLvalue(doc, "capHeight", "%g"%self.capHeight))
    n.appendChild(smds.createXMLvalue(doc, "capWidth", "%g"%self.capWidth))
    n.appendChild(smds.createXMLvalue(doc, "channelRadius",
                                      "%g"%self.channelRadius))
    n.appendChild(smds.createXMLvalue(doc, "hitPoint", "%g"%self.hitPoint))
    return n
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.setAttribute("dT", str(self.dT))
    n.setAttribute("binWidth", str(self.binWidth))
    n.setAttribute("num", str(self.num))
    n.appendChild(smds.createXMLvalue(doc, "timeToMiss", "%g"%self.timeToMiss))
    n.appendChild(smds.createXMLvalue(doc, "initHeight", "%g"%self.initHeight))
    n.appendChild(smds.createXMLvalue(doc, "D", "%g"%self.D))
    return n
  def interpretXMLchannel(self, n) :
    smds.interpretXMLvalues(self, n)
  def interpretXMLsim(self, n) :
    self.dT = float(n.getAttribute("dT"))
    self.binWidth = float(n.getAttribute("binWidth"))
    self.num = int(n.getAttribute("num"))
    smds.interpretXMLvalues(self, n)
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'timeToMiss' :
      self.timeToMiss = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'initHeight' :
      self.initHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'D' :
      self.D = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'capHeight' :
      self.capHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'capWidth' :
      self.capWidth = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'channelRadius' :
      self.channelRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'hitPoint' :
      self.hitPoint = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False
  def okay(self) :
    return True
  def strip(self) :
    return self

class cAHL_fp (cAHL):
  """Simulation parameters, Alpha-Hemolycin Core, Variable step size"""
  pass
  
class cAHL_lef:
  """Simulation parameters, Alpha-Hemolycin Core, Electrophoretic Flow"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    self.dur = 100				# s
    self.macro_dT = 500				# ns
    self.micro_dT = 0.001			# ns
    self.binWidth = 1				# us
    self.numMolecs = 100
    self.concentration = 0.1			# molecs / um^3
    self.D = 3e-6				# cm^2/s
    self.capHeight = 7.				# nm above membrane
    self.capWidth = 3.3				# nm
    self.channelRadius = 1.			# nm
    self.hitPoint = 7.5				# nm below cap height
    self.vestibuleRadius = 3.			# nm
    self.vestibuleHeight = 3.			# nm tall
    self.vestibulePos = 2.			# nm below cap height
    self.potExtent = 10.		# nm, radial from cap height center
    self.potA = 1.				# mV at cap height
    self.potB = 5.				# mV at bottom of vestibule
    self.potC = 100.				# mV at hitpoint
    self.mobility = 15e9			# nm^2/Vs
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = None
  Unity = False
  def copy(self) :
    return copy.deepcopy(self)
  def getCore(cls) :
    return smds.Cores.paramsToCore[cls]
  getCore = classmethod(getCore)
  def initCore(self) :
    return self.getCore()(self)
  def getNumBins(self) :
    return int(self.dur * 1e6 / self.binWidth + 0.5)
  def toXML(self, doc) :
    "Returns these Params as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Params")
    n.appendChild(self.createXMLchannel(doc))
    n.appendChild(self.createXMLsim(doc))
    ## add extra values with smds.createXMLvalue() here ##
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Params created from the xml.dom node n."""
    if (n == None or n.tagName != 'Params') : return None
    r = cls()
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE :
        if c.tagName == 'Channel' : r.interpretXMLchannel(c)
        if c.tagName == 'Sim' : r.interpretXMLsim(c)
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def createXMLchannel(self, doc) :
    n = doc.createElement("Channel")
    n.appendChild(smds.createXMLvalue(doc, "capHeight", "%g"%self.capHeight))
    n.appendChild(smds.createXMLvalue(doc, "capWidth", "%g"%self.capWidth))
    n.appendChild(smds.createXMLvalue(doc, "channelRadius",
                                      "%g"%self.channelRadius))
    n.appendChild(smds.createXMLvalue(doc, "hitPoint", "%g"%self.hitPoint))
    n.appendChild(smds.createXMLvalue(doc, "vestibuleRadius",
                                      "%g" % self.vestibuleRadius))
    n.appendChild(smds.createXMLvalue(doc, "vestibuleHeight",
                                      "%g" % self.vestibuleHeight))
    n.appendChild(smds.createXMLvalue(doc, "vestibulePos",
                                      "%g" % self.vestibulePos))
    n.appendChild(smds.createXMLvalue(doc, "vestibuleRadius",
                                      "%g" % self.vestibuleRadius))
    n.appendChild(smds.createXMLvalue(doc, "potExtent",
                                      "%g" % self.potExtent))
    n.appendChild(smds.createXMLvalue(doc, "potA",
                                      "%g" % self.potA))
    n.appendChild(smds.createXMLvalue(doc, "potB",
                                      "%g" % self.potB))
    n.appendChild(smds.createXMLvalue(doc, "potC",
                                      "%g" % self.potC))
    return n
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.setAttribute("macro_dT", str(self.macro_dT))
    n.setAttribute("micro_dT", str(self.micro_dT))
    n.setAttribute("binWidth", str(self.binWidth))
    n.setAttribute("dur", str(self.dur))
    n.setAttribute("numMolecs", str(self.numMolecs))
    n.appendChild(smds.createXMLvalue(doc, "concentration",
                                      "%g" % self.concentration))
    n.appendChild(smds.createXMLvalue(doc, "D", "%g" % self.D))
    n.appendChild(smds.createXMLvalue(doc, "mobility", "%g" % self.mobility))
    return n
  def interpretXMLchannel(self, n) :
    smds.interpretXMLvalues(self, n)
  def interpretXMLsim(self, n) :
    self.macro_dT = float(n.getAttribute("macro_dT"))
    self.micro_dT = float(n.getAttribute("micro_dT"))
    self.binWidth = float(n.getAttribute("binWidth"))
    self.numMolecs = float(n.getAttribute("numMolecs"))
    self.dur = int(n.getAttribute("dur"))
    smds.interpretXMLvalues(self, n)
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'concentration' :
      self.concentration = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'D' :
      self.D = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'mobility' :
      self.mobility = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'capHeight' :
      self.capHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'capWidth' :
      self.capWidth = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'channelRadius' :
      self.channelRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'hitPoint' :
      self.hitPoint = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'vestibuleRadius' :
      self.vestibuleRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'vestibuleHeight' :
      self.vestibuleHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'vestibulePos' :
      self.vestibulePos = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potExtent' :
      self.potExtent = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potA' :
      self.potA = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potB' :
      self.potB = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potC' :
      self.potC = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False
  def okay(self) :
    return True
  def strip(self) :
    return self

class cAHL_elef (cAHL_lef):
  """Simulation parameters, Alpha-Hemolycin Core, Electrophoretic Flow, Exponential Field"""
  pass

class cAHL_chelef (cAHL_elef):
  """Simulation parameters, Alpha-Hemolycin Core, Electrophoretic Flow, Exponential Field, Point Charges"""
  def __init__(self, p=None) :
    self.chargeRadius = 3.3			# nm
    self.charge = 1.0				# effective, in electrons
    cAHL_elef.__init__(self, p)
  def createXMLchannel(self, doc):
    n = cAHL_elef.createXMLchannel(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "chargeRadius",
                                      "%g" % self.chargeRadius))
    n.appendChild(smds.createXMLvalue(doc, "charge", "%g" % self.charge))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'chargeRadius' :
      self.chargeRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'charge' :
      self.charge = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return cAHL_elef.interpretXMLvalue(self, c)

class Base_AHL_sdx :
  """Simulaton parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distributioon"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    self.num = 100				# iterations
    self.dT = 500				# ns
    self.binWidth = 1				# us
    self.D = 3e-6				# cm^2/s
    self.capHeight = 7.				# nm above membrane
    self.capWidth = 3.3				# nm
    self.channelRadius = 1.			# nm
    self.hitPoint = 7.5				# nm below cap height
    self.timeToMiss = 10			# us
    self.vestibuleRadius = 3.			# nm
    self.vestibuleHeight = 3.			# nm tall
    self.vestibulePos = 2.			# nm below cap height
    self.initHeight = 2.                        # nm above cap height
    self.offset = 0.                            # nm from center
    self.potExtent = 10.		# nm, radial from cap height center
    self.potA = 1.				# mV at cap height
    self.potB = 5.				# mV at bottom of vestibule
    self.potC = 100.				# mV at hitpoint
    self.mobility = 15e9			# nm^2/Vs
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = None
  Unity = False
  def copy(self) :
    return copy.deepcopy(self)
  def getCore(cls) :
    return smds.Cores.paramsToCore[cls]
  getCore = classmethod(getCore)
  def initCore(self) :
    return self.getCore()(self)
  def getNumBins(self) :
    return int(self.num)
  def toXML(self, doc) :
    "Returns these Params as a node for insertion into xml.dom.Document doc."
    n = doc.createElement("Params")
    n.appendChild(self.createXMLchannel(doc))
    n.appendChild(self.createXMLsim(doc))
    ## add extra values with smds.createXMLvalue() here ##
    return n
  def fromXML(cls, n) :
    """Returns a new smds.Params created from the xml.dom node n."""
    if (n == None or n.tagName != 'Params') : return None
    r = cls()
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE :
        if c.tagName == 'Channel' : r.interpretXMLchannel(c)
        if c.tagName == 'Sim' : r.interpretXMLsim(c)
    smds.interpretXMLvalues(r, n)
    return r
  fromXML = classmethod(fromXML)
  def createXMLchannel(self, doc) :
    n = doc.createElement("Channel")
    n.appendChild(smds.createXMLvalue(doc, "capHeight", "%g"%self.capHeight))
    n.appendChild(smds.createXMLvalue(doc, "capWidth", "%g"%self.capWidth))
    n.appendChild(smds.createXMLvalue(doc, "channelRadius",
                                      "%g"%self.channelRadius))
    n.appendChild(smds.createXMLvalue(doc, "hitPoint",
                                      "%g"%self.hitPoint))
    n.appendChild(smds.createXMLvalue(doc, "vestibuleHeight",
                                      "%g" % self.vestibuleHeight))
    n.appendChild(smds.createXMLvalue(doc, "vestibulePos",
                                      "%g" % self.vestibulePos))
    n.appendChild(smds.createXMLvalue(doc, "vestibuleRadius",
                                      "%g" % self.vestibuleRadius))
    n.appendChild(smds.createXMLvalue(doc, "initHeight",
                                      "%g" % self.initHeight))
    n.appendChild(smds.createXMLvalue(doc, "offset",
                                      "%g" % self.offset))
    n.appendChild(smds.createXMLvalue(doc, "potExtent",
                                      "%g" % self.potExtent))
    n.appendChild(smds.createXMLvalue(doc, "potA",
                                      "%g" % self.potA))
    n.appendChild(smds.createXMLvalue(doc, "potB",
                                      "%g" % self.potB))
    n.appendChild(smds.createXMLvalue(doc, "potC",
                                      "%g" % self.potC))
    return n
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.setAttribute("dT", str(self.dT))
    n.setAttribute("binWidth", str(self.binWidth))
    n.setAttribute("num", str(self.num))
    n.appendChild(smds.createXMLvalue(doc, "timeToMiss", "%g"%self.timeToMiss))
    n.appendChild(smds.createXMLvalue(doc, "D", "%g" % self.D))
    n.appendChild(smds.createXMLvalue(doc, "mobility", "%g" % self.mobility))
    return n
  def interpretXMLchannel(self, n) :
    smds.interpretXMLvalues(self, n)
  def interpretXMLsim(self, n) :
    self.dT = float(n.getAttribute("dT"))
    self.binWidth = float(n.getAttribute("binWidth"))
    self.num = int(n.getAttribute("num"))
    smds.interpretXMLvalues(self, n)
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'timeToMiss' :
      self.timeToMiss = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'D' :
      self.D = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'mobility' :
      self.mobility = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'capHeight' :
      self.capHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'capWidth' :
      self.capWidth = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'channelRadius' :
      self.channelRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'hitPoint' :
      self.hitPoint = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'vestibuleRadius' :
      self.vestibuleRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'vestibuleHeight' :
      self.vestibuleHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'vestibulePos' :
      self.vestibulePos = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'initHeight' :
      self.initHeight = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'offset' :
      self.offset = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potExtent' :
      self.potExtent = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potA' :
      self.potA = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potB' :
      self.potB = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'potC' :
      self.potC = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False
  def okay(self) :
    return True
  def strip(self) :
    return self

class Base_AHL_1Tape :
  """Simulation paramters with tape record of particle coordinates"""
  def __init__(self, p = None):
    self.stepsPerTape = 1                             # ns
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'stepsPerTape' :
      self.stepsPerTape = int(smds.strFromTxtNodes(c.childNodes))
      return True
    return False

class Base_AHL_2Tape :
  """Simulation paramters with hit and miss tape records of coordinates"""
  def __init__(self, p = None):
    self.stepsPerTape = 1           # ns until record
    self.tapeSelect = 1             # 0 for just hit tape, 1 for both, 2 for just miss tape
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    n.appendChild(smds.createXMLvalue(doc, "tapeSelect", "%g" % self.tapeSelect))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'stepsPerTape' :
      self.stepsPerTape = int(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'tapeSelect' :
      self.tapeSelect = int(smds.strFromTxtNodes(c.childNodes))
      return True
    return False

class Base_AHL_Spherical :
  """Simulation parameters with spherical particle"""
  def __init__(self, p = None) :
    self.particleRadius = 0.4       # nm
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLchannel(self, doc) :
    n = doc.createElement("Channel")
    n.appendChild(smds.createXMLvalue(doc, "particleRadius",
                                          "%g" % self.particleRadius))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'particleRadius' :
      self.particleRadius = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False

class Base_AHL_Electroosmotic :
  """Simulation paramerters with electroosmotic flow"""
  def __init__(self, p = None):
    self.conductance = 295          # pS
    self.ionicSelectivity = 0.55    # unitless K/Cl
    self.Nw = 10                    # water molecules per ion
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = doc.createElement("Sim")
    n.appendChild(smds.createXMLvalue(doc, "conductance", "%g" % self.conductance))
    n.appendChild(smds.createXMLvalue(doc, "ionicSelectivity", "%g" % self.ionicSelectivity))
    n.appendChild(smds.createXMLvalue(doc, "Nw", "%g" % self.Nw))
    return n
  def interpretXMLvalue(self, c) :
    name = c.getAttribute('Name')
    if name == 'conductance' :
      self.conductance = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'ionicSelectivity' :
      self.ionicSelectivity = float(smds.strFromTxtNodes(c.childNodes))
      return True
    if name == 'Nw' :
      self.Nw = float(smds.strFromTxtNodes(c.childNodes))
      return True
    return False

class cAHL_sdx (Base_AHL_sdx) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distribution"""
  pass

class cAHL_sdxp (Base_pot, Base_AHL_sdx) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Custom Potential Field, Static Launch Point, Time to Hit Distribution"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    Base_pot.__init__(self, p)
    Base_AHL_sdx.__init__(self,p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = None
  Unity = False
  def createXMLchannel(self, doc) :
    n = Base_AHL_sdx.createXMLchannel(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "PotProfile", self.pot.name))
    n.appendChild(smds.createXMLvalue(doc, "potA",
                                      "%g" % self.potA))
    n.appendChild(smds.createXMLvalue(doc, "potB",
                                      "%g" % self.potB))
    return n
  def interpretXMLvalue(self, c) :
    if Base_pot.interpretXMLvalue(self, c) : return True
    if Base_AHL_sdx.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sphdx (Base_AHL_sdx, Base_AHL_Spherical) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distribution"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    Base_AHL_sdx.__init__(self, p)
    Base_AHL_Spherical.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = None
  Unity = False
  def createXMLchannel(self, doc) :
    n = Base_AHL_sdx.createXMLchannel(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "particleRadius",
                                      "%g" % self.particleRadius))
    return n
  def interpretXMLvalue(self, c) :
    if Base_AHL_sdx.interpretXMLvalue(self, c) : return True
    if Base_AHL_Spherical.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sdxt (Base_AHL_sdx, Base_AHL_1Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distribution, Tape Record"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    Base_AHL_sdx.__init__(self, p)
    Base_AHL_1Tape.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  SpeciesType = None
  Unity = False
  def createXMLsim(self, doc) :
    n = Base_AHL_sdx.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    return n
  def interpretXMLvalue(self, c) :
    if(Base_AHL_sdx.interpretXMLvalue(self, c)) : return True
    if(Base_AHL_1Tape.interpretXMLvalue(self,c)) : return True
    return False

class cAHL_sdx2t (Base_AHL_sdx, Base_AHL_2Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distribution, Hit & Miss Tape Records"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    Base_AHL_sdx.__init__(self, p)
    Base_AHL_2Tape.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = Base_AHL_sdx.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    n.appendChild(smds.createXMLvalue(doc, "tapeSelect", "%g" % self.tapeSelect))
    return n
  def interpretXMLvalue(self, c) :
    if Base_AHL_sdx.interpretXMLvalue(self, c) : return True
    if Base_AHL_2Tape.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sphdx2t (Base_AHL_sdx, Base_AHL_Spherical, Base_AHL_2Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Static Launch Point, Time to Hit Distribution, Hit & Miss Tape Records, Spherical Particle"""
  def __init__(self, p = None) :
    """Constructor takes an optional dictionary of parameter values."""
    Base_AHL_sdx.__init__(self, p)
    Base_AHL_Spherical.__init__(self, p)
    Base_AHL_2Tape.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLchannel(self, doc) :
    n = Base_AHL_sdx.createXMLchannel(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "particleRadius",
                                      "%g"%self.channelRadius))
    return n
  def createXMLsim(self, doc) :
    n = Base_AHL_sdx.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    n.appendChild(smds.createXMLvalue(doc, "tapeSelect", "%g" % self.tapeSelect))
    return n
  def interpretXMLvalue(self, c) :
    if Base_AHL_sdx.interpretXMLvalue(self, c) : return True
    if Base_AHL_2Tape.interpretXMLvalue(self, c) : return True
    if Base_AHL_Spherical.interpretXMLvalue(self, c) : return True
    return False
  
class cAHL_sdxeo (cAHL_sdx, Base_AHL_Electroosmotic) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution"""
  def __init__(self, p = None) :
    cAHL_sdx.__init__(self, p)
    Base_AHL_Electroosmotic.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = Base_AHL_sdx.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "conductance", "%g" % self.conductance))
    n.appendChild(smds.createXMLvalue(doc, "ionicSelectivity", "%g" % self.ionicSelectivity))
    n.appendChild(smds.createXMLvalue(doc, "Nw", "%g" % self.Nw))
    return n
  def interpretXMLvalue(self, c) :
    if Base_AHL_sdx.interpretXMLvalue(self, c) : return True
    if Base_AHL_Electroosmotic.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sdxeot (cAHL_sdxeo, Base_AHL_1Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Tape Record"""
  def __init__(self, p = None) :
    cAHL_sdxeo.__init__(self, p)
    Base_AHL_1Tape.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = cAHL_sdxeo.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    return n
  def interpretXMLvalue(self, c) :
    if cAHL_sdxeo.interpretXMLvalue(self, c) : return True
    if Base_AHL_1Tape.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sdxeo2t (cAHL_sdxeo, Base_AHL_2Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Hit & Miss Tape Records"""
  def __init__(self, p = None) :
    cAHL_sdxeo.__init__(self, p)
    Base_AHL_2Tape.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = cAHL_sdxeo.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    n.appendChild(smds.createXMLvalue(doc, "tapeSelect", "%g" % self.tapeSelect))
    return n
  def interpretXMLvalue(self, c) :
    if cAHL_sdxeo.interpretXMLvalue(self, c) : return True
    if Base_AHL_2Tape.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sphdxeo (cAHL_sdxeo, Base_AHL_Spherical) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Spherical Particle"""
  def __init__(self, p = None):
    cAHL_sdxeo.__init__(self, p)
    Base_AHL_Spherical.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = cAHL_sdxeo.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "particleRadius", "%g" % self.particleRadius))
    return n
  def interpretXMLvalue(self, c):
    if cAHL_sdxeo.interpretXMLvalue(self, c) : return True
    if Base_AHL_Spherical.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sphdxeot (cAHL_sphdxeo, Base_AHL_1Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Spherical Particle, Tape Record"""
  def __init__(self, p = None):
    cAHL_sphdxeo.__init__(self, p)
    Base_AHL_1Tape.__init__(self, p)
    if p != None : self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc) :
    n = cAHL_sphdxeo.createXMLsim(self, doc)
    n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
    return n
  def interpretXMLvalue(self, c):
    if cAHL_sphdxeo.interpretXMLvalue(self, c) : return True
    if Base_AHL_1Tape.interpretXMLvalue(self, c) : return True
    return False

class cAHL_sphdxeo2t (cAHL_sphdxeo, Base_AHL_2Tape) :
  """Simulation parameters, Alpha-Hemolysin Core, Electrophoretic Flow, Exponential Field, Electroosmotic Flow, Static Launch Point, Time to Hit Distribution, Spherical Particle, Tape Record"""
  def __init__(self, p = None):
      cAHL_sphdxeo.__init__(self,p)
      Base_AHL_2Tape.__init__(self,p)
      if p != None :self.__dict__.update(copy.deepcopy(p))
  def createXMLsim(self, doc):
      n = cAHL_sphdxeo.createXMLsim(self, doc)
      n.appendChild(smds.createXMLvalue(doc, "stepsPerTape", "%g" % self.stepsPerTape))
      return n
  def interpretXMLvalue(self, c):
      if cAHL_sphdxeo.interpretXMLvalue(self, c) : return True
      if Base_AHL_1Tape.interpretXMLvalue(self, c) : return True
      return False