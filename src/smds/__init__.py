# smds module
#  Single Molecule Diffusion Simulator
#  $Id: __init__.py,v 1.24 2009/04/16 21:58:04 mculbert Exp $

# To prevent attempting to load c:\windows\system32\zlib.dll on windows
def fixPath() :
  from sys import path
  from os import getcwd, chdir
  for x in path :
    if x.lower() == 'c:\\windows\\system32' or x.lower() == 'c:\\winnt\\system32' :
      path.remove(x)
  x = getcwd().lower()
  if x == 'c:\\windows\\system32' or x == 'c:\\winnt\\system32' :
    chdir('c:\\')
fixPath()

import xml.dom
from xml.dom.minidom import parse
from time import strftime
from sys import stderr
import Analysis, Params, Results, Task
import Cores
import Engine, Dispatch

verbosity = 2		# 0 = silent, 1 = Warnings, 2 = Info, 3 = Debug
MSG_WARN = 1
MSG_INFO = 2
MSG_DBUG = 3

addToQueue = queueIsEmpty = clearQueue = None
initialized = False

from Correlator import Correlator
Analysis.Correlator = Correlator
del Correlator

from Tabulator import Tabulator
Analysis.Tabulator = Tabulator
del Tabulator
#from Tabulator_PCH import Tabulator_PCH
#Analysis.Tabulator_PCH = Tabulator_PCH
#del Tabulator_PCH
from Tabulator_IPCH import Tabulator_IPCH
Analysis.Tabulator_IPCH = Tabulator_IPCH
del Tabulator_IPCH

def notInit(*args) :
  raise "SMDS is not initialized."

addToQueue = notInit
queueIsEmpty = notInit
clearQueue = notInit
theCache = None

def shutdown() :
  global addToQueue, queueIsEmpty, clearQueue
  msg("SMDS Shutting down.")
  Engine.shutdown()
  Dispatch.shutdown()
  if theCache : theCache.close()
  addToQueue = notInit
  queueIsEmpty = notInit
  clearQueue = notInit

def initialize(dispatch = True, seed = None, cache = None, autoClean = True) :
  """initialize(dispatch = True, seed = None, cache = None, autoClean = True)
Initialize the Single-Molecule Diffusion Simulator.
  If dispatch is True, SMDS will attempt to contact the local Foundry to
  operate in distributed mode.  If this contact fails or if dispatch is
  False, SMDS will operate in single-processor mode.  An optional seed for
  the random number generator may be provided, otherwise the generator is
  initialized with the current time.\n"""
  global initialized, addToQueue, queueIsEmpty, clearQueue, theCache
  if initialized : return
  if dispatch and Dispatch.initialize() :
    addToQueue = Dispatch.addToQueue
    queueIsEmpty = Dispatch.queueIsEmpty
    clearQueue = Dispatch.clearQueue
    msg("SMDS Initialized in Dispatch mode.")
  else :
    Engine.initialize()
    addToQueue = Engine.addToQueue
    queueIsEmpty = Engine.queueIsEmpty
    clearQueue = Engine.clearQueue
    msg("SMDS Initialized in Single-processor mode.")
  if seed : Cores.c2D.sRnd(int(seed))
  else    : Cores.c2D.sRnd()
  initialized = True
  if cache : theCache = Cache(cache, autoClean)
  # Ensure proper shutdown()
  from atexit import register
  from sys import exit
  from signal import signal, SIGTERM
  register(shutdown)
  signal(SIGTERM, lambda x,y : exit())
  try :
    from signal import SIGBREAK
    signal(SIGBREAK, lambda x,y : exit())
  except :
    pass

def xmlRead(filename = None) :
  """xmlRead(<filename>) -> [smds.Batch]
  Reads the specified SMDS XML file or stdin if none is supplied and returns
  a list of batches of tasks."""
  r = []
  acc = Batch()
  if filename == None :
    import sys
    doc = parse(sys.stdin)
  else :
    doc = parse(filename)
  if doc.firstChild.tagName != 'SMDS' :
    raise TypeError, "Root node of " + (filename or '<stdin>') + \
    			" is not <SMDS>,"
  for c in doc.firstChild.childNodes :
    if c.nodeName == 'Task' : acc.addTask(Task.Task.fromXML(c))
    if c.nodeName == 'Batch' :
      if acc.numTasks() > 0 :
        r.append(acc)
        acc = Batch()
      b = Batch.fromXML(c)
      if b : r.append(b)
  if acc.numTasks() > 0 : r.append(acc)
  doc.unlink()
  return r

def xmlWrite(data, filename = None) :
  """xmlWrite([smds.Task.Task | smds.Batch], <filename>)
  Writes an SMDS XML file composed of the supplied tasks and batches of
  tasks to the specified file or stdout."""
  doc = xml.dom.getDOMImplementation().createDocument(None, 'SMDS', None)
  if type(data) != type(()) and type(data) != type([]) : data = (data,)
  for x in data :
    doc.firstChild.appendChild(x.toXML(doc))
  if filename == None :
    import sys
    out = sys.stdout
  else :
    out = open(filename, 'w')
  out.write(doc.toprettyxml("  ") + "\n")
  out.close()
  doc.unlink()

def strFromTxtNodes (n, CDATA = False) :
  r = ''
  for c in n: 
    if (not CDATA and c.nodeType == xml.dom.Node.TEXT_NODE) or		\
       (    CDATA and c.nodeType == xml.dom.Node.CDATA_SECTION_NODE) :
      r += c.data
  return r

def createXMLvalue (doc, name, v, CDATA = False) :
  c = doc.createElement("Value")
  c.setAttribute("Name", name)
  if CDATA : c.appendChild(doc.createCDATASection(v))
  else     : c.appendChild(doc.createTextNode(v))
  return c

def interpretXMLvalues (o, n) :
  for c in n.childNodes :
    if c.nodeType == xml.dom.Node.ELEMENT_NODE and    \
       c.tagName == 'Value' and not o.interpretXMLvalue(c) :
      from warnings import warn
      warn("Unknown XML Value %s!" % c.getAttribute('Name'),      \
                      stacklevel = 2)

def msg(s, level = MSG_INFO) :
  if verbosity >= level :
    stderr.write("[%s] %s\n" % (strftime("%H:%M.%S"),s))
    stderr.flush()

def defaultCallback (task) :
  if task.status == Task.Status_Error :
    msg("Task <%s> had an error." % task.ID, MSG_WARN)
  else :
    msg("Task <%s> completed." % task.ID)

class socketIO :
  def __init__(self, socket) :
    self.socket = socket
  def read(self, n) :
    r = ''
    while n > 0 :
      buf = self.socket.recv(n)
      if not buf : return r
      r += buf
      n -= len(buf)
    return r
  def readline(self) :
    r = ''
    while True :
      c = self.read(1)
      r += c
      if not c or c == '\n' : return r
  def write(self, str) :
#    self.socket.sendall(str)
    pos = 0
    rem = len(str)
    while rem > 0 :
      sent = self.socket.send(str[pos:pos+1024])
      pos += sent
      rem -= sent

# Takes an open socket, an integer message, and parameters
# Returns True on successful send
def sendMsg(socket, msg, *args) :
  from cPickle import dump
  try    : dump((msg,args), socketIO(socket), 2)
  except : return False
  else   : return True

# Takes on open socket with data to be read
# Returns None on error or (msg, (args)) on success.
def recvMsg(socket) :
  from cPickle import load
  try    : packet = load(socketIO(socket))
  except : return None
  return packet

def junk(socket):
  (r, buf) = ([], '')
  while True :
    if len(buf) < 8 : 			# get the remainder of the header
      str = buf
      size = 8-len(str)
      try: buf = socket.recv(4096)
      except : return None
      if len(buf) < size : return None
      str += buf[:size]
      buf = buf[size:]
    else :
      str = buf[:8]
      buf = buf[8:]
    (msg, size) = unpack("LL", str)	# unpack header
    str = ''
    while True :
      if len(buf) == size :		# just enough, we're done.
        r.append((msg, loads(str+buf)))
        return r
      if len(buf) < size :		# not enough data yet
        str += buf
        size -= len(buf)
      else :				# too much data
        r.append((msg, loads(str+buf[:size])))
        buf = buf[size:]
        break
        
messages = {
		'QUIT'   : 1,
		'START'  : 2,
		'WORK'   : 3,
		'INTENS' : 4,
		'ADD'    : 5,
		'REMOVE' : 6,
		'FOUNDRY': 7,
		'ETOC'   : 8,
		'RELEASE': 9,
		'CALC'   : 10,
	   }

def addHost(h) :
  from smds.Foundry import FOUNDRY_PORT
  try :
    from socket import socket
    s = socket()
    s.connect(('localhost', FOUNDRY_PORT))
    sendMsg(s, messages['ADD'], h)
  except :
    pass

def removeHost(h) :
  from smds.Foundry import FOUNDRY_PORT
  try :
    from socket import socket
    s = socket()
    s.connect(('localhost', FOUNDRY_PORT))
    sendMsg(s, messages['REMOVE'], h)
  except :
    pass

class Batch :
  """A collection of tasks."""
  def __init__(self, data = None, ID = None) :
    self.ID = ID or ""
    self.data = data or []
    self.idMap = {}
    for i in range(len(self.data)) :
      self.idMap[data[i].ID] = i
  def __getitem__ (self, ID) :
    """Get a task by ID."""
    return self.data[self.idMap[ID]]
  def toXML(self, doc) :
    """Returns this Batch as a node for insertion into xml.dom.Document doc."""
    n = doc.createElement("Batch")
    if self.ID != "" :
      n.setAttribute("ID", self.ID)
    for t in self.data :
      n.appendChild(t.toXML(doc))
    return n
  def fromXML(n) :
    """Returns a new smds.Batch created from the xml.dom node n."""
    if (n == None or n.tagName != 'Batch') : return None
    r = Batch()
    r.ID = n.getAttribute('ID')
    for c in n.childNodes :
      if c.nodeType == xml.dom.Node.ELEMENT_NODE :
        if c.tagName == 'Task' : r.addTask(Task.Task.fromXML(c))
    return r
  fromXML = staticmethod(fromXML)
  def addTask(self, task) :
    """Add the given task to the batch."""
    if task == None : return
    self.idMap[task.ID] = len(self.data)
    self.data.append(task)
  def addToQueue(self, callback = defaultCallback) :
    """addToQueue(<callback>)
    Add the tasks in this batch to the process queue."""
    for t in self.data : addToQueue(t, callback)
  def numTasks(self) :
    return len(self.data)

class Cache :
  def __init__(self, filename, autoClean = True) :
    import anydbm
    self.db = anydbm.open(filename, 'c')
    self.file = filename
    self.autoClean = autoClean
    self.doc = xml.dom.getDOMImplementation().createDocument(None,
                          'SMDS', None)
    msg("Opened cache: %s" % filename, MSG_DBUG)
  def fetch(self, task) :
    if not self.db.has_key(task.ID) : return False
    from cPickle import loads
    task.__dict__.update(Task.Task.fromXML(loads(self.db[task.ID])).__dict__)
    msg("Fetched <%s> from the cache." % task.ID, MSG_DBUG)
    return True
  def put(self, task) :
    from cPickle import dumps
    self.db[task.ID] = dumps(task.toXML(self.doc), 2)
    msg("Stored <%s> in the cache." % task.ID, MSG_DBUG)
  def close(self) :
    self.db.close()
    if self.autoClean :
      msg("Autocleaning cache: %s" % self.file, MSG_DBUG)
      from os import remove
      for ext in ['', '.db', '.dir', '.dat', 'pag'] :
        try    : remove(self.file+ext)
        except : pass

def numCPUs() :
  from sys import platform
  try :
    if platform == 'linux2' :
      from os import sysconf
      return sysconf('SC_NPROCESSORS_ONLN')
    if platform == 'win32' :
      from os import environ
      return int(environ['NUMBER_OF_PROCESSORS'])
    if platform == 'darwin' :
      from os import popen
      return int(popen('/usr/sbin/sysctl -n hw.ncpu').read())
    if platform == 'netbsd3' :
      from os import popen
      return int(popen('/sbin/sysctl -n hw.ncpu').read())
  except :
    msg("Error identifying number of CPUs.  Defaulting to 1.", MSG_WARN)
    return 1
  msg("Unknown platform.  Defaulting to 1 CPU.", MSG_WARN)
  return 1
