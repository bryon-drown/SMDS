# smds.Calc module
#   Generalized distributed computation of a supplied function
# $Id: Calc.py,v 1.4 2008/02/28 16:29:02 mculbert Exp $

__all__ = [ 'registerFunction', 'registerData' ]
import smds, socket, numpy, math, numpy.oldnumeric.fft as FFT #FFTf
from sys import platform, exit
from select import select
from time import time
from types import FunctionType
from cPickle import dumps as pDump, loads as pLoad
from marshal import dumps as mDump, loads as mLoad
from zlib import compress, decompress

if platform != 'win32' :
  from resource import getrlimit, setrlimit, RLIMIT_DATA
  limit = getrlimit(RLIMIT_DATA)
  setrlimit(RLIMIT_DATA, (limit[1], limit[1]))

################################
###  Public interface
################################

def registerFunction(func) :
  """registerFunction(func) -> distributedFunc
  Returns the function wrapped for distributed computation.
  The numpy and math modules are automatically imported, and data stored
  with registerData() can be retrieved by calling fetchData(name).

  The distributed function that is returned may be called with a list of
  tuples that provide the sets of parameters with which to evaluate the
  function, and it returns a list of the return values for each tuple."""
  global funcID
  if not initialized : initialize()
  id = '_func'+str(funcID)
  funcID += 1
  smds.msg("Registering function <%s>" % id, smds.MSG_DBUG)
  store[id] = mDump(func.func_code, 1)
  return lambda args : execute(id, args)

def registerData(name, data) :
  """registerData(name, data)
  Stores data for retrieval by remote processors.  Names must be unique.
  In registered functions, data can be retrived by calling fetchData(name)."""
  if name in store :
    raise ValueError, "An object named <%s> already exists." % name
  smds.msg("Registering data <%s>" % name, smds.MSG_DBUG)
  store[name] = compress(pDump(data, 2))

################################
###  Private variables
################################

initialized = False
store = {}
funcID = 0
workID = 0
theSocket = None
theFoundry = None
myPort = None
hosts = []
waiting = {}
bees = []
queue = []
working = 0
done = None
theFunc = None

################################
###  Remote routines
################################

def retrieveData(name) :
  smds.msg("Requesting <%s> from parent." % name, smds.MSG_DBUG)
  try :
    smds.sendMsg(theSocket, smds.messages['INTENS'], name)
    packet = smds.recvMsg(theSocket)
    if not packet :
      theSocket.close()
      exit("Error reading from parent.")
    (msg, args) = packet
    if msg == smds.messages['QUIT'] : quit()
    if msg != smds.messages['INTENS'] :
      smds.msg("Unknown message from parent: %s.\nGiving up." 
      			% str(msg, args), smds.MSG_WARN)
      quit()
    if args[0] != name :
      smds.msg("Unknown data <%s> from parent.  Giving up.", smds.MSG_WARN)
      quit()
  except :
    from traceback import format_exc
    smds.msg("Error retrieving data <%s> from master." % name, smds.MSG_WARN)
    print format_exc()
    exit()
  return args[1]

def fetchData(name) :
  if name in store : return store[name]
  data = pLoad(decompress(retrieveData(name)))
  store[name] = data
  return data

def moreWork(workID, funcID, args) :
  smds.msg("moreWork(%s, %s)" % (funcID, str(args)), smds.MSG_DBUG)
  if funcID not in store :
    store[funcID] = FunctionType(mLoad(retrieveData(funcID)), env)
  func = store[funcID]
  try : ret = func(*args) ; err = False
  except SystemExit : raise
  except :
    from traceback import format_exc
    ret = format_exc()
    err = True
  smds.msg("Returning %s to the master." % str(ret), smds.MSG_DBUG)
  if err : r = smds.sendMsg(theSocket, smds.messages['RELEASE'], workID, ret)
  else   : r = smds.sendMsg(theSocket, smds.messages['WORK'], workID, ret)
  if not r :
    smds.msg("Error sending work to master.", smds.MSG_WARN)
    exit()

def quit():
  smds.msg("Quit!")
  exit()

env = { 'fetchData' : fetchData, 'abs' : abs, 'apply' : apply, 'chr' : chr,
        'cmp' : cmp, 'complex' : complex, 'delattr' : delattr, 'dict' : dict,
        'True' : True, 'False' : False, 'divmod' : divmod,
        'enumerate' : enumerate, 'filter' : filter, 'float' : float,
        'hasattr' : hasattr, 'hash' : hash, 'hex' : hex, 'int' : int,
        'isinstance' : isinstance, 'issubclass' : issubclass, 'iter' : iter,
        'len' : len, 'list' : list, 'long' : long, 'map' : map, 'max' : max,
        'min' : min, 'oct' : oct, 'ord' : ord, 'pow' : pow, 'range' : range,
        'reduce' : reduce, 'round' : round, 'slice' : slice, 'str' : str,
        'sum' : sum, 'tuple' : tuple, 'type' : type, 'vars' : vars,
        'xrange' : xrange, 'zip' : zip, 'ValueError' : ValueError,
        'TypeError' : TypeError, 'ArithmeticError' : ArithmeticError,
        'AssertionError' : AssertionError, 'AttributeError' : AttributeError,
        'FloatingPointError' : FloatingPointError, 'IndexError' : IndexError,
        'KeyError' : KeyError, 'NotImplementedError' : NotImplementedError,
        'OverflowError' : OverflowError, 'RuntimeError' : RuntimeError, }
for m in (FFT, math, numpy) :
  for k in dir(m) :
    if k[0] != '_' : env[k] = m.__dict__[k]
for k in [ 'math', 'copy', 'oldnumeric', 'rec', 'lib', 'core' ] :
  del env[k]

#### Remote Initialization

def nullHandler(*args) :
  smds.msg("Unknown message: " + str(args), smds.MSG_WARN)

remoteHandlers = [
		nullHandler,		# 0
		quit,			# 1
		nullHandler,		# START, 2
		moreWork,		# WORK, 3
	   ]
maxRemoteHandler = 3

def run(host, port) :
  global theSocket
  if platform == 'win32' :
    from win32process import GetCurrentProcess, SetPriorityClass, \
  				IDLE_PRIORITY_CLASS
    SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS)
  else :
    from os import nice
    nice(15)

  from atexit import register
  from signal import signal, SIGTERM
  from sys import exit
  register(shutdown)
  signal(SIGTERM, lambda x,y : exit() )
  if platform == 'win32' :
    from signal import SIGBREAK
    signal(SIGBREAK, lambda x,y : exit() )

  # Connect to parent
  smds.verbosity = smds.MSG_DBUG
  theSocket = socket.socket()
  theSocket.connect((host,port))

  while (True) :
    (r,w,e) = select([theSocket],[],[],0)
    if len(r) :
      packet = smds.recvMsg(theSocket)
      if not packet :
        theSocket.close()
        exit("Error reading from parent.")
      (msg, args) = packet
      if (msg > maxRemoteHandler) :
        h = nullHandler
        args = (msg,)+args
      else : h = remoteHandlers[msg]
      h(*args)

if __name__ == "__main__" :
  from sys import argv
  if len(argv) != 3 :
    exit("Usage: Calc.py host port")
  run(argv[1], int(argv[2]))

####################################
###  Local routines
####################################

def initialize() :
  global theFoundry, theSocket, myPort, funcID, workID, initialized
  from os import getpid
  from errno import EADDRINUSE
  from atexit import register
  from signal import signal, SIGTERM
  from sys import exit
  from random import randint
  register(shutdown)
  signal(SIGTERM, lambda x,y : exit() )
  try :
    from signal import SIGBREAK
    signal(SIGBREAK, lambda x,y : exit())
  except :
    pass
  funcID = randint(0, 1000000)
  workID = randint(0, 1000)
  initialized = True
  if theFoundry == None :
    # Attempt to connect to local foundary
    s = socket.socket()
    s.connect(('localhost', smds.Foundry.FOUNDRY_PORT))
    if not smds.sendMsg(s, smds.messages['WORK'], getpid()) :
      raise IOError, "Could not register with local Foundry."
    theFoundry = s

    p = 8740
    while myPort == None :
      try :
        s = socket.socket()
        s.bind(('', p))
      except socket.error, (errno, errstr) :
        if errno == EADDRINUSE : p += 1
        else : raise
      else :
        myPort = p
    s.setblocking(0)
    s.listen(30)
    theSocket = s
    smds.msg("smds.Calc initialized.", smds.MSG_DBUG)

def shutdown() :
  global theFoundry, theSocket, myPort, bees, hosts, waiting, initialized
  if theFoundry : theFoundry.close()
  if theSocket : theSocket.close()
  (theFoundry, myPort, theSocket) = (None, None, None)
  for b in bees :
    smds.sendMsg(b['socket'], smds.messages['QUIT'])
    b['socket'].close()
  bees = []
  hosts = []
  waiting = {}
  initialized = False

def addHost(host) :
  for x in hosts + waiting.keys() + [ b['host'] for b in bees ] :
    if x == host :
      return
  hosts.append(host)
  smds.msg("Adding host %s" % host, smds.MSG_DBUG)

def removeHost(host) :
  smds.msg("Host %s down" % host, smds.MSG_DBUG)
  for b in bees :
    if b['host'] == host :
      smds.sendMsg(b['socket'], smds.messages['QUIT'])
      remoteDead(b)
  if host in waiting : del waiting[host]
  if host in hosts   : hosts.remove(host)

def workDone(b, workID, ret) :
  global working
  if not b['work'] or b['work'][2] != workID :
    smds.msg("Spurious message from " + b['host'], smds.MSG_DBUG)
    return
  done[b['work'][0]] = ret
  b['work'] = None
  working -= 1
  smds.msg("Work from %s (%d/%d avail)" %
		(b['host'], len(bees)-working, len(bees)), smds.MSG_DBUG)
  if working and working < 4 :
    smds.msg("Waiting on: " + ' '.join([ b['host']
                        for b in bees if b['work'] ]), smds.MSG_DBUG)

def workError(b, workID, ret) :
  global working
  if not b['work'] or b['work'][2] != workID :
    smds.msg("Spurious message from " + b['host'], smds.MSG_DBUG)
    return
  smds.msg("Error from %s on argument %d:\n%s" %
  		(b['host'], b['work'][0], ret), smds.MSG_WARN)
  b['work'] = None
  working -= 1

def dataRequest(b, name) :
  if name not in store :
    smds.msg("Unnown data request <%s> by %s." % (name, b['host']))
    remoteDead(b)
  else :    
    smds.msg("Data [%s] request by %s" % (name, b['host']), smds.MSG_DBUG)
    if not smds.sendMsg(b['socket'], smds.messages['INTENS'],
    				name, store[name]) :
      remoteDead(b)

def remoteDead(b) :
  global working
  smds.msg("Remote %s dead." % b['host'], smds.MSG_DBUG)
  if b['work'] : 
    queue.append(b['work'][:-1])
    working -= 1
  bees.remove(b)
  hosts.append(b['host'])

handlers = [
                nullHandler,            # 0
                nullHandler,            # QUIT, 1
                nullHandler,            # START, 2
                workDone,               # WORK, 3
                dataRequest,            # INTENS, 4
                addHost,                # ADD, 5
                removeHost,             # REMOVE, 6
                nullHandler,            # FOUNDRY, 7
                nullHandler,            # ETOC, 8
                workError,              # RELEASE, 9
           ]
maxHandler = 9

def processMessages() :
  (r,w,e) = select([theFoundry]+[ b['socket'] for b in bees ], [],[], 1.0)
  for s in r :
    packet = smds.recvMsg(s)
    b = None
    if s != theFoundry : b = [ b for b in bees if b['socket'] == s ][0]
    if not packet :
      if s == theFoundry :
        raise "Local foundry went down!"
      else :
        remoteDead(b)
      continue
    (msg, args) = packet
    if msg > maxHandler :
      h = nullHandler
      args = (msg,)+args
    else :
      h = handlers[msg]
      if b : args = (b,)+args
    h(*args)

def sendWork() :
  global working, workID
  avail = [ b for b in bees if b['work'] == None ]
  while (len(queue) > 0 and len(avail) > 0) :
    b = avail.pop()
    (i, args) = queue.pop()
    working += 1
    workID += 1
    b['work'] = (i, args, workID)
    if not smds.sendMsg(b['socket'], smds.messages['WORK'],
			workID, theFunc, args) :
      remoteDead(b)
    else :
      smds.msg("Sent %s%s <%d> to %s (%d remaining in queue)" %
			(theFunc, str(args), workID, b['host'], len(queue)),
			smds.MSG_DBUG)

def touchHosts() :
  global waiting
  for host in waiting.keys() :
    if time()-waiting[host] > 60.0 :
      del waiting[host]
      hosts.append(host)
  for host in [ x for x in hosts ] :
    s = socket.socket()
    s.settimeout(2.0)
    try : s.connect((host, smds.Foundry.FOUNDRY_PORT))
    except : pass
    else :
      if smds.sendMsg(s, smds.messages['CALC'],
                                socket.gethostname(), myPort) :
        waiting[host] = time()
        hosts.remove(host)
      s.close()

def incomingConnections() :
  go = True
  while go :
    try : (c, addr) = theSocket.accept()
    except : go = False
    else :
      r = True
      c.setblocking(1)
      (host, aliases, addrs) = socket.gethostbyaddr(addr[0])
      if host in waiting : del waiting[host]
      if host in hosts   : hosts.remove(host)
      if host in [ b['host'] for b in bees ] :
        smds.sendMsg(c, smds.messages['QUIT'])
      else :
        bees.append({ 'host' : host, 'socket' : c,
                      'work' : None, })
        smds.msg("New Bee: %s (%d/%d avail)" %
                        (host, len(bees)-working, len(bees)), 
			smds.MSG_DBUG)

def execute(funcID, args) :
  global queue, working, done, theFunc
  if not initialized or funcID not in store :
    raise ValueError, "Unknown function."
  smds.msg("Executing <%s> with %d parameter sets." % (funcID, len(args)),
				smds.MSG_DBUG)
  theFunc = funcID
  done = [None]*len(args)
  queue = [ (i,args[i]) for i in range(len(args)) ]
  while (len(queue) or working) :
    processMessages()
    touchHosts()
    incomingConnections()
    sendWork()
  ret = done
  done = None
  theFunc = None
  return ret
