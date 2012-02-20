# smds.Foundry module
#  Provides connection initialization between different nodes of the VM.
#  $Id: Foundry.py,v 1.18 2009/04/16 21:58:04 mculbert Exp $
 
import smds, socket, threading
from select import select
from os import environ
from sys import platform
from time import sleep

FOUNDRY_PORT = 8738
REFRESH_TIME = 5	# s
#REMOTE_CMD = 'python2.4 -c "from smds.Remote import run; run(%s, %d)"'
if environ.has_key('SMDS_HOSTS') :
  HOSTS_FILE = environ['SMDS_HOSTS']
else :
  try :
    HOSTS_FILE = environ['HOME'] + '/.smds_hosts'
  except :
    HOSTS_FILE = '.smds_hosts'

hosts = []		# { name, socket, usable }
myMasters = {}		# socket --> id
myName = ''
includeMe = False
numChildren = 0
theSocket = socket.socket()
theSocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
theThread = None
mx = threading.RLock()
work = threading.Event()

#### Spawn routines

def REMOTE_CMD(master, port, seed) :
  if platform == 'win32' :
    return 'python -c "from smds.Remote import run; run(\'%s\', %d, %d)"' \
    		% (master, port, seed)
  else :
    return ('python', 'python', '-c', 
  	  "from smds.Remote import run; run('%s', %d, %d)" % 
		(master, port, seed))

def CALC_REMOTE_CMD(master, port) :
  if platform == 'win32' :
    return 'python -c "from smds.Calc import run; run(\'%s\', %d)"' \
    		% (master, port)
  else :
    return ('python', 'python', '-c', 
  	  "from smds.Calc import run; run('%s', %d)" % 
		(master, port))

def spawnChild(master, port, seed) :
  if platform == 'win32' :
    from win32api import SearchPath
    from win32process import CreateProcess, STARTUPINFO
    (path, d) = SearchPath(None, 'python', '.exe')
    CreateProcess(path, REMOTE_CMD(master, port, seed), None, None, 0,
    			0, None, None, STARTUPINFO())
  else :
    from os import spawnlp, P_NOWAIT
    spawnlp(P_NOWAIT, *(REMOTE_CMD(master, port, seed)))

def spawnCalcChild(master, port) :
  if platform == 'win32' :
    from win32api import SearchPath
    from win32process import CreateProcess, STARTUPINFO
    (path, d) = SearchPath(None, 'python', '.exe')
    CreateProcess(path, CALC_REMOTE_CMD(master, port), None, None, 0,
    			0, None, None, STARTUPINFO())
  else :
    from os import spawnlp, P_NOWAIT
    spawnlp(P_NOWAIT, *(CALC_REMOTE_CMD(master, port)))

def cleanChildren() :
  global numChildren
  import os
  if numChildren > 0 :
    if os.waitpid(0, os.WNOHANG) != (0, 0) : numChildren -= 1
    if len(myMasters) == 0 and numChildren == 0 : theSocket.settimeout(None)

def refreshCodeImage() :
  if platform == 'win32' :
    from win32api import SearchPath
    from win32process import CreateProcess, STARTUPINFO
    (path, d) = SearchPath(None, 'python', '.exe')
    CreateProcess(path, 
		'python -c "from smds.starter import restart; restart()"', 
			None, None, 0,
    			0, None, None, STARTUPINFO())
  else :
    from os import spawnlp, P_NOWAIT
    spawnlp(P_NOWAIT, ('python', 'python', '-c',
          "from smds.starter import restart; restart()") )
  
#### Interface routines

def quit() :
  smds.msg("Quit!", smds.MSG_DBUG)
  from sys import exit
  exit()

def shutdown() :
  global theSocket
  if theSocket : theSocket.close()
  theSocket = None
  for h in hosts :
    if h['socket'] : h['socket'].close(); h['socket'] = None
  for s in myMasters.keys() :
    s.close()

def addHost(name) :
  global includeMe
  (host, aliases, addrs) = socket.gethostbyaddr(socket.gethostbyname(name))
  smds.msg("addHost(%s) as %s " % (name, host))
  if host == myName :
    includeMe = True
    return
  for h in hosts :
    if h['name'] == host :
      if ((not h['usable']) and h['socket'])  :
        h['usable'] = True
        activateFoundry(h)
      else : h['usable'] = True
      return
  hosts.append({ 
  		 'name' : host,
                 'socket' : None,
                 'usable' : True,
               })

def removeHost(name) :
  global includeMe
  (host, aliases, addrs) = socket.gethostbyaddr(socket.gethostbyname(name))
  smds.msg("removeHost(%s) as %s" % (name, host))
  if host == myName :
    includeMe = False
    deactivateFoundry({'name':myName})
    return
  for h in hosts :
    if h['name'] == host :
      mx.acquire()
      h['usable'] = False
      mx.release()
      deactivateFoundry(h)
      return

def registerFoundry(s, name) :
  (host, aliases, addrs) = socket.gethostbyaddr(socket.gethostbyname(name))
  smds.msg("registerFoundry(%s) as %s" % (name, host))
  x = [ h for h in hosts if h['name'] == host ]
  if len(x) :
    h = x[0]
    h['socket'] = s
  else :
    h = { 'name' : host, 'socket' : s, 'usable' : False }
    hosts.append(h)
  if h['usable'] : activateFoundry(h)

def registerMaster(s, id) :
  global includeMe
  id = str(id)
  mx.acquire()
  smds.msg("registerMaster(%s)" % id)
  refreshFoundries(True)
  if includeMe:
    if not smds.sendMsg(s, smds.messages['ADD'], myName) :
      masterDead(s, id)
      return
  for h in [ x for x in hosts if x['usable'] and x['socket'] ] :
    if not smds.sendMsg(s, smds.messages['ADD'], h['name']) :
      masterDead(s, id)
      return
  myMasters[s] = id
  theSocket.settimeout(300.0)
  mx.release()

def startRemote(master, port, seed) :
  global numChildren
  smds.msg("startRemote(%s, %d, %d)" % (master, port, seed))
  spawnChild(master, port, seed)
  numChildren += 1
  theSocket.settimeout(300.0)

def calc(master, port) :
  global numChildren
  smds.msg("calc(%s, %d)" % (master, port))
  spawnCalcChild(master, port)
  numChildren += 1
  theSocket.settimeout(300.0)

#### Utility routines

def activateFoundry(h) :
  mx.acquire()
  if not h['usable'] :
    mx.relase()
    return
  smds.msg("activateFoundry(%s)" % h['name'])
  defunct = []
  for m in myMasters :
    if not smds.sendMsg(m, smds.messages['ADD'], h['name']) :
      defunct.append(m)
  for m in defunct :
    masterDead(m)
  mx.release()

def deactivateFoundry(h) :
  smds.msg("deactivateFoundry(%s)" % h['name'])
  defunct = []
  for m in myMasters :
    if not smds.sendMsg(m, smds.messages['REMOVE'], h['name']) :
      defunct.append(m)
  for m in defunct :
    masterDead(m)

# masters and foundries don't send data after their initial registration,
# so, if select() says they're ready for reading, the connection has died.
def checkConnections() :
  if len(myMasters) :
    (r,w,e) = select(myMasters.keys(), [], [], 0)
    for m in r :
      masterDead(m)
  
  l = [ h['socket'] for h in hosts if h['socket'] ]
  if len(l) :
    (r,w,e) = select(l, [],[],0)
    if len(r) :
      for h in hosts :
        if h['socket'] and h['socket'] in r :
          h['socket'] = None
          deactivateFoundry(h)

def refreshFoundries(fast = False) :
  checkConnections()
  for h in [x for x in hosts if (x['usable'] and not x['socket']) ] :
    s = socket.socket()
    if fast : s.settimeout(1.0)
    else    : s.settimeout(10.0)
    try : s.connect((h['name'], FOUNDRY_PORT))
    except : pass
    else :
      if smds.sendMsg(s, smds.messages['FOUNDRY'], socket.gethostname()) :
        h['socket'] = s
        activateFoundry(h)

def masterDead(s, name = None) :
  if not name : name = myMasters[s]
  smds.msg("Master %s died." % name)
  if s in myMasters : del myMasters[s]
  if len(myMasters) == 0 and numChildren == 0 : theSocket.settimeout(None)

#### Initialization routines

def loadHosts() :
  global myName
  (myName, aliases, addrs) = socket.gethostbyaddr(
				socket.gethostbyname(socket.gethostname()))
  try :
    for hostname in open(HOSTS_FILE).readlines() :
      hostname = hostname.strip()
      if len(hostname) == 0 or hostname[0] == "#" : continue
      try : addHost(hostname)
      except : pass
  except :
    pass

def nullHandler(*args) :
  smds.msg("Unknown message: " + str(args), smds.MSG_WARN)

handlers = [
		nullHandler,		# 0
		quit,			# QUIT, 1
		startRemote,		# START, 2
		registerMaster,		# WORK, 3
		nullHandler,		# INTENS, 4
		addHost,		# ADD, 5
		removeHost,		# REMOVE, 6
		registerFoundry,	# FOUNDRY, 7
		nullHandler,		# 8
		refreshCodeImage,	# RELEASE, 9
		calc,			# CALC, 10
	   ]
maxHandler = 10

def bkg() :
  from time import sleep
  while True :
    work.wait()
    sleep(300)
    refreshFoundries()

def run() :
  global theSocket, theThread
  from atexit import register
  from sys import exit
  from signal import signal, SIGTERM
  register(shutdown)
  signal(SIGTERM, lambda x,y: exit() )
  try :
    from signal import SIGBREAK
    signal(SIGBREAK, lambda x,y : exit())
  except :
    pass
                
  loadHosts()
  s = theSocket
  s.bind(('', FOUNDRY_PORT))
  s.listen(5)
  
  theThread = threading.Thread(None, bkg, "SMDS Foundry Background Thread")
  theThread.setDaemon(True)
  theThread.start()

  while True:
    try : (c, addr) = s.accept()
    except : pass
    else :			# got a new connection
      c.setblocking(1)
      packet = smds.recvMsg(c)
      if packet :
        (msg, args) = packet
        if msg == smds.messages['QUIT'] : quit()
        if msg == 3 or msg == 7 : args = (c,)+args
        if msg > maxHandler :
          h = nullHandler
          args = (msg,)+args
        else : h = handlers[msg]
        try :
          h(*args)
        except :
          smds.msg("Bad message from %s" % addr[0], smds.MSG_WARN)
          from traceback import format_exc
          print format_exc() 
        if msg == smds.messages['RELEASE'] : quit()
      else :
        smds.msg("Bad connection from %s" % addr[0], smds.MSG_WARN)
    
    if len(myMasters) : work.set()
    else              : work.clear()
    if platform != 'win32' : cleanChildren()

if __name__ == "__main__" : run()
