# smds.Remote module
#  The remote processor interface
#  $Id: Remote.py,v 1.8 2008/08/25 16:20:21 mculbert Exp $
 
import smds, socket, threading
from errno import EADDRINUSE
from sys import platform, exit
from select import select
from time import sleep, time

if platform == 'win32' :
  from win32process import GetCurrentProcess, SetPriorityClass, \
				IDLE_PRIORITY_CLASS
  SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS)
else :
  from os import nice
  nice(15)
  from resource import getrlimit, setrlimit, RLIMIT_DATA
  limit = getrlimit(RLIMIT_DATA)
  setrlimit(RLIMIT_DATA, (limit[1], limit[1]))

intensCache = {}
mx = threading.RLock()
waiting = None

#### Calls from the master

def recvIntens(names, intens) :
  global waiting
  from cPickle import loads
  from zlib import decompress
  for i in range(len(names)) :
    intensCache[names[i]] = loads(decompress(intens[i]))
  moreWork(*waiting)
  waiting = None

def moreWork(t, stint, solo, worktag) :
  global waiting
  t.status = smds.Task.Status_Ready
  t.solo = solo
  t.rem = stint
  names = []
  if isinstance(t.p, smds.Params.Base_intens) :
    name = t.p.intens.name
    if (name not in intensCache) : names.append(name)
    else : t.p.intens = intensCache[name]
  if isinstance(t.p, smds.Params.Base_triplet) :
    name = t.p.cef.name
    if (name not in intensCache) : names.append(name)
    else : t.p.cef = intensCache[name]
  if isinstance(t.p, smds.Params.Base_pot) :
    name = t.p.pot.name
    if (name not in intensCache) : names.append(name)
    else : t.p.pot = intensCache[name]
  if len(names) > 0 :
      smds.sendMsg(parent, smds.messages['INTENS'], names)
      waiting = (t, stint, solo, worktag)
      return
  smds.addToQueue(t, lambda task : returnTask(task, worktag))

def quit():
  mx.release()
  smds.msg("Quit!")
  parent.close()
  smds.shutdown()
  exit()

# requested seconds from now to complete
# sends RELEASE msg to master with number of bins released
def releaseWork(etoc) :
  r = smds.Engine.releaseWork(etoc)
  if r > 0 :
    if not smds.sendMsg(parent, smds.messages['RELEASE'], r) :
      mx.release()
      smds.msg("Error writing to master.")
      parent.close()
      smds.shutdown()
      exit()

#### Data sent to the master

def returnTask(t, worktag) :
  mx.acquire()
  try : smds.sendMsg(parent, smds.messages['WORK'], worktag, t.strip(True))
  except :
    mx.release()
    smds.msg("Error sending work to master.")
    parent.close()
    smds.shutdown()
    exit()
  else :
    mx.release()
  return None

lastEtoc = 0

def checkTime() :
  global lastEtoc
  etoc = smds.Engine.etoc() - time()
  if etoc <= 0 and lastEtoc > 0 :
    lastEtoc = 0
    return
  if etoc > 0 and abs(smds.Engine.etoc()-lastEtoc) > 90.0 :
    smds.msg("Notifying master of etoc %f" % etoc, smds.MSG_DBUG)
    if not smds.sendMsg(parent, smds.messages['ETOC'], 
					etoc, smds.Engine.rate()) :
      mx.release()
      smds.msg("Error writing to master.", smds.MSG_WARN)
      parent.close()
      smds.shutdown()
      exit()
    lastEtoc = smds.Engine.etoc()

#### Initialization routines

def nullHandler(*args) :
  smds.msg("Unknown message: " + str(args), smds.MSG_WARN)

parent = socket.socket()

handlers = [
		nullHandler,		# 0
		quit,			# 1
		nullHandler,		# START, 2
		moreWork,		# WORK, 3
		recvIntens,		# INTENS, 4
                nullHandler,            # ADD, 5
                nullHandler,            # REMOVE, 6
                nullHandler,            # FOUNDRY, 7
                nullHandler,            # ETOC, 8
                releaseWork,            # RELEASE, 9
	   ]
maxHandler = 9

def run(host, port, seed) :
  smds.initialize(False, seed=seed)

  # Connect to parent
  parent.connect((host,port))

  while (True) :
    mx.acquire()
    (r,w,e) = select([parent],[],[],0)
    if len(r) :
      packet = smds.recvMsg(parent)
      if not packet :
        mx.release()
        parent.close()
        smds.shutdown()
        exit("Error reading from parent.")
      (msg, args) = packet
      if (msg > maxHandler) :
        h = nullHandler
        args = (msg,)+args
      else : h = handlers[msg]
      h(*args)
    checkTime()
    mx.release()
    sleep(1)

if __name__ == "__main__" :
  from sys import argv
  if len(argv) != 4 :
    stderr.write("Usage: Remote.py host port seed\n");
    exit("Usage: Remote.py host port seed")
  run(argv[1], int(argv[2]), int(argv[3]))
