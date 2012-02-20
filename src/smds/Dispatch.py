# smds.Dispatch module
#  Remote processing engine
#  $Id: Dispatch.py,v 1.27 2009/04/17 14:42:55 mculbert Exp $

import threading, smds, socket
from smds.Task import intensCache
from smds.Foundry import FOUNDRY_PORT
from numpy import zeros
from select import select
from time import time, localtime, strftime
from traceback import format_exc
from errno import EADDRINUSE

# Private
mx = threading.RLock()
clear = threading.Event()
clear.set()
work = threading.Condition(mx)
process = threading.Condition(mx)
working = 0
unassigned = 0
queue = []
out_queue = []
callbacks = {}
bees = []		# { host, socket, task }
hosts = []		# name
waiting = {}		# name -> time; waiting for startRemote completion
worktag = 0
workWaitTime = 0
theThread = None
theProcThread = None
theFoundry = None	# socket
myPort = None
mySocket = None
halt = threading.Event()
intensCompressedCache = {}

#### Initialization functions

def initialize() :
  global theThread, theProcThread, theFoundry, myPort, mySocket
  from os import getpid
  if theFoundry == None :
    # Attempt to connect to local foundary
    s = socket.socket()
    try : s.connect(('localhost', FOUNDRY_PORT))
    except : return False
    if not smds.sendMsg(s, smds.messages['WORK'], getpid()) :
      return False
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
    mySocket = s

    theProcThread = threading.Thread(None, theProcessor, 
    					"SMDS Dispatch Postprocessor")
    theProcThread.setDaemon(True)
    theProcThread.start()
    theThread = threading.Thread(None, theDispatcher, "SMDS Dispatch Engine")
    theThread.setDaemon(True)
    theThread.start()
  return True

def shutdown() :
  global theFoundry, intensCompressedCache
  if theFoundry :
    halt.set()
    mx.acquire()
    work.notify()
    process.notify()
    for b in bees :
      smds.sendMsg(b['socket'], smds.messages['QUIT'])
      b['socket'].close()
    mx.release()
    theThread.join()
    theProcThread.join()
    mySocket.close()
    theFoundry.close()
    theFoundry = None
    intensCompressedCache = {}

#### External interface

def addToQueue(task, callback) :
  """addToQueue(task, callback)\nAdds the given task to the work queue.
  callback is the function accepting a single argument (the task) to be
  called when the Task is complete.  Note that the task's status may indicate
  an error occured when the callback is executed. """

  global callbacks, unassigned
  if not isinstance(task, smds.Task.Task) :
    raise TypeError, 		\
    	"Only objects of type smds.Task.Task can be added to the queue."
  if not callable(callback) :
    raise TypeError, "Callback is not callable in Dispatch.addToQueue."
  if (task.status == smds.Task.Status_Done) :
    callback(task)
    return
  if (task.status != smds.Task.Status_Ready) :
    raise ValueError, "Task <%s> is not ready." % (task.ID,)
  if smds.theCache and smds.theCache.fetch(task) :
    callback(task)
    return

  mx.acquire()
  if task in queue or task in out_queue :
    smds.msg("Task <%s> already on queue." % task.ID)
  else :
    smds.msg("Adding task <%s> to queue." % task.ID)
    task.status = smds.Task.Status_Waiting
    task.assigned = 0
    queue.append(task)
    callbacks[task] = callback
    unassigned += task.rem
    clear.clear()
    work.notify()
  mx.release()

def queueIsEmpty() :
  """Returns true if the queue is empty."""
  return clear.isSet()

def clearQueue() :
  """Waits for the queue to clear."""
  # Ideally we'd just clear.wait(), but then the main thread can't handle
  # keyboard interrupts and has to wait for the queue to clear.
  while not clear.isSet() : clear.wait(1)

#### Dispatch functions

seed = int(time())
def touchHosts() :
  global waiting, seed, workWaitTime
  r = False
  for host in waiting.keys() :
    if time()-waiting[host] > 60.0 :
      del waiting[host]
      hosts.append(host)
  for host in [ x for x in hosts ] :
    s = socket.socket()
    s.settimeout(2.0)
    try : s.connect((host, FOUNDRY_PORT))
    except : pass
    else :
      seed += 1
      if smds.sendMsg(s, smds.messages['START'], 
      				socket.gethostname(), myPort, seed) :
        waiting[host] = time()
        hosts.remove(host)
        workWaitTime = time()+10.0
        r = True
      s.close()
  return r

def incomingConnections() :
  r = False
  go = True
  while go :
    try : (c, addr) = mySocket.accept()
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
		      'task' : None, 'etoc'   : 0, 'rate' : 0  })
        smds.msg("New Bee: %s (%d/%d avail)" % 
        		(host, len(bees)-working, len(bees)), smds.MSG_DBUG)
  return r

def nullHandler(*args) :
  smds.msg("Unknown message: " + str(args), smds.MSG_WARN)

# mx must be locked before calling sendWork()
def sendWork() :
  global queue, working, unassigned
  if time() < workWaitTime : return False
  r = False
  avail = [ b for b in bees if b['task'] == None ]
  availPower = 0.0
  i = 0
  for b in avail :
    if b['rate'] : availPower += b['rate'] ; i += 1
  if i == 0 : availPower  = float(len(avail))
  else      : availPower *= float(len(avail))/float(i)
  if len(avail) : avgPower = availPower / float(len(avail))
    
  i = 0
  while (len(avail) > 0 and unassigned > 0) :
    r = True

    # Get next bee
    b = avail.pop()
    b['etoc'] = 0
    if b['rate'] :
      power = b['rate']/availPower
      availPower -= b['rate']
    else         : 
      power = avgPower/availPower
      availPower -= avgPower
    
    # Get next work unit
    while (queue[i].rem-queue[i].assigned == 0) : i += 1
    t = queue[i]
    if t.p.Unity : stint = t.rem
    else :
      stint = int(unassigned * power + 0.5)
      if stint < 1 : stint = 1
      if t.rem-t.assigned-stint <= 10 : stint = t.rem-t.assigned
    b['solo'] = (stint == t.rem and t.completed == 0 and t.assigned == 0)
    b['stint'] = stint

    if smds.sendMsg(b['socket'], smds.messages['WORK'],
    			t.strip(False), stint, b['solo'], worktag) :
      t.status = smds.Task.Status_Running
      t.assigned += stint
      unassigned -= stint
      b['task'] = t
      b['worktag'] = worktag
      working += 1
      smds.msg("Sent %d from <%s> to %s, %d/%d assigned; %d in queue." % 
		(stint, t.ID, b['host'], t.assigned+t.completed, 
			t.rem+t.completed, unassigned), smds.MSG_DBUG)
    else :
      remoteDead(b)
  return r

nextAdjustTime = 0

# mx must be locked before calling adjustLoads()
# Once all the work has been doled out, check the completion time estimates.
# Have machines estimated to complete after the average time relinquish
#   excess work.
def adjustLoads() :
  global nextAdjustTime
  if nextAdjustTime > time() : return False
  if unassigned > 0 or working == 0 : return False
  times = [ b['etoc'] for b in bees if b['task'] ]
  if 0 in times : return False	# must have an estimate from everyone.
  nextAdjustTime = time()+90.0
  sum = 0.0
  for x in times : sum += x
  avg = sum/len(times)
  if avg < time() : proj = 0
  else : proj = (avg-time())*working/len(bees)+time()
  a = False
  for b in bees :
    if b['task'] and b['task'].p.Unity : continue
    if b['etoc'] and b['etoc']+75 < time() :  # remote is 75s past its ETOC
      if 'bail' not in b :	# request remote to terminate
        a = True
        b['bail'] = time()
        smds.msg("Request %s to terminate now." % b['host'], smds.MSG_DBUG)
        smds.sendMsg(b['socket'], smds.messages['RELEASE'], 0)
      elif b['bail']+75 < time() :
        a = True	# remote did not respond to request to terminate
        smds.msg("%s is unresponsive." % b['host'], smds.MSG_DBUG)
        remoteDead(b)
    elif proj and b['etoc'] > proj+120.0 :
      r = proj - time()		# remote will take much longer than avg
      if r > 120.0 :
        a = True
        smds.msg("Request %s terminate in %f s." % (b['host'], r),
                          smds.MSG_DBUG)
        smds.sendMsg(b['socket'], smds.messages['RELEASE'], r)
  return a

def theDispatcher() :
  """Private."""
  mx.acquire()
  while (1) :
    if halt.isSet() :
      mx.release()
      return
    activity = touchHosts()
    activity = incomingConnections() or activity
    activity = processMessages() or activity
    activity = sendWork() or activity
    activity = adjustLoads() or activity
    if not activity : work.wait(1.0)

def remoteDead(b) :
  hosts.append(b['host'])
  reapWork(b)
  bees.remove(b)
  smds.sendMsg(b['socket'], smds.messages['QUIT'])
  b['socket'].close()
  smds.msg("Remote %s dead (%d available)." % (b['host'], len(bees)), 
		smds.MSG_DBUG)

# Given a dead bee, unassign the work assigned to it
def reapWork(b) :
  global working, unassigned
  t = b['task']
  if t == None : return
  t.assigned -= b['stint']
  smds.msg("Reaping %d from <%s>, %d/%d assigned, %d completed" % \
		(b['stint'], t.ID, t.assigned+t.completed, t.rem+t.completed, 
		t.completed), smds.MSG_DBUG)
  if t.assigned == 0 : t.status = smds.Task.Status_Waiting
  unassigned += b['stint']
  working -= 1
  b['task'] = None

def workDone(b, worktag, task) :
  global working
  if (worktag != b['worktag']) :
    smds.msg("Spurious work from remote %s" % b['host'], smds.MSG_WARN)
    return False
  if 'bail' in b : del b['bail']
  b['etoc'] = 0
  if task.status == smds.Task.Status_Error :
    smds.msg("Error executing <%s> on %s:\n%s" % 
		(b['task'].ID, b['host'], task.error), smds.MSG_DBUG)
    releaseWork(b, b['stint'])
    b['worktag'] = None
    b['task'] = None
    working -= 1
    return False

  t = b['task']
  t.assigned -= b['stint']

  if b['solo'] :
    t.run = task.run
    t.rtr = task.rtr
    t.results = task.results
    t.anal = task.anal
    t.status = smds.Task.Status_Done
  else :
    if (t.completed == 0) :
      t.run = smds.Task.Run()
      t.run.numSegments = task.run.numSegments
      if t.rtr : t.rtr.data = zeros(t.rem, 'h')
      t.results = task.results
      t.anal = task.anal
    else :
      t.run.incSegments(task.run.numSegments)
      t.run.setTime()
      t.results.cat(task.results)
      for (A, B) in [ (t.anal[i], task.anal[i]) for i in range(len(t.anal))
      				if t.anal[i].online ] :
      	A.combine(B)
    if t.rtr :
      t.rtr.data[t.completed:t.completed+b['stint']] = task.rtr.data
  
  t.completed += b['stint']
  t.rem -= b['stint']
  b['worktag'] = None
  b['task'] = None
  working -= 1
  smds.msg("Work for <%s> from %s, %d/%d completed." % 
	   (t.ID, b['host'], t.completed, t.rem+t.completed), smds.MSG_DBUG)
  
  if (t.rem == 0) :
    queue.remove(t)
    if not b['solo'] : t.status = smds.Task.Status_Processing
    out_queue.append(t)
    process.notify()
    smds.msg("Processing <%s>; %d tasks rem, %d/%d bees working" %
    		(t.ID, len(queue), working, len(bees)), smds.MSG_DBUG)
  if working and working < 4 :
    smds.msg("Waiting on: " + ' '.join([ b['host']
    			for b in bees if b['task'] ]), smds.MSG_DBUG)
  return True

def fetchIntens(b, names) :
  global intensCompressedCache
  intens = []
  for name in names :
    if name not in intensCache :
      smds.msg("Unknown intens request [%s] from %s." % (name, b['host']),
      		smds.MSG_WARN)
      return
    smds.msg("Intens [%s] request by %s" % (name, b['host']))
    if name not in intensCompressedCache :
      from cPickle import dumps
      from zlib import compress
      from exceptions import MemoryError
      try :
        intensCompressedCache[name] = compress(dumps(intensCache[name],2))
      except MemoryError :
        intensCompressedCache = {}
        intensCompressedCache[name] = compress(dumps(intensCache[name],2))
    intens.append(intensCompressedCache[name])
  if not smds.sendMsg(b['socket'], smds.messages['INTENS'], names, intens) :
    remoteDead(b)

# number of seconds until estimated completion
def reportETOC(b, etoc, rate) :
  if not b['task'] :
    smds.msg("%s reported a time to completion but isn't working on anything!"
    		% b['host'], smds.MSG_WARN )
    return
  b['etoc'] = time() + etoc
  b['rate'] = sum(rate)
  if 'bail' in b : del b['bail']
  smds.msg("ETOC for %s at %s (%s bins/s)." %
  	    (b['host'], strftime("%H:%M.%S", localtime(b['etoc'])),
  	     ', '.join([ "%f" % x for x in rate]) ),
  	    smds.MSG_DBUG)

# number of bins released
def releaseWork(b, num) :
  global working, unassigned, workWaitTime
  if 'bail' in b : del b['bail']
  t = b['task']
  if t == None :
    smds.msg("%s offered to release %d but isn't working on anything!" % 
    		(b['host'], num), smds.MSG_WARN )
    return
  if num > b['stint'] :
    smds.msg("%s offered to release %d but was only given %d!" %
    		(b['host'], num, b['stint']), smds.MSG_WARN )
    return
  t.assigned -= num
  b['stint'] -= num
  unassigned += num
  b['solo'] = False
  b['etoc'] = 0
  workWaitTime = time() + 10.0
  smds.msg(
    "Released %d of <%s> from %s, %d/%d assigned, %d completed" %
    (num, t.ID, b['host'], t.assigned+t.completed, t.rem+t.completed, 
	t.completed), smds.MSG_DBUG )

##### Processing functions

def theProcessor() :
  global out_queue
  mx.acquire()
  while (1) :
    if (len(out_queue) == 0 and len(queue) == 0) : clear.set()
    if (len(out_queue) == 0) : process.wait()
    if halt.isSet() :
      mx.release()
      return
    t = out_queue[0]
    out_queue = out_queue[1:]
    callback = callbacks.pop(t)
    mx.release()
    
    if (t.status == smds.Task.Status_Processing) :
      for a in t.anal :
        if a.online : a.finalize()
        else        : a.analyze(t)
    if smds.theCache : smds.theCache.put(t)
    try : callback(t)
    except :
      smds.msg(
        "SMDS Engine caught exception executing callback %s for Task <%s>:\n"
        % (str(callback), t.ID) + format_exc(), smds.MSG_WARN )
    
    mx.acquire()

##### Foundry functions

def addHost(host) :
  # Look to see if we already know about this one
  for x in hosts + waiting.keys() + [ b['host'] for b in bees ] :
    if x == host :
      return
  hosts.append(host)
  smds.msg("Adding host %s" % host, smds.MSG_DBUG)

def removeHost(host) :
  smds.msg("Host %s down" % host, smds.MSG_DBUG)
  for b in bees :
    if b['host'] == host :
      remoteDead(b)
  if host in waiting : del waiting[host]
  if host in hosts   : hosts.remove(host)

handlers = [
		nullHandler,		# 0
		nullHandler,		# QUIT, 1
		nullHandler,		# START, 2
		workDone,		# WORK, 3
		fetchIntens,		# INTENS, 4
		addHost,		# ADD, 5
		removeHost,		# REMOVE, 6
		nullHandler,		# FOUNDRY, 7
                reportETOC,		# ETOC, 8
                releaseWork,		# RELEASE, 9
	   ]
maxHandler = 9

def processMessages() :
  activity = False
  while 1 :
    (r,w,e) = select([theFoundry]+[ b['socket'] for b in bees ], [],[],0)
    for s in r :
      packet = smds.recvMsg(s)
      b = None
      if s != theFoundry :
        try : b = [ b for b in bees if b['socket'] == s ][0]
        except :
          smds.msg("Trashing message from %s." % s.getpeername()[0],
          		smds.MSG_WARN)
          continue
      if not packet :
        if s == theFoundry :
          smds.msg("Local foundry went down!", smds.MSG_WARN)
          from thread import interrupt_main
          interrupt_main()
          clear.set()
          halt.set()
          mx.release()
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
    if len(r) : activity = True
    else      : break
  return activity
