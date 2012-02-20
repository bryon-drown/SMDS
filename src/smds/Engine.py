# smds.Engine module
#  The simulation engine
#  $Id: Engine.py,v 1.16 2009/04/16 22:44:52 mculbert Exp $

import smds, numpy, threading, traceback
from time import time

# Private
mx = threading.RLock()
clear = threading.Event()
clear.set()
work = threading.Condition(mx)
queue = []
working = False
theThread = None
pistons = []
halt = threading.Event()
numRates = 3

STINT_SIZE = 50000

class Piston :
  "Private."
  def __init__(self) :
    self.wake = threading.Condition(mx)
    self.work = None
    self.rates = None
    self.rate = 0
    self.etoc = 0
    self.thread = threading.Thread(None, lambda : self.go(),
                                         "SMDS Engine %d" % len(pistons))
    self.thread.setDaemon(True)
    self.thread.start()

  def go(self) :
    mx.acquire()
    while (1) :
      if halt.isSet() :
        mx.release()
        return
      if not self.work : self.wake.wait()
      if halt.isSet() :
        mx.release()
        return
      task = self.work
      mx.release()
      
      onlineAnal = [ x for x in task.anal if x.online ]
      offlineAnal = [ x for x in task.anal if not x.online ]
      for a in onlineAnal :
        a.prepare(task)
      if (task.rtr) :
        task.rtr.data = numpy.zeros(task.rem, 'h')
        workSpace = task.rtr.data
        workSpaceSize = task.rem
      else :
        workSpaceSize = min(task.rem,STINT_SIZE)
        workSpace = numpy.zeros(workSpaceSize, 'h')
      pos = 0
      trialSize = 1
      try:
        core = task.p.initCore()
        mx.acquire()
        while (task.rem > 0) :
          if halt.isSet() :
            mx.release()
            return
          if self.rate == 0 :  stint = min(task.rem,trialSize)
          else :          stint = min(task.rem,int(60.0*self.rate+0.5))
          if stint < 1 : stint = 1
          if workSpaceSize < stint :
            workSpaceSize = stint
            workSpace = numpy.zeros(workSpaceSize, 'h')
          task.rem -= stint
          task.completed += stint
          mx.release()

          start = time()
          core.run(stint, workSpace, pos)
          for a in onlineAnal :
            a.analyze(stint, workSpace, pos)
          if (task.rtr) : pos += stint
          stop = time()
          if stop > start+1 :
            thisRate = float(stint)/(stop-start)
            if self.rates != None :
              self.rates[:numRates-1] = self.rates[1:numRates]
              self.rates[numRates-1] = thisRate
            else :
              self.rates = numpy.array([thisRate]*numRates)
            self.rate = numpy.sum(self.rates)/numRates
            self.etoc = time() + float(task.rem)/self.rate
          else :	# Stint not long enough to get a good rate
            trialSize *= 10

          mx.acquire()

        self.etoc = 0
        self.rate = 0
        self.rates = None
        mx.release()
        task.results = core.getResults()
        task.status = smds.Task.Status_Processing
        for a in onlineAnal :
          a.finalize()
        
      except:
        task.status = smds.Task.Status_Error
        task.error = traceback.format_exc()
        smds.msg(
          "SMDS Engine caught exception while processsing Task <%s>:\n"
          % (task.ID,) + traceback.format_exc(), smds.MSG_WARN )
      
      task.run = smds.Task.Run()
      mx.acquire()
      self.done = self.work
      self.work = None
      work.notify()

def theEngine() :
  "Private"
  global queue, working
  for i in range(smds.numCPUs()) : pistons.append(Piston())
  mx.acquire()
  while (1) :
    # Get a new task to work on
    if halt.isSet() :
      mx.release()
      return
    if len(queue) == 0 : work.wait()
    if halt.isSet() :
      mx.release()
      return
    (task, callback) = queue[0]
    queue = queue[1:]
    working = True
    
    # Divvy task
    task.status = smds.Task.Status_Running
    if task.p.Unity :
      pistons[0].work = task
      pistons[0].wake.notify()
    else :
      stint = int(task.rem/len(pistons))
      for e in pistons :
        if task.rem > 0 :
          if stint > task.rem : stint = task.rem
          e.work = task.copy()
          e.work.rem = stint
          e.wake.notify()
          task.rem -= stint
      if task.rem > 0 :
        pistons[0].work.rem += task.rem
        task.rem = 0
    while len([ e for e in pistons if e.work ]) : work.wait()
    mx.release()
    
    # Combine results
    if task.p.Unity :
      if task.status != smds.Task.Status_Error :
        task.run = smds.Task.Run()
        for a in task.anal :
          if a.online : a.finalize()
          else        : a.analyze(task)
        task.status = smds.Task.Status_Done
    else :
      for e in pistons :
        if e.done.status == smds.Task.Status_Error :
          task.status = smds.Task.Status_Error
          task.error = e.done.error
      if task.status != smds.Task.Status_Error :
        task.status = smds.Task.Status_Processing
        task.run = smds.Task.Run(len(pistons))
        if task.rtr :
          task.rtr.data = numpy.concatenate([e.done.rtr.data for e in pistons])
        task.results = pistons[0].done.results
        task.anal = pistons[0].done.anal
        for e in pistons[1:] :
          task.solo = task.solo and e.done.solo
          task.results.cat(e.done.results)
          for (A, B) in [ (task.anal[i], e.done.anal[i])
                          for i in range(len(task.anal))
                          if task.anal[i].online ] :
            A.combine(B)
        for e in pistons : e.done = None
        if task.solo :
          for a in task.anal :
            if a.online : a.finalize()
            else        : a.analyze(task)
          task.status = smds.Task.Status_Done
    
    # Submit completed task
    if smds.theCache : smds.theCache.put(task)
    smds.msg("Task <%s> complete." % task.ID)
    try :
      callback(task)
    except :
      smds.msg(
        "SMDS Engine caught exception executing callback %s for Task <%s>:\n"
        % (str(callback), task.ID) + traceback.format_exc(), smds.MSG_WARN )

    # Done
    mx.acquire()
    if len(queue) == 0 : clear.set()
    working = False

def initialize() :
  global theThread
  if theThread == None :
    theThread = threading.Thread(None, theEngine, "SMDS Engine")
    theThread.setDaemon(True)
    theThread.start()

def shutdown() :
  global theThread
  return		# Don't bother with all this.
  if theThread :
    mx.acquire()
    halt.set()
    work.notify()
    for e in pistons : e.wake.notify()
    mx.release()
    for e in pistons : e.thread.join()
    theThread.join()
    theThread = None
  
# requested seconds from now to complete
# Used only by Remote.py
def releaseWork(rtoc) :
  if sum(rate()) == 0 : return
  released = 0
  mx.acquire()
  for e in pistons :
    t = e.work
    if not t : continue
    if t.p.Unity : return
    r = t.rem - int(rtoc*e.rate+0.5)
    if r > 0 :
      t.rem -= r
      t.solo = False
      released += r
  mx.release()
  return released

def rate() :
  return ([ e.rate for e in pistons ])

def etoc() :
  if len(pistons) : return max([ e.etoc for e in pistons ])
  else            : return 0

def addToQueue(task, callback) :
  """addToQueue(task, callback)\nAdds the given task to the work queue.
  callback is the function accepting a single argument (the task) to be
  called when the Task is complete.  Note that the task's status may indicate
  an error occured when the callback is executed. """
  if not isinstance(task, smds.Task.Task) :
    raise TypeError, 		\
    	"Only objects of type smds.Task.Task can be added to the queue."
  if not callable(callback) :
    raise TypeError, "Callback is not callable in Engine.addToQueue."
  if (task.status == smds.Task.Status_Done) :
    callback(task)
    return
  if (task.status != smds.Task.Status_Ready) :
    raise ValueError, "Task <%s> is not ready." % (task.ID,)
  if smds.theCache and smds.theCache.fetch(task) :
    callback(task)
    return
  smds.msg("Adding task <%s> to queue." % task.ID)
  task.status = smds.Task.Status_Waiting
  mx.acquire()
  queue.append((task, callback))
  clear.clear()
  work.notify()
  mx.release()

def queueIsEmpty() :
  """Returns true if the queue is empty."""
  return clear.isSet()

def clearQueue() :
  """Waits for the queue to clear."""
  while not clear.isSet() : clear.wait(1)

