#! python

# Checks to see if the SMDS Foundry is running and starts it if it isn't.
# $Id: starter.py,v 1.10 2009/04/16 21:58:04 mculbert Exp $

import smds, socket
from sys import platform, exit
from os import environ

def check_version() :
  if environ.has_key('SMDS_URL') :
    from smds.version import version
    from urllib import urlopen, urlretrieve
    import tarfile
    try:
      repository_version = urlopen(environ['SMDS_URL']+'.version').readline()
    except: return
    if version != repository_version.strip() :
      stop()
      if platform == 'win32' :
        from win32api import SearchPath
        from win32process import CreateProcess, STARTUPINFO
        (path, d) = SearchPath(None, 'python', '.exe')
        CreateProcess(path, 'python %s\\smds\\downloader.py' % 
        		environ['SMDS_LOC'],
    			None, None, 0,
                            0, None, None, STARTUPINFO())
        exit()
      else:
        try :
          fn = urlretrieve(environ['SMDS_URL'])[0]
          f = tarfile.open(fn, 'r:gz')
          for m in f.getmembers() :
            f.extract(m, environ['SMDS_LOC'])
        except : return

def start(check = True) :
  if check : check_version()
  # Attempt to connect to local foundary
  s = socket.socket()
  try : s.connect(('localhost', smds.Foundry.FOUNDRY_PORT))
  except :
    # Foundry not running, start it.
    if platform == 'win32' :
      from win32api import SearchPath
      from win32process import CreateProcess, STARTUPINFO
      (path, d) = SearchPath(None, 'python', '.exe')
      CreateProcess(path, 'python -c "from smds.Foundry import run; run()"', 
  			None, None, 0,
                          0, None, None, STARTUPINFO())
    else :
      from os import system
      system("python -m smds/Foundry > /dev/null 2>&1 &")

def stop() :
  s = socket.socket()
  try : s.connect(('localhost', smds.Foundry.FOUNDRY_PORT))
  except :
    pass
  else :
    smds.sendMsg(s, smds.messages['QUIT'])

def restart() :
  stop()
  start()

if __name__ == "__main__" : start()
