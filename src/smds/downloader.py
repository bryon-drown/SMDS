from time import sleep
from os import environ
from sys import path
from urllib import urlretrieve
import tarfile

# To prevent attempting to load c:\windows\system32\zlib.dll on windows
for x in path :
  if x.lower() == 'c:\\windows\\system32' :
    path.remove(x)
    break

sleep(5)
try :
  fn = urlretrieve(environ['SMDS_URL'])[0]
  f = tarfile.open(fn, 'r:gz')
  for m in f.getmembers() :
    f.extract(m, environ['SMDS_LOC'])
except : pass

from smds.starter import start
start(False)
