# smds.genacf Module
#   Computes an ACF from a given intensity profile
#   The only public function is genacf().
# $Id: genacf.py,v 1.14 2008/01/31 03:34:11 mculbert Exp $

__all__ = [ 'genacf' ]
import smds
from smds.Calc import registerFunction, registerData

################################
###  Private variables
################################

Gdist = None
registeredIntens = []
intensPrep = {}

################################
###  Utility functions
################################

def initialize() :
  global Gdist
  try : Gdist = registerFunction(G)
  except : return False
  return True

def pad(data, (Nx, Ny, Nz)) :
  from numpy import zeros
  (x, y, z) = data.shape
  r = zeros((Nx, Ny, Nz), data.dtype())
  r[(Nx-x)/2:(Nx+x)/2, (Ny-y)/2:(Ny+y)/2, (Nz-z)/2:(Nz+z)/2] = data
  return r

def prepIntens(intens, Pad) :
  from math import ceil, log
  from numpy import repeat, transpose, reshape
  from numpy.oldnumeric.fft import fftnd
  if Pad :
    dyadicize = lambda N: 2**int(ceil(log(N)/log(2)))
    name = "__padded_" + intens.name
  else :
    dyadicize = lambda N: N+1
    name = intens.name
  if name in intensPrep : return (name, intensPrep[name])

  smds.msg("Preparing intens <%s>" % name, smds.MSG_DBUG)
  (Nx, Ny, Nz) = map(lambda x : dyadicize(2*x+1),
  		     (intens.numBins_x, intens.numBins_x, intens.numBins_z))
  r2 =  (transpose(repeat([repeat([range(-Nx/2, Nx/2)], Ny, 0)], Nz, 0),
  			(2, 1, 0))*intens.binWidth_x) ** 2	# X
  r2 += (transpose(repeat([repeat([range(-Ny/2, Ny/2)], Nx, 0)], Nz, 0),
  			(1, 2, 0))*intens.binWidth_x) ** 2	# Y
  r2 += (repeat([repeat([range(-Nz/2, Nz/2)], Ny, 0)], Nx, 0) * 	\
  			intens.binWidth_z) ** 2			# Z
  r2 = r2.astype('f')

  O = pad(intens.getData().astype('f'), (Nx, Ny, Nz))
  del intens
  Ohat = fftnd(O)

  intensPrep[name] = (O, Ohat, r2)
  return (name, intensPrep[name])

## Dummy function used in non-distributed mode
def fetchData(name) :
  return intensPrep[name]

################################
###  Core computations
################################

from numpy import zeros, exp, sum, absolute, asarray
from numpy.oldnumeric.fft import fftnd, inverse_fftnd
from math import pi

def G(intensName, D, tau) :
  def real_fftnd(a) :
    a = asarray(a)
    for i in range(len(a.shape)) :
      a = real_fft(a, axis=i)
    return a
                                                                                  
  def inverse_real_fftnd(a) :
    a = asarray(a)
    for i in range(len(a.shape)) :
      a = inverse_real_fft(a, axis=i)
    return a

  def unwrap(data) :
    (Nx,Ny,Nz) = data.shape
    r = zeros((Nx,Ny,Nz), data.dtype())
    (r[:Nx/2,:,:],r[Nx/2:,:,:]) = (data[Nx/2:,:,:],data[:Nx/2,:,:])
    (r[:,:Ny/2,:],r[:,Ny/2:,:]) = (r[:,Ny/2:,:],r[:,:Ny/2,:].copy())
    (r[:,:,:Nz/2],r[:,:,Nz/2:]) = (r[:,:,Nz/2:],r[:,:,:Nz/2].copy())
    return r

  Psi = lambda r2, D, tau : \
  	asarray(1./(4.*pi*D*tau)**(1.5),'f') * \
  		exp(asarray(-1./(4.*D*tau),'f')*r2)

  (O, Ohat, r2) = fetchData(intensName)
  return sum(sum(sum(O*unwrap(absolute(inverse_fftnd(
  			fftnd(Psi(r2, D, tau))*Ohat ))) )))

####################################
###  Interface
####################################

def genacf(intens, D, acf, pad=True, parallel=True) :
  """genacf(intens, D, acf, pad=True, parallel=True) -> acf
  Generates the predicted autocorrelation function for a species with
  diffusion coefficient D (cm^2/s) in the detection profile
  (smds.Task.Intens) given by intens.  The generated curve is stored in the
  smds.Analysis.ACF object acf using the values for tau in acf.t.  The acf
  is also returned.  Set parallel false if you want a single-processor."""

  global Gdist
  D *= 1e14
  if parallel and Gdist == None : parallel = initialize()
  (intensName, intens) = prepIntens(intens, pad)
  if parallel and intensName not in registeredIntens :
    registerData(intensName, intens)

  if parallel : acf.G = Gdist([ (intensName, D, tau) for tau in acf.t ])
  else        : acf.G = [ G(intensName, D, tau) for tau in acf.t ]

  x = [ repr(g) for g in acf.G ]
  while 'nan' in x :
    i = x.index('nan')
    del x[i]
    del acf.G[i]
    del acf.t[i]
  while None in acf.G :
    i = x.index(None)
    del acf.G[i]
    del acf.t[i]
  return acf
