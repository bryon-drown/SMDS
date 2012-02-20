# $Id: opt.py,v 1.1 2006/05/25 18:19:12 mculbert Exp $
#
# Modified from:
#   optimize.py module by Travis E. Oliphant
#   part of SciPy-0.3
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

import numpy.oldnumeric as Numeric
import numpy.oldnumeric.mlab as MLab
#from scipy_base import atleast_1d, eye, mgrid, argmin, zeros, shape, \
#     squeeze, isscalar, vectorize, asarray, absolute, sqrt
#import scipy_base
from numpy import asarray	# changed from scipy_base.asarray
Num = Numeric
max = MLab.max
min = MLab.min
abs = Numeric.absolute		# changed from scipy_base.absolute
import __builtin__
pymin = __builtin__.min
pymax = __builtin__.max
#_epsilon = sqrt(scipy_base.limits.double_epsilon)
_epsilon = asarray(1.4901161193847656e-08)

def fmin(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None, 
         full_output=0, disp=1, retall=0):
    """Minimize a function using the downhill simplex algorithm.

    Description:
    
      Uses a Nelder-Mead simplex algorithm to find the minimum of function
      of one or more variables.

    Inputs:

      func -- the Python function or method to be minimized.
      x0 -- the initial guess.
      args -- extra arguments for func.

    Outputs: (xopt, {fopt, iter, funcalls, warnflag})

      xopt -- minimizer of function

      fopt -- value of function at minimum: fopt = func(xopt)
      iter -- number of iterations
      funcalls -- number of function calls
      warnflag -- Integer warning flag:
                  1 : 'Maximum number of function evaluations.'
                  2 : 'Maximum number of iterations.'
      allvecs  -- a list of solutions at each iteration

    Additional Inputs:

      xtol -- acceptable relative error in xopt for convergence.
      ftol -- acceptable relative error in func(xopt) for convergence.
      maxiter -- the maximum number of iterations to perform.
      maxfun -- the maximum number of function evaluations.
      full_output -- non-zero if fval and warnflag outputs are desired.
      disp -- non-zero to print convergence messages.
      retall -- non-zero to return list of solutions at each iteration
      
      """
    x0 = asarray(x0)
    N = len(x0)
    rank = len(x0.shape)
    if not -1 < rank < 2:
        raise ValueError, "Initial guess must be a scalar or rank-1 sequence."
    if maxiter is None:
        maxiter = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    one2np1 = range(1,N+1)

    if rank == 0:
        sim = Num.zeros((N+1,),x0.dtype())
    else:        
        sim = Num.zeros((N+1,N),x0.dtype())
    fsim = Num.zeros((N+1,),'d')
    sim[0] = x0
    if retall:
        allvecs = [sim[0]]
    fsim[0] = apply(func,(x0,)+args)
#    nonzdelt = 0.05
    nonzdelt = 0.5
    zdelt = 0.00025
    for k in range(0,N):
        y = Num.array(x0,copy=1)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt

        sim[k+1] = y
        f = apply(func,(y,)+args)
        fsim[k+1] = f

    ind = Num.argsort(fsim)
    fsim = Num.take(fsim,ind)  # sort so sim[0,:] has the lowest function value
    sim = Num.take(sim,ind,0)
    
    iterations = 1
    funcalls = N+1
    
    while (funcalls < maxfun and iterations < maxiter):
#        if (max(Num.ravel(abs(sim[1:]-sim[0]))) <= xtol \
#            and max(abs(fsim[0]-fsim[1:])) <= ftol):
        if (max(Num.ravel(abs(sim[1:]-sim[0])/sim[0])) <= xtol \
            and max(abs(fsim[0]-fsim[1:]))/fsim[0] <= ftol):
            break

        xbar = Num.add.reduce(sim[:-1],0) / N
        xr = (1+rho)*xbar - rho*sim[-1]
        fxr = apply(func,(xr,)+args)
        funcalls = funcalls + 1
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1]
            fxe = apply(func,(xe,)+args)
            funcalls = funcalls + 1

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else: # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else: # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1+psi*rho)*xbar - psi*rho*sim[-1]
                    fxc = apply(func,(xc,)+args)
                    funcalls = funcalls + 1

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1]
                    fxcc = apply(func,(xcc,)+args)
                    funcalls = funcalls + 1

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma*(sim[j] - sim[0])
                        fsim[j] = apply(func,(sim[j],)+args)
                    funcalls = funcalls + N

        ind = Num.argsort(fsim)
        sim = Num.take(sim,ind,0)
        fsim = Num.take(fsim,ind)
        iterations = iterations + 1
        if retall:
            allvecs.append(sim[0])        

    x = sim[0]
    fval = min(fsim)
    warnflag = 0

    if funcalls >= maxfun:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of function evaluations has "\
                  "been exceeded."
    elif iterations >= maxiter:
        warnflag = 2
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iterations
            print "         Function evaluations: %d" % funcalls


    if full_output:
        retlist = x, fval, iterations, funcalls, warnflag
        if retall:
            retlist += (allvecs,)
    else: 
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist

