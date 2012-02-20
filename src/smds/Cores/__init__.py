# smds.Cores module
#  $Id: __init__.py,v 1.7 2009/06/16 22:36:11 mculbert Exp $

CoreTypes = {}
paramsToCore = {}

from smds.Cores import c2D, c2Di, c2Dph, c2Dphi, c3D, c3Di, c3Dph, c3Dphi
from smds.Cores import c3Df, c3Db, c3Dit, c3Dv, c2D_ccd, c2D_sp
from smds.Cores import cAHL, cAHL_fp, cAHL_lef, cAHL_elef, cAHL_chelef
from smds.Cores import cAHL_sdx, cAHL_sphdx, cAHL_sdxt
from smds.Cores import cAHL_sdxp
from smds.Cores import cAHL_sdx2t, cAHL_sphdx2t, cAHL_sdxeo, cAHL_sdxeot, cAHL_sdxeo2t
from smds.Cores import cAHL_sphdxeo, cAHL_sphdxeot, cAHL_sphdxeo2t
from smds.Cores import c3Dv_static

# coretype.getParamsType() = paramsclass
# coretype.getResultsType() = resultsclass
# core.__init__(self, Params)
# core.Run(self, numBins, array, pos)
		# writes trace length numBins to array begining at pos
		# array must be at least numBins+pos in length and of type
		#   PyArray_SHORT.
		# Core will not check these initial conditions!!
# core.GetResults(self) 
