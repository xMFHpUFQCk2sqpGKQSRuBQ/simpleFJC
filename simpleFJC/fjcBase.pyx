# cython: language_level=3
cimport cython

import numpy as np
cimport numpy as np
from gaussxw import gaussxwab
import simpleFJC.integral as gi

from numba.decorators import jit, autojit
from numba import vectorize, float64, float32
from functools import partial

# turn on logging
import logging
logging.basicConfig(level=logging.DEBUG)


cdef int KBMODE_STANDARD = 0
cdef int KBMODE_1D       = 1
cdef int KBMODE_2D       = 2

cdef double BOLTZMANN_CONSTANT_kB = 1.3806504e-23

cdef class constantsFJC:
    cdef public double kuhn
    cdef public double force
    cdef public double kuhn2
    
    cdef public double pi4
    cdef public double pi2kuhn2
    cdef public double twoKuhn2
    
    cdef public int NF
    cdef public int NT
    cdef public int ND
    cdef public double LD
    
    cdef public int NTD
    cdef public int NS
    
    cdef public double KBT
    
    cdef public double tetherExclude
    cdef public double tetherExclude2
    cdef public object excludeCentre 
    
    cdef public double L
    cdef public double L2
    
    cdef public double LU
    cdef public double LU2
    
    cdef public double domain
    
    cdef public double LF
    cdef public double LF2
    
    cdef public double LS
    cdef public double LS2
       
    cdef public double LT
    cdef public double LT2
    
    cdef public double LTD
    cdef public double LTD2
        
    cdef public double limX1T
    cdef public double limX2T
    cdef public double limY1T
    cdef public double limY2T
    cdef public double limZ1T
    cdef public double limZ2T
        
    cdef public double limX1S
    cdef public double limX2S
    cdef public double limY1S
    cdef public double limY2S
    cdef public double limZ1S
    cdef public double limZ2S
    
    cdef public double limX1SU
    cdef public double limX2SU
    cdef public double limY1SU
    cdef public double limY2SU
    cdef public double limZ1SU
    cdef public double limZ2SU
    
    cdef public double phiC1_T
    cdef public double phiC1_TD
    cdef public double phiC1_F
    cdef public double phiC1_S
    
    cdef public double phiC2_T
    cdef public double phiC2_TD
    cdef public double phiC2_F
    cdef public double phiC2_S
    
    cdef public double forceX
    cdef public double forceY
    cdef public double forceZ
    
    cdef public double kb
    cdef double kbBase
    cdef double kbMolMod
    cdef public int kbMode
    cdef public double concConv
    
    cdef public double[:] RD
    cdef public double[:] RTarget
    cdef public double[:] RTarget_Back
    
    cdef public double[:] RDockableStart
    cdef public double[:] RDockableEnd
    cdef public double dockPercent
    
    cdef public double[:,:] limLTD
    cdef public double[:,:] limLT
    cdef public double[:,:] limLTH
    cdef public double[:,:] limLU
    cdef public double[:,:] limLUH
    cdef public double[:,:] limVLTD
    cdef public double[:,:] limVLT
    cdef public double[:,:] limVLTDH
    cdef public double[:,:] limVLTH
    
    cdef public double PDocked
    cdef public double PUndocked
    cdef public double dockDG
    
    def __init__(self):
        self.RD = np.asarray([0.0, 0.0, 0.0], dtype=float)
        self.concConv = 10 / 6.02214129
        self.setExcludeCentre()
        self.setTetherExclude()
        self.setKon()
        self.setDockDG()
        
        self.setRTarget()
        self.setRTarget_Back()
        
        self.setRDockableStart()
        self.setRDockableEnd()
        
        self.recalcConstants()
        
        
        
        
    def setKon(self, double kb=20, double molMod=1e-6, object mode='standard'):
        self.kbBase = kb
        self.kbMolMod = molMod
        
        self.kb = self.kbBase * (1 / self.kbMolMod) * 6.02214129 / 10
        
        if mode == 'standard' :
            self.kbMode = KBMODE_STANDARD
        elif mode == 'fibril' :
            self.kbMode = KBMODE_1D
        elif mode == 'sheet' :
            self.kbMode = KBMODE_2D
        else :
            raise ValueError("Please provide binding rate mode: 'standard', 'fibril', 'sheet' ( {} )".format(mode))
        
    def setDockDG(self, double dockDG=-2):
        self.dockDG = -1 * dockDG
        
        self.PUndocked = 1.0 / (np.exp(self.dockDG) + 1.0)
        self.PDocked = 1.0 - self.PUndocked
        
    def setRTarget(self, double x=8.0, double y=0.0, double z=0.0):
        self.RTarget = np.array([x, y, z], dtype=np.double) - self.RD
    
    def setRTarget_Back(self, double x=-8.0, double y=0.0, double z=0.0):
        self.RTarget_Back = np.array([x, y, z], dtype=np.double) - self.RD
        
    def setRDockableStart(self, double x=0.0, double y=0.0, double z=0.0):
        self.RDockableStart = np.array([x, y, z], dtype=np.double)
        
    def setRDockableEnd(self, double x=3.5, double y=0.0, double z=0.0):
        self.RDockableEnd = np.array([x, y, z], dtype=np.double)
        
    def recalcExcludeCentre(self):
        self.excludeCentre = self.excludeCentre - self.RD
        
    def setTetherExclude(self, double r=2.0):
        self.tetherExclude  = r
        self.tetherExclude2 = self.tetherExclude ** 2.0 
    
    def excludeMultiplier(self, double Rx, double Ry, double Rz):
        cdef double[:] R = np.array([Rx, Ry, Rz], dtype=np.double)
        cdef double[:] RE = np.add(R , self.RD) 
        cdef double dRE2 = np.sum(np.multiply(RE , RE))
    
        return (dRE2 >= self.tetherExclude2)
    
    def setExcludeCentre(self,double x=0.0, double y=0.0, double z=0.0):
            self.excludeCentre = np.asarray([x,y,z])
        
    cdef _calcRD(self):
        cdef double dockingLength = self.ND * self.kuhn * self.dockPercent
        
#         print("{} - {}".format(self.RDockableEnd, self.RDockableStart))
        cdef double[:] dockVector = np.substract(self.RDockableEnd, self.RDockableStart)
        cdef double maxDockLength = np.sqrt(np.square(dockVector).sum())
        
        if (dockingLength > maxDockLength):
            dockingLength = maxDockLength
        
        if maxDockLength > 0.0 :
            dockVector = np.multiply(np.divide(dockVector , maxDockLength) , dockingLength)            
        
        self.RD = np.add(self.RDockableStart , dockVector)
        self.LD = dockingLength
        
        return maxDockLength
        
        
    def recalcRD(self, double dockPercent=1.0):
        self.dockPercent = dockPercent
        self._calcRD()
        
        
    cdef double l(self):
        cdef double length = (self.NF + self.NTD) * self.kuhn
        return length
    
    cdef double lU(self):
        cdef double length = (self.NF + self.NT) * self.kuhn
        return length
       
    cdef _calculate(self):
        self.L = self.l();
        self.L2 = self.L ** 2.0
        
        self.LU = self.lU();
        self.LU2 = self.LU ** 2.0
        
        self.domain = (2 * self.L) ** 3
           
        self.LF = self.NF * self.kuhn
        self.LF2 = self.LF ** 2.0
        
        self.LS = self.NS * self.kuhn
        self.LS2 = self.LS ** 2.0
        
        self.LT = self.NT * self.kuhn
        self.LT2 = self.LT ** 2.0
        
        self.LTD = self.NTD * self.kuhn
        self.LTD2 = self.LTD ** 2.0
        
        self.limX1T = -1.0 * self.LT
        self.limX2T = self.LT
        self.limY1T = -1.0 * self.LT
        self.limY2T = self.LT
        self.limZ1T = -1.0 * self.LT
        self.limZ2T = self.LT
            
        self.limX1S = -1.0 * self.L
        self.limX2S = self.L
        self.limY1S = -1.0 * self.L
        self.limY2S = self.L
        self.limZ1S = -1.0 * self.L
        self.limZ2S = self.L
        
        self.limLTD   = np.array([ [-self.LTD, self.LTD], [-self.LTD, self.LTD], [-self.LTD, self.LTD] ], dtype=np.double)
        self.limLT    = np.array([ [-self.LT, self.LT], [-self.LT, self.LT], [-self.LT, self.LT] ], dtype=np.double)
        self.limLTH   = np.array([ [-self.LT, self.LT], [0.0, self.LT], [-self.LT, self.LT] ], dtype=np.double)
        self.limLU    = np.array([ [-self.LU, self.LU], [-self.LU, self.LU], [-self.LU, self.LU] ], dtype=np.double)
        self.limLUH   = np.array([ [-self.LU, self.LU], [0.0, self.LU], [-self.LU, self.LU] ], dtype=np.double)
        
        self.limVLTD  = np.array([ [-self.LTD, self.LTD], [-self.LTD, self.LTD], [-self.LTD, self.LTD] ], dtype=np.double)
        self.limVLT   = np.array([ [-self.LT, self.LT], [-self.LT, self.LT], [-self.LT, self.LT] ], dtype=np.double)
        self.limVLTDH = np.array([ [-self.LTD, self.LTD], [0.0, self.LTD], [-self.LTD, self.LTD] ], dtype=np.double)
        self.limVLTH  = np.array([ [-self.LT, self.LT], [0.0, self.LT], [-self.LT, self.LT] ], dtype=np.double)
        
        
        self.phiC1_T = ((3.0 / (self.pi2kuhn2 * self.NT)) ** (1.5))
        self.phiC1_TD = ((3.0 / (self.pi2kuhn2 * self.NTD)) ** (1.5)) 
        if (self.NF > 0):
            self.phiC1_F = ((3.0 / (self.pi2kuhn2 * self.NF)) ** (1.5))
        else:
            self.phiC1_F = 0
        self.phiC1_S = ((3.0 / (self.pi2kuhn2 * self.NS)) ** (1.5)) 
        
        self.phiC2_T = (self.NT * self.twoKuhn2) / -3.0 
        self.phiC2_TD = (self.NTD * self.twoKuhn2) / -3.0
        self.phiC2_F = (self.NF * self.twoKuhn2) / -3.0 
        self.phiC2_S = (self.NS * self.twoKuhn2) / -3.0 
        
        self.forceX = -1.0 * self.force
        self.forceY = self.force
        self.forceZ = 0.0
            
    cdef recalcConstants(self, int NT=8, int NF=8, int ND=0, double LD=0.0, double kuhn=0.88, double force = 0.0, double dockPercent=1.0):
        cdef double LDm
        self.kuhn = kuhn
        self.force = force
        self.kuhn2 = self.kuhn ** 2.0
        
        self.pi4 = 4.0 * np.pi
        self.pi2kuhn2 = 2.0 * np.pi * self.kuhn2
        self.twoKuhn2 = 2.0 * self.kuhn2
        
        self.NF = NF
        self.NT = NT
        self.ND = ND
        self.LD = LD
        
        
        self.dockPercent = dockPercent
        self._calcRD()
        
        self.setRTarget()
        self.setRTarget_Back()
        
        LDm = self.ND * self.kuhn
        if (self.LD > LDm) :
            logging.warn("WARNING: dock length too large!")
            self.LD = LDm
        
        self.NTD = self.NT - self.ND
        self.NS = self.NTD + self.NF
        
        # compensating with 1e18 for the nm and nN
        self.KBT = BOLTZMANN_CONSTANT_kB * 300.0 * 1e18
        
        self._calculate()

cdef public constantsFJC c = constantsFJC()

cdef double doter(np.ndarray R1, np.ndarray R2) :
    cdef object dot = np.tensordot(R1, R2, 1)
    return dot.item()


cdef calcPhi(double c1, double c2, np.ndarray R, double ML):
    cdef double[:] radE  = R - c.excludeCentre
    cdef double    radE2 = doter( radE, radE)
    cdef double    r2    = doter(R, R) 
    
#     print('radE2:', radE2 , 'TE2:', c.tetherExclude)
#     print('r2:', r2 , 'c2:', c2)
    
    return  ( radE2 > c.tetherExclude2 ) * ( ML > r2 ) *  c1 * np.exp(r2 / c2 )

cdef calcVolume(double[:] R, double ML):
    cdef double[:] radE  = R - c.excludeCentre
    cdef double    radE2 = doter( radE, radE)
    cdef double    r2    = doter(R, R) 
    
#     print('radE2:', radE2 , 'TE2:', c.tetherExclude)
#     print('r2:', r2 , 'c2:', c2)
    
    return  ( radE2 > c.tetherExclude2 ) * ( ML > r2 ) *  1.0


cdef double calcPhi_1D(double c1, double c2, np.ndarray R, double ML):
    cdef double[:] radE  = R - c.excludeCentre
    cdef double    radE2 = doter(radE, radE)
    cdef double    r2    = doter(R, R)
    
    return ( radE2 > c.tetherExclude2 ) * ( ML > r2 ) *  c1 * np.exp(r2 / c2 )

cdef phiRR(np.ndarray R)    : return c.phiC1_S * np.exp(R * R / c.phiC2_S)
cdef phiRT(np.ndarray R)    : return calcPhi(c.phiC1_T, c.phiC2_T, R, c.LT2)
cdef phiRTD(np.ndarray R)   : return calcPhi(c.phiC1_TD, c.phiC2_TD, R, c.LTD2)
cdef phiRF(np.ndarray R)    : return calcPhi(c.phiC1_F, c.phiC2_F, R, c.LF2)
cdef phiRS(np.ndarray R)    : return calcPhi(c.phiC1_S, c.phiC2_S, R, c.LS2)

cdef phiVT(np.ndarray R)    : return calcVolume(R, c.LT2)
cdef phiVTD(np.ndarray R)   : return calcVolume(R, c.LTD2)

cdef phiRR_1D(np.ndarray R) : return c.phiC1_S * np.exp(R * R / c.phiC2_S)
cdef phiRT_1D(np.ndarray R) : return calcPhi_1D(c.phiC1_T, c.phiC2_T, R, c.LT2)
cdef phiRTD_1D(np.ndarray R): return calcPhi_1D(c.phiC1_TD, c.phiC2_TD, R, c.LTD2)
cdef phiRF_1D(np.ndarray R) : return calcPhi_1D(c.phiC1_F, c.phiC2_F, R, c.LF2)
cdef phiRS_1D(np.ndarray R) : return calcPhi_1D(c.phiC1_S, c.phiC2_S, R, c.LS2)



def getC():
    return c

cpdef testC():
    cdef double _RTDx = 2.3 - c.LD 
    
    cdef double[:] _RTD = np.array([_RTDx, 0.0, 0.0], dtype=np.double)
    cdef double[:] _RF = np.array([0.7, 0.0, 0.0], dtype=np.double)
    cdef double[:] fff = phiRTD(_RTD) * phiRF(_RF)
    
    logging.info("TC: {} {} {} {} {} Fx: {}".format(c.NF, c.NT, c.NTD, c.ND , c.LD, c.forceX) ) 
    logging.info("{} - {} - {} ".format(c.phiC1_T, c.phiC1_TD , c.phiC1_F))
    logging.info("{} - {} - {} ".format(c.phiC2_T, c.phiC2_TD , c.phiC2_F))
    
    logging.info("FR: {} RTDx: {}".format(fff , _RTDx) ) 



def integrateP(f, bounds, steps=10 ** 5, method="trap", error=0, args=[]):
    pf = f
    
    if (len(args) > 0) :
        pf = partial(f, args,)
    
    return gi.integrate(pf, bounds, steps, method, error)
