#!/usr/bin/python3


import numpy as np
import time
from timeit import default_timer as timer
import numba as nu
from numba import vectorize, float64, float32

import matplotlib.pyplot as plot
import matplotlib as mpl
from matplotlib import gridspec

from ctypes import * 

from functools import partial


import fjcBase as gn
import plinkstate as ps
from math import inf


v = ps.getVars()
c = gn.getC()

debug = False

def setPresetGeometry(preset='cellulase'):
    if preset == 'cellulase' :
        c.setRDockableStart( -4.2 , 2.5 , 0.0 )
        c.setRDockableEnd(0.0, 1.0, 0.0)
        c.setExcludeCentre(-1.8,0.3,0.0)
        c.setTetherExclude(2.25)
    elif preset == 'null' :
        c.setRDockableStart( 0.0 , 0.0 , 0.0 )
        c.setRDockableEnd(0.0, 0.0, 0.0)
        c.setExcludeCentre(0.0,0.0,0.0)
        c.setTetherExclude(0.0001)
    else :
        print('ERROR: Unknown preset!')
        exit(1)

def version():  
    print("Version: FJCQ-3.0-b127")

def initVars():    
    
    global c
    c = gn.getC()
    
    global v
    v = ps.getVars()
    
    v.initVars()

    global debug
    debug = False
    
    

# def phiRT(R): return c.phiC1_T * np.exp(np.dot(R, R) / c.phiC2_T)
# def phiRTD(R): return c.phiC1_TD * np.exp(np.dot(R, R) / c.phiC2_TD)
# def phiRF(R): return c.phiC1_F * np.exp(np.dot(R, R) / c.phiC2_T)
# def phiRS(R): return c.phiC1_S * np.exp(np.dot(R, R) / c.phiC2_T)

# def forceRT(R): 
#    F = nar([c.forceX, c.forceY, c.forceZ], dtype=nda)
#    forcedot = np.dot(R, F) / c.KBT 
#    return np.exp(forcedot)

# def forceRTD(R): 
#    F = nar([c.forceX, c.forceY, c.forceZ], dtype=nda)
#    RD = nar([-c.LD, c.LD, 0.0], dtype=nda)
#    R = R - RD
#    forcedot = np.dot(R, F) / c.KBT 
#    return np.exp(forcedot)

# def forceD(R):
#    # this is a dumm
#    return 1.0

'''
Cellulose binding targets:

1.037 nm    [GLC]19:c.O1 #4034    [GLC]17:c.O1 #3992        
2.068 nm    [GLC]19:c.O1 #4034    [GLC]15:c.O1 #3950        
3.105 nm    [GLC]19:c.O1 #4034    [GLC]13:c.O1 #3908        
4.133 nm    [GLC]19:c.O1 #4034    [GLC]11:c.O1 #3866        
5.158 nm    [GLC]19:c.O1 #4034    [GLC]9:c.O1 #3824        
6.158 nm    [GLC]19:c.O1 #4034    [GLC]7:c.O1 #3782        
7.177 nm    [GLC]19:c.O1 #4034    [GLC]5:c.O1 #3740        
8.208 nm    [GLC]19:c.O1 #4034    [GLC]3:c.O1 #3698        
9.244 nm    [GLC]19:c.O1 #4034    [GLC]1:c.O1 #3655        

---------------------------------------------------
CBH1 - Ry0     2.5
CBH1 - Rx0     4.2
CBH1 - RyDock  1.0
CBH1 - DOMx    4.2
CBH1 - DOMy    2.5
CBH1 - DOMz    2.5

'''

def configureMPL():
    mpl.rcParams['axes.titlesize']   = 24
    mpl.rcParams['axes.titleweight'] = 'bold'
    mpl.rcParams['axes.labelweight'] = 'bold'
    mpl.rcParams['xtick.labelsize']  = 'large'
    mpl.rcParams['ytick.labelsize']  = 'large'
    mpl.rcParams['axes.labelsize']   = 'large'
    mpl.rcParams['legend.fontsize']  = 'large'
    
#     print(mpl.rcParams)
    
    pass

def calcMarkSize( ms = 7):
    return int( 2 ** ms )

def getLimits():
    LIMITS = c.limLU
    if v.IS_SHEET_SUBSTRATE :
        LIMITS = c.limLUH
    return LIMITS

def getTimeDockBack(concUndocked_Back, concDocked_Back):
    timeDock_Back = inf
    if ( concUndocked_Back > 0.0 ) | (concDocked_Back > 0.0 ) :
        timeDock_Back = 1 / (c.kb * c.concConv * (c.PUndocked * concUndocked_Back + c.PDocked * concDocked_Back))
    
    return timeDock_Back
def plotTest0B():
    v.initVars()
    
    step = 2 * c.L / 2000.0
    colours = ['skyblue', 'green', 'tomato', 'purple', 'slategrey', 'darkolivegreen', 'brown', 'darkslateblue' ]
    
    plot.scatter(8.0,-1e-6,color='red', alpha=1.0, edgecolors='black', s=128, marker="^" )
    plot.scatter(16.0,-1e-6,color='red', alpha=1.0, edgecolors='black', s=128, marker="^" )
    
    force = 0.0
    for ND in range(0,5,1) :
        
        NTD = 8 - ND
        LD = ND * v.KUHN
        NT = 8
        NF = 8
        
        c.recalcConstants(NT, NF, ND, LD, v.KUHN, force) 
        
        
        N = c.NS
        print("N:", N, flush=True)
                
        labelS = "S_8_" + str(ND)
        
        addLabel = True
        
        xCoords = []
        yCoords = []
        colour = colours[ND]
        
        for R in (np.arange(-c.L, c.L, step)):
            
            x = R + LD
            y = gn.phiRR(R)
            
            xCoords.append(x)
            yCoords.append(y)
            
            
        plot.scatter(xCoords, yCoords , color=colour, alpha=1.0, edgecolors='none', label=labelS)
    
    plot.xlim(xmax = 17, xmin = -17)
    plot.ylim(ymax = 0.015, ymin = -0.003)
    plot.grid(True)
    plot.legend()
    plot.show()

def plotTest5B():
    v.initVars()
    
    step = 2 * c.L / 2000.0
    colours = ['skyblue', 'green', 'tomato', 'purple', 'slategrey', 'darkolivegreen', 'brown', 'darkslateblue' ]
    
    plot.scatter(8.0,-1e-6,color='red', alpha=1.0, edgecolors='black', s=128, marker="^" )
    plot.scatter(16.0,-1e-6,color='red', alpha=1.0, edgecolors='black', s=128, marker="^" )
    
    force = 0.0
    
    colourStep = -1
    for NT in range(10,50,8) :
        ND = 0
        colourStep += 1
        
        NTD = NT - ND
        LD = ND * v.KUHN
        NF = 0
        
        c.recalcConstants(NT, NF, ND, LD, v.KUHN, force) 
        
        
        N = c.NS
        print("N:", N, flush=True)
                
        labelS = "S_"+str(NT) + "_" + str(ND)
        
        addLabel = True
        
        xCoords = []
        yCoords = []
        colour = colours[colourStep]
        
        for R in (np.arange(-c.L, c.L, step)):
            
            x = R + LD
            y = gn.phiRR(R)
            
            xCoords.append(x)
            yCoords.append(y)
            
            
        plot.scatter(xCoords, yCoords , color=colour, alpha=1.0, edgecolors='none', label=labelS)
    
    plot.xlim(xmax = 17, xmin = -17)
    plot.ylim(ymax = 0.015, ymin = -0.003)
    plot.grid(True)
    plot.legend()
    plot.show()    

def plDoubleQF(RTDx, RTDy, RTDz, Rx, Ry, Rz, skew):
    R = nar([Rx, Ry, Rz])
    RTD = nar([RTDx, RTDy, RTDz])
    
    RF = R - RTD  # - c.RD
 
    distF = gn.forceRTD(RTD)
    
    fResult = (distF * gn.phiRTD_1D(RTD) / skew) * gn.phiRF_1D(RF)
    
    return fResult

def gDoubleQF(A, RTDx, RTDy, RTDz):
    
    (Rx, Ry, Rz, skew) = A
    
    R = nar([Rx, Ry, Rz], dtype=double)
    RTD = nar([RTDx, RTDy, RTDz], dtype=double)
    
    RF = R - RTD  # - c.RD
 
    distF = gn.forceRTD(RTD)
    
    fResult = (distF * gn.phiRTD(RTD) / skew) * gn.phiRF(RF)

    return fResult

# probability density at R with pulling force present
def plSkewQF(RTDx, RTDy, RTDz):
    RTD = nar([RTDx, RTDy, RTDz])
        
    fResult = gn.forceRTD(RTD) * gn.phiRTD(RTD) 
    
    return fResult

# probability density at R with undocked linker
def plSingleQU(RTx, RTy, RTz):
    RT = nar([RTx, RTy, RTz])
        
    fResult = gn.phiRT(RT) 
    
    return fResult

# probability density at R with docked linker
def plSingleQD(RTDx, RTDy, RTDz):
    RTD = nar([RTDx, RTDy, RTDz])
        
    fResult = gn.phiRTD(RTD) 
    
    return fResult


def gcalcD3B(integrand, A):
    Rx   = A[0]
    Ry   = A[1]
    Rz   = A[2]
    x    = A[3]
    skew = A[4]
            
    y    = gn.integrateP(integrand, c.limLTD, steps=30, method="gauss", args=[Rx, Ry, Rz, skew])
       
    return [ x, y ]   


@vectorize([float64(float64, float64, float64, float64)])
def gcalcD3B_O(skew, Rx, Ry, Rz):
    
    y = c.excludeMultiplier(Rx, Ry, Rz) * gn.integrateP(gDoubleQF, c.limLU, steps=35, method="gauss", args=[Rx, Ry, Rz, skew])
       
    return y




def calcDockedTime(params, debug=False):
    (NT, NF, ND, LD, v.KUHN, force, forcePick, concUndocked, concUndocked_Back, LDPerc, cUndocked1, cUndocked2, cUndocked1_Back, cUndocked2_B, colour) = params
    
    ldp = LDPerc / 100.0
    
    c.recalcConstants(NT, NF, ND, LD, v.KUHN, force, ldp)
            
    skew = gn.integrate(plSkewQF, c.limLU, 35, method="gauss")    
    print("Skew:" , skew)
    
    calcD3B_O_Skew = partial(calcD3B_O, skew)
    gcalcD3B_O_Skew = partial(gcalcD3B_O, skew)
    
    dockTarget = c.RTarget
    dockTarget_Back = c.RTarget_Back
        
    skew = 1.0
    
    print("dockTarget:" , dockTarget , " | RD:", c.RD)
    
    cDocked1 = gn.integrateP(gDoubleQF, c.lumLU, steps=35, method="gauss", args=[dockTarget[0], dockTarget[1], dockTarget[2] , skew])
    cDocked1_Back = gn.integrateP(gDoubleQF, c.lumLU, steps=35, method="gauss", args=[dockTarget_Back[0], dockTarget_Back[1], dockTarget_Back[2] , skew])
    
    timerStart2 = timer()
    _area2 = gn.integrate(gcalcD3B_O_Skew, c.lumLU, 35, method="gauss")
    timerEnd2 = timer()
    time2 = timerEnd2 - timerStart2
    cDocked2_B = _area2
    
    cDocked2 = gn.integrate(plSingleQD, c.limLTH, 35, method="gauss")
    
    concDocked = cDocked1 / cDocked2_B
    concDocked_Back = cDocked1_Back / cDocked2_B
    
    print(" ", flush=True)
    print("NTD:", c.NTD, " NF:" , c.NF, " F:", c.force, " LDP:", LDPerc)
    print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
    print("Conc. Docked       :" , concDocked, " | " , cDocked1, " ", cDocked2, " | ", cDocked2_B)
    print("Conc. Undocked     :", concUndocked, " | " , cUndocked1, " " , cUndocked2, " | ", cUndocked2_B)
    print("---------------------------------")
    print("Conc. Docked Back  :" , concDocked_Back, " | " , cDocked1_Back, " ", cDocked2_B)
    print("Conc. Undocked Back:", concUndocked_Back, " | " , cUndocked1_Back, " " , cUndocked2_B, flush=True)
    print("---------------------------------")
      
    
    timeDock = 0.0
    timeDock_Back = 0.0
    
    if (forcePick == 0) :
        print("PUndocked:", PUndocked, " PDocked:" , PDocked)
        
        timeDock = 1 / (c.kb * c.concConv * (c.PUndocked * concUndocked + c.PDocked * concDocked))
        timeDock_Back = getTimeDockBack(concUndocked_Back,concDocked_Back)
        rrr = timeDock_Back / timeDock
        
        print("TIME Dock     :" , timeDock, " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
        print("TIME Dock Back:" , timeDock_Back, " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
        print("r: ", rrr)
     
     
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print(" ", flush=True)
    
    result = [ c.RD[0], LDPerc, colour, timeDock, timeDock_Back ]
    
    print(" ", flush=True)
    
    return result

def calcDockedTimeS(params,debug=False):
    (NT, NF, ND, LD, kuhn, force, forcePick, concUndocked, concUndocked_Back, LDPerc, cUndocked1, cUndocked2, cUndocked1_Back, colour, targetStepSizeX, minTargetRangeX, maxTargetRangeX, minBackRangeX, maxBackRangeX) = params
    
    ldp = LDPerc / 100.0
    
    c.recalcConstants(NT, NF, ND, LD, kuhn, force, ldp)
    
    LIMITS = getLimits()
    
    skew = gn.integrate(plSkewQF, LIMITS, 35, method="gauss")    
    print("Skew:" , skew)
    
    calcD3B_O_Skew = partial(calcD3B_O, skew)
    gcalcD3B_O_Skew = partial(gcalcD3B_O, skew)
    
    dockTarget = c.RTarget - c.RD
    dockTarget_Back = c.RTarget_Back - c.RD
        
    #skew = 1.0
    
    cDocked1 = 0.0
    cDocked1_Back = 0.0
    for target in range( minTargetRangeX, maxTargetRangeX,  1):
        c.setRTarget( targetStepSizeX * target, 0.0 , 0.0 )
        #print("dockTarget:", c.RTarget, " | RD:", c.RD)
        cDocked1 = cDocked1 + gn.phiRTD(c.RTarget)
    
    for target in range(minBackRangeX, maxBackRangeX, 1) :
        #print("CALCULATING BACK STEP")
        c.setRTarget_Back( targetStepSizeX * target, 0.0, 0.0 )
        cDocked1_Back = cDocked1_Back + gn.phiRTD(c.RTarget_Back)
                                         
    cDocked2_B = gn.integrate(plSingleQD, LIMITS, 35, method="gauss")
    
    concDocked = cDocked1 / cDocked2_B
    concDocked_Back = cDocked1_Back / cDocked2_B
    
    if debug:
        print(" ", flush=True)
        print("NT:", NT, "ND:",ND, "LD:",c.LD, "kuhn", v.KUHN, "RD:", c.RD)
        #print("LIMITS:",LIMITS, "LU:",c.limLU)
        print("NTD:", c.NTD, " NF:" , c.NF, " F:", c.force, " LDP:", LDPerc)
        print("FORW. MIN:",minTargetRangeX,"MAX:",maxTargetRangeX," | BACK MIN:",minBackRangeX, "MAX:", maxBackRangeX )
        print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
        print("Conc. Docked       :" , concDocked, " | " , cDocked1, " ",  cDocked2_B)
        print("Conc. Undocked     :", concUndocked, " | " , cUndocked1, " " , cUndocked2)
        print("---------------------------------")
        print("Conc. Docked Back  :" , concDocked_Back, " | " , cDocked1_Back, " ", cDocked2_B)
        print("Conc. Undocked Back:", concUndocked_Back, " | " , cUndocked1_Back, " " , cUndocked2, flush=True)
        print("---------------------------------")
        
    PUndocked = 1.0 / (np.exp(2) + 1.0)
    PDocked = 1.0 - PUndocked
    
    timeDock = 0.0
    timeDock_Back = 0.0
    
    if (forcePick == 0) :
        if debug:
            print("PUndocked:", PUndocked, " PDocked:" , PDocked)
        
        timeDock = 1 / (c.kb * c.concConv * (PUndocked * concUndocked + PDocked * concDocked))
        timeDock_Back = getTimeDockBack(concUndocked_Back,concDocked_Back)
        rrr = timeDock_Back / timeDock
        
        if debug:
            print("TIME Dock     :" , timeDock, " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
            print("TIME Dock Back:" , timeDock_Back, " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
            print("r: ", rrr)
     
    if debug: 
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print(" ", flush=True)
    
    result = [ c.RD[0], LDPerc, colour, timeDock, timeDock_Back, c.NS ]
    
    if debug:
        print(" ", flush=True)
    
    return result

def plotTest3B():
    #plot.yscale('log')   
    
    stepXCount = stepXHNum
    stepYCount = stepYZHNum
    stepZCount = stepYZNum
    
    global stepTotal
    
    stepTotal = stepXCount * stepYCount * stepZCount
    
    stepx = 2 * c.L / stepXCount
    stepy = 2 * c.L / stepYCount
    stepz = 2 * c.L / stepZCount
    # mmm = 4 * c.L * c.L * c.L * np.pi / 3
    
    colours = ['skyblue', 'green', 'tomato', 'purple', 'slategrey', 'darkolivegreen', 'brown', 'darkslateblue' ]
    colours_Back = ['slategrey', 'darkolivegreen', 'brown', 'darkslateblue' ]
    fff = [ 0.000112947292395 , 0.000284316130786, 0.000569715608088, 0.000786799913971, 'NA' , 'NA' , 'NA' , 'NA' ]
    f3 = [ 0.0025444706759403645, 0.004664181408475181, 0.007471084814595963, 0.010131357362568566, 'NA' , 'NA' , 'NA' , 'NA' ]
    forces = [0.0 , 7.0e-3]
    
    for yLineVal in [10.0, 8.0, 6.0, 4.0, 2.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.08, 0.06, 0.04, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.001, 0.0008, 0.0006, 0.0004, 0.0002, 0.0001, 0.00008, 0.00006, 0.00004, 0.00002] :
#         yLineVal = 10**yLineExp
        plot.axhline(y=yLineVal, lw=1.0, ls=':', color='black')
    
    markSize = 128
    xPaper_7 = [ 0.0 ]
    yPaper_7 = [ 0.015 ]
    
    xPaper_7_1 = [ (1*v.KUHN*0.4), (1*v.KUHN) ]
    yPaper_7_1 = [ 0.015, 0.003 ]
    
    xPaper_7_2 = [ (2*v.KUHN*0.4) , (2*v.KUHN)]
    yPaper_7_2 = [ 0.012, 0.0006 ]
    
    xPaper_7_3 = [ (3*v.KUHN*0.4) , (3*v.KUHN) ]
    yPaper_7_3 = [ 0.012 , 0.0002 ]
    
    xPaper_8 = [ 0.0 ]
    yPaper_8 = [ 0.004 ]
    
    xPaper_8_1 = [ (1*v.KUHN*0.4) , (1*v.KUHN) ]
    yPaper_8_1 = [ 0.0035, 0.001 ]
    
    xPaper_8_2 = [ (2*v.KUHN*0.4) , (2*v.KUHN) ]
    yPaper_8_2 = [ 0.0030, 0.0003 ]
    
    xPaper_8_3 = [ (3*v.KUHN*0.4) , (3*v.KUHN) ]
    yPaper_8_3 = [ 0.0025, 0.0001]
    
    xPaper_8_4 = [ (4*v.KUHN*0.4) , (4*v.KUHN) ]
    yPaper_8_4 = [ 0.0020, 0.00003 ]
    
    plot.plot(xPaper_7_1, yPaper_7_1, color="black", lw=1.0)
    plot.plot(xPaper_7_2, yPaper_7_2, color="black", lw=1.0)
    plot.plot(xPaper_7_3, yPaper_7_3, color="black", lw=1.0)
    
    plot.plot(xPaper_8_1, yPaper_8_1, color="black", lw=1.0)
    plot.plot(xPaper_8_2, yPaper_8_2, color="black", lw=1.0)
    plot.plot(xPaper_8_3, yPaper_8_3, color="black", lw=1.0)
    plot.plot(xPaper_8_4, yPaper_8_4, color="black", lw=1.0)
    
    plot.scatter(xPaper_7, yPaper_7, color="lightgrey", s=markSize, edgecolors="black", label="Paper 7", marker="d")
    plot.scatter(xPaper_7_1, yPaper_7_1, color="lightgrey", s=markSize, edgecolors="black", label="Paper 7_1", marker="d")
    plot.scatter(xPaper_7_2, yPaper_7_2, color="lightgrey", s=markSize, edgecolors="black", label="Paper 7_2", marker="d")
    plot.scatter(xPaper_7_3, yPaper_7_3, color="lightgrey", s=markSize, edgecolors="black", label="Paper 7_3", marker="d")
    
    plot.scatter(xPaper_8, yPaper_8, color="dimgrey", s=markSize, edgecolors="black", label="Paper 8", marker="s")
    plot.scatter(xPaper_8_1, yPaper_8_1, color="dimgrey", s=markSize, edgecolors="black", label="Paper 8_1", marker="s")
    plot.scatter(xPaper_8_2, yPaper_8_2, color="dimgrey", s=markSize, edgecolors="black", label="Paper 8_2", marker="s")
    plot.scatter(xPaper_8_3, yPaper_8_3, color="dimgrey", s=markSize, edgecolors="black", label="Paper 8_3", marker="s")
    plot.scatter(xPaper_8_4, yPaper_8_4, color="dimgrey", s=markSize, edgecolors="black", label="Paper 8_4", marker="s")
    
    plot.grid(True, ls="-")
    #plot.show()
    
    colourSelect = 0
    for N in range(5, 9, 1) :
        NT = N
        NF = N
        
        colour = colours[colourSelect]
        colour_Back = colours_Back[colourSelect]
        
        colourSelect = colourSelect + 1
        xRes = []
        yRes = []
        xRes_Back = []
        yRes_Back = []
        
        for forcePick in [0] :
            force = forces[forcePick]
            
            ND = 0
            LD = 0.0
            
            print(" ")
            print("###########################################")
            print("### Calculating Undocked Concentration  ###")
            print("###########################################", flush=True)
            
            c.recalcConstants(NT, NF, ND, LD, v.KUHN, force)
#             gn.testC()
            
            # ## nullSkew ia undocked Z(F)
            nullSkew = gn.integrate(plSkewQF, c.limLU, 35, method="gauss")
            _nullSkew2 = gn.integrate(plSkewQF, c.limLT, 35, method="gauss")
            
            print(" ")
            print("NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
            print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
            print("nullSkew:", nullSkew, " ", _nullSkew2, flush=True)
            gcalcD3B_O_Skew = partial(gcalcD3B_O, nullSkew)
            
            nullSkew = 1.0
            print("RTarget:", c.RTarget, " | RD:", c.RD)
            
            cUndocked1 = gn.integrateP(gDoubleQF, c.lumLU, steps=35, method="gauss", args=[c.RTarget[0], c.RTarget[1], c.RTarget[2] , nullSkew])
            cUndocked1_Back = gn.integrateP(gDoubleQF, c.lumLU, steps=35, method="gauss", args=[c.RTarget_Back[0], c.RTarget_Back[1], c.RTarget_Back[2] , nullSkew])
            
            cUndocked1_E = cUndocked1
                        
            cUndocked2_B = gn.integrate(gcalcD3B_O_Skew, c.limLU, 35, method="gauss")
            cUndocked2 = gn.integrate(plSingleQU, c.limLU, 35, method="gauss")
            
            concUndocked = cUndocked1 / cUndocked2_B
            concUndocked_Back = cUndocked1_Back / cUndocked2_B
            print("Conc. Undocked     :", concUndocked, " | " , cUndocked1, " " , cUndocked2, " | " , cUndocked2_B)
            print("Conc. Undocked Back:", concUndocked_Back, " | " , cUndocked1_Back, " " , cUndocked2_B)
            
            timeNoDock = 1 / (c.kb * c.concConv * concUndocked)
            timeNoDock_Back = 1 / (c.kb * c.concConv * concUndocked_Back)
            
            print("TIME noDock     :" , timeNoDock , " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
            print("TIME noDock Back:" , timeNoDock_Back , " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
            print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            print(" ", flush=True)
            
            labelS = "Forw - " + str(N) + " F:" + str(force)
            labelS_Back = "Back - " + str(N) + " F:" + str(force)
            
            xRes.append(0.0)
            xRes_Back.append(0.0)
            
            yRes.append(timeNoDock)
            yRes_Back.append(timeNoDock_Back)
            
            dockedParams = []
            
            for ND in range(1, 5, 1):
                
                NTD = N - ND
                if (NTD < 4):
                        print("Not calculating due to increased inaccuracy")
                        break 
                
                for LDPerc in [40, 60, 80, 100] :
                    labelS = "3D DD - " + str(N) + "_" + str(ND) + " F:" + str(force)
                    
                    LD = ND * v.KUHN * LDPerc / 100.0
                    
                    dockedParams.append([NT, NF, ND, LD, v.KUHN, force, forcePick, concUndocked, concUndocked_Back, LDPerc, cUndocked1, cUndocked2, cUndocked1_Back, cUndocked2_B, colour ])
                    
                    
                    
                    
                                 
            with Pool(8) as pD3:
                resultsD3 = pD3.map(calcDockedTime, dockedParams)
            pD3.close()
            pD3.join()
            
            xRes_Lines = []
            yRes_Lines = []
            for result in resultsD3:
                (x, LDPerc, colour , timeDock, timeDock_Back) = result
                
                # we can do this because map guarantees in-order results
                if (LDPerc == 40) :
                    xRes_Lines = []
                    yRes_Lines = []
                    
                xRes.append(x)
                yRes.append(timeDock)
                
                xRes_Lines.append(x)
                yRes_Lines.append(timeDock)
                
                if LDPerc == 100 :
                    xRes_Back.append(x)
                    yRes_Back.append(timeDock_Back)
                    
                    plot.plot(xRes_Lines, yRes_Lines , color=colour, alpha=1.0, lw=1.5)
                   
            #################################################
            ### PLOT ######################################## 
            #################################################
                        
            
            plot.plot(xRes_Back, yRes_Back , color=colour_Back, alpha=0.5, lw=2)
            
            plot.scatter(xRes_Back, yRes_Back , color=colour_Back, alpha=1.0, s=markSize, marker="^", edgecolors='none', label=labelS_Back)
            plot.scatter(xRes, yRes , color=colour, alpha=1.0, s=markSize, marker="o", edgecolors='black', label=labelS)
                
        #         with Pool(8) as pD3:
        #             resultsD3 = pD3.map(calcD3B, coordList)
        #         pD3.close()
        #         pD3.join()
            print(" ")
            print(" ")
            
def getCMVal(vmin=None,vmax=None, val=None):
    if (vmin is None) or (vmax is None) or (val is None) or (vmin == vmax) :
        return None
    return (val-vmin)/(vmax-vmin)

def plotBackground(kcat_span=True):
    
    if len(v.Lower_LL) == 2:
        plot.axvspan(xmin=v.Lower_LL[0], xmax=v.Lower_LL[1] ,facecolor='lightsalmon',edgecolor='orangered', alpha = 1.0, zorder=1)
    
    if kcat_span is True :
        plot.axhspan(ymin=v.Kcat_TIMES_LOWER[0], ymax=v.Kcat_TIMES_LOWER[1] ,facecolor='darkgreen',edgecolor='darkgreen', alpha = 0.5, zorder=1)
        plot.axhline(y=v.Kcat_TIMES_LOWER[0], lw=1.0, c='darkgreen', alpha = 1.0, zorder=2)
        plot.axhline(y=v.Kcat_TIMES_LOWER[1], lw=1.0, c='darkgreen', alpha = 1.0, zorder=2)
    else :    
        for h in v.Kcat_TIMES:
            plot.axhline(y=h, lw=1.0, c='darkgreen', alpha = 1.0, zorder=2)
        
def plotCB(left=0.90,bottom=0.1,width=0.01,height=0.3, ticks=[0.0,-2.0,-4.0], ticklabels=[' 0.0', '-2.0', '-4.0']):
    
    cax = plot.axes([left, bottom, width, height])
    cbar = plot.colorbar(cax=cax,ticks=ticks)
    cbar.ax.set_yticklabels(ticklabels)
    
def plotHLine():
    plot.axhline(y=v.HLINE, lw=v.HLINE_THICKNESS, c='firebrick', alpha = 1.0, zorder = 3)

def plotTest4B(preset = None,  colourBar = True, marker='o', axs=None, ax=None, label=None, colourVal=None, cmVMIN=-4.2, cmVMAX=0.0, cmName='RdYlBu', debug=False):
    if preset is not None:
        v.initVars(preset)

    plot.yscale('log')
    
    colours = ['skyblue', 'green', 'tomato', 'purple', 'slategrey', 'darkolivegreen', 'brown', 'darkslateblue' ]
    colours_Back = ['slategrey', 'darkolivegreen', 'brown', 'darkslateblue' ]
    forces = [0.0 , 7.0e-3]
    
    
#     for yLineVal in [10.0, 8.0, 6.0, 4.0, 2.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.08, 0.06, 0.04, 0.02, 0.01, 0.008, 0.006, 0.004, 0.002, 0.001, 0.0008, 0.0006, 0.0004, 0.0002, 0.0001, 0.00008, 0.00006, 0.00004, 0.00002] :
# #         yLineVal = 10**yLineExp
#         #ax.axhline(y=yLineVal, lw=1.0, ls=':', color='black')
#         pass
    
    markSize = calcMarkSize(v.MARKSIZE)
    
    dockPerc = 0.8
    
    targetStepSizeX = v.STEP_X
    targetStepSizeY = v.STEP_Y
    
    setPresetGeometry(v.PRESETGEOM)
    
    c.setKon(v.KON, v.KON_M)
    
    colourSelect = 2
    cm = plot.cm.get_cmap(cmName)
    
    scatter = None
    
    minSegments = int( (v.MIN_TARGET_STEP_X * targetStepSizeX) / v.KUHN ) + 1 
    if v.MIN_SEGMENT_NUM > minSegments :
        minSegments = v.MIN_SEGMENT_NUM
    
    for N in range(minSegments, v.MAX_SEGMENT_NUM, v.SEGMENT_STEP) :
        NT = N
        NF = 0
        
        maxTargetRangeX = int( N * v.KUHN / targetStepSizeX ) + 1
        if maxTargetRangeX > v.MAX_TARGET_STEP_X :
            maxTargetRangeX = v.MAX_TARGET_STEP_X
        minTargetRangeX = v.MIN_TARGET_STEP_X
        minBackRangeX   = v.MIN_BACK_STEP_X
        maxBackRangeX   = v.MAX_BACK_STEP_X
        
        maxTargetRangeY = v.MAX_TARGET_STEP_Y
        if targetStepSizeY > 0 :
            maxTargetRangeY = int( N * v.KUHN / targetStepSizeY ) + 1
            
            if maxTargetRangeY > v.MAX_TARGET_STEP_Y :
                maxTargetRangeY = v.MAX_TARGET_STEP_Y
        
        
        minTargetRangeY = v.MIN_TARGET_STEP_Y
        if targetStepSizeY > 0 :
            minTargetRangeY = -1 * ( int( N * v.KUHN / targetStepSizeY ) + 1 )
        
            if minTargetRangeY < v.MIN_TARGET_STEP_Y :
                minTargetRangeY = v.MIN_TARGET_STEP_Y
            
        minTargetRangeY = v.MIN_TARGET_STEP_Y
        minBackRangeY   = v.MIN_BACK_STEP_Y
        maxBackRangeY   = v.MAX_BACK_STEP_Y
        
        colour      = colours[colourSelect]
        colour_Back = colours_Back[colourSelect]
        
#         colourSelect = colourSelect + 1
        xRes      = []
        xRes_NS   = []
        yRes      = []
        xRes_Back = []
        yRes_Back = []
        
        for forcePick in [0] :
            force = forces[forcePick]
            
            ND = 0
            LD = 0.0
            
            if debug:
                print(" ")
                print("###########################################")
                print("### Calculating Undocked Concentration  ###")
                print("###########################################", flush=True)
            
            c.recalcConstants(NT, NF, ND, LD, v.KUHN, force, v.DOCKPERCENT)
            
            LIMITS = getLimits()
            
            nullSkew = 1.0
            
            if debug:
                print(" ")
                print("NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
                print("FORW. MIN:",minTargetRangeX,"MAX:",maxTargetRangeX," | BACK MIN:",minBackRangeX, "MAX:", maxBackRangeX )
                print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
  
            
            cUndocked1 = 0.0
            for target in range( minTargetRangeX, maxTargetRangeX,  1):
                #print('YRANGE:', minTargetRangeY,maxTargetRangeY)
                for targetY in range( minTargetRangeY, maxTargetRangeY, 1):
                    TX = targetStepSizeX * target
                    TY = targetStepSizeY * targetY
                    c.setRTarget( TX, TY , 0.0 )
                    
                    if (c.RTarget[0] > c.LU) | (c.RTarget[1] > c.LU) :
                        continue
                    
                    #print("RTarget:", c.RTarget, " | RD:", c.RD)
                    
                    cUndocked1 = cUndocked1 + gn.phiRT(c.RTarget)
            
            cUndocked1_Back = 0.0
            if (maxBackRangeX - minBackRangeX) > 0 :
                for target in range( minBackRangeX, maxBackRangeX, 1 ):
                    print("CALCULATING BACK STEP")
                    for targetY in range( minBackRangeY, maxBackRangeY, 1 ):
                        TX = targetStepSizeX * target
                        TY = targetStepSizeY * targetY
                        c.setRTarget( TX, TY , 0.0 )
                        
                        if (c.RTarget[0] > c.LU) | (c.RTarget[1] > c.LU) :
                            continue
                        
                        cUndocked1_Back = cUndocked1 + gn.phiRT(c.RTarget_Back)
                                                 
            cUndocked2 = gn.integrate(plSingleQU, LIMITS, 35, method="gauss")            
            concUndocked = cUndocked1 / cUndocked2
            concUndocked_Back = cUndocked1_Back / cUndocked2
            
            if debug:
                print("Conc. Undocked     :", concUndocked, " | " , cUndocked1, " " , cUndocked2)
                print("Conc. Undocked Back:", concUndocked_Back, " | " , cUndocked1_Back, " " , cUndocked2)
            
            timeNoDock = 1 / (c.kb * c.concConv * concUndocked)
            timeNoDock_Back = getTimeDockBack(concUndocked_Back, 0.0)
            
            if debug:
                print("TIME noDock     :" , timeNoDock , " NTD:", c.NTD, " NF:" , c.NF, " F:", c.force)
                print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
                print(" ", flush=True)
            
            labelS = "Forw - " + str(N) + " F:" + str(force)
            labelS_Back = "Back - " + str(N) + " F:" + str(force)
            
#             xRes.append(c.RD[0])
#             print('XRES', xRes, c.RD[0])
#             
#             xRes_Back.append(c.RD[0])
#             
#             if ( v.USE_AA ):
#                 aa = v.calcResiduesFromSegments(c.NS)
#                 xRes_NS.append(aa)
#                 print("res:", aa, "NS:", c.NS)
#             else :
#                 xRes_NS.append(c.NS)
#             
#             yRes.append(timeNoDock)
#             yRes_Back.append(timeNoDock_Back)
            
            dockedParams = []
            resultsD3 = []
            
            #the possible range for kuhn segments docking, given the size of the
            #domain and 0.8 dockPerc efficiency
            if debug:
                print("DOCKLIST:", v.DOCKLIST) 
            for ND in v.DOCKLIST:
                
                NTD = N - ND
                if (NTD < 4):
                        print("Not calculating due to increased inaccuracy")
                        break 
                
                for LDPerc in [80] :
                    labelS = "3D DD - " + str(N) + "_" + str(ND) + " F:" + str(force)
                    
                    LD = ND * v.KUHN * LDPerc / 100.0
                    #                   (NT, NF, ND, LD, v.KUHN, force, forcePick, concUndocked, concUndocked_Back, LDPerc, cUndocked1, cUndocked2, cUndocked1_Back, colour) = params
                    #dockedParams.append([NT, NF, ND, LD, v.KUHN, force, forcePick, concUndocked, concUndocked_Back, LDPerc, cUndocked1, cUndocked2, cUndocked1_Back, colour, targetStepSizeX ])
                    dockedParams = [NT, NF, ND, LD, v.KUHN, force, forcePick, concUndocked, concUndocked_Back, LDPerc, cUndocked1, cUndocked2, cUndocked1_Back, colour, targetStepSizeX, minTargetRangeX, maxTargetRangeX, minBackRangeX, maxBackRangeX ]
                    resultsD3.append(calcDockedTimeS(dockedParams))
                    
            
            xRes_Lines = []
            yRes_Lines = []
            if debug:
                print("RES:",len(resultsD3), resultsD3)
            for result in resultsD3:
                (x, LDPerc, colour , timeDock, timeDock_Back, x_ns) = result
                
                # we can do this because map guarantees in-order results
                if (LDPerc == 40) :
                    xRes_Lines = []
                    yRes_Lines = []
                    
                xRes.append(x)
                if debug:
                    print('XRES_D', xRes, x)
                if ( v.USE_AA ):
                    aa = v.calcResiduesFromSegments(x_ns)
                    xRes_NS.append(aa)
                    if debug:
                        print("res:", aa, "NS:", x_ns)
                else :
                    xRes_NS.append(x_ns)
                
                yRes.append(timeDock)
                
                xRes_Lines.append(x)
                yRes_Lines.append(timeDock)
                
                if LDPerc == 100 :
                    xRes_Back.append(x)
                    yRes_Back.append(timeDock_Back)
                    
                    #ax.plot(xRes_Lines, yRes_Lines , cmap=cm, vmin=1e-10, vmax = 1e4, alpha=1.0, lw=1.5)
                    #ax.plot(xRes_Lines, yRes_Lines , cmap=cm, alpha=1.0, lw=1.5)
                   
            #################################################
            ### PLOT ######################################## 
            #################################################
                        
            
            #ax.plot(xRes_Back, yRes_Back , color=colour_Back, alpha=0.5, lw=2)
            
            #ax.scatter(xRes_Back, yRes_Back , color=colour_Back, alpha=1.0, s=markSize, marker="^", edgecolors='none', label=labelS_Back)
            
            if debug:
                print('XRES_P', xRes)
            
            norm = mpl.colors.Normalize(vmin=-4.2, vmax=0.0)
            
            if label is not None :
                labelS = label
            
            scatter = None
            if colourVal is not None:
#                 print("Using colour: {}".format(cm(colourVal)))
                if debug:
                    print("Using colour: {}".format(cm(colourVal)))
                xRes = np.asarray(xRes) + colourVal
                scatter = plot.scatter(xRes_NS, yRes , c=xRes, cmap=cm, vmin=cmVMIN, vmax=cmVMAX, alpha=1.0, s=markSize, marker=marker, edgecolors='black', label=labelS, zorder=10)
            else :
                scatter = plot.scatter(xRes_NS, yRes , cmap=cm, c=xRes, norm=norm, alpha=1.0, s=markSize, marker=marker, edgecolors='black', label=labelS, zorder=10)
                
        #         with Pool(8) as pD3:
        #             resultsD3 = pD3.map(calcD3B, coordList)
        #         pD3.close()
        #         pD3.join()
            if debug:
                print(" ")
                print(" ")
                
    
    #cb1 = mpl.colorbar.ColorbarBase(xRes, cmap=cm, norm=norm, orientation='vertical')
    
    if colourBar :
        plot.colorbar(scatter, ax=axs, cax=ax)
    
    if v.YLIM_MIN is None or v.YLIM_MAX is None :
        plot.ylim(auto=True)
    else :
        plot.ylim((v.YLIM_MIN, v.YLIM_MAX),auto=False)
    if v.XLIM_MIN is None or v.XLIM_MAX is None :
        plot.xlim(auto=True)
    else :
        plot.xlim((v.XLIM_MIN, v.XLIM_MAX),auto=False)
    
    if v.DRAW_HLINE :
        plotHLine()
    
    plot.xlabel("linker length (AA)")
    plot.ylabel("average binding time (s)")
    
    
    
    return scatter
        
    
def plotLegend(plots = [], legendText=(), bbox_to_anchor=(1.05, 1), loc=2, apad=0.):
    legend = plot.legend(plots, legendText , bbox_to_anchor=bbox_to_anchor, loc=loc, borderaxespad=apad)    

def genPlotLegend(markers=None,legends=None,colourVals=None, cmName=None, vmin=None, vmax=None, loc=None, bbox_to_anchor=(1.05, 1), apad=0.0):
    
    markSize = calcMarkSize(v.MARKSIZE)
    cm = plot.cm.get_cmap(cmName)
    
    plots = []
    for marker, legend, colourVal in zip( markers, legends, colourVals) :
        plots.append(mpl.lines.Line2D([0], [0], marker=marker
                                      ,color='w'
                                      , label=legend
                                      , markerfacecolor=cm(getCMVal(vmin=vmin,vmax=vmax,val=colourVal))
                                      , markeredgecolor='black'
                                      , markersize=v.MARKSIZE * 2.0)) 
    plotLegend(plots=plots, legendText=legends, loc=loc, bbox_to_anchor=bbox_to_anchor, apad=apad)
    
def generateFig2():
    fig = plot.figure() 
    
    plot.subplots_adjust(bottom=0.1, right=0.88, top=0.9, hspace=0.3)
    
    #### FIG2 A
    plot.subplot(2,2,1)
    
    plot.title('A', loc='left')
    
    v.initVars(v.PRESET_001)
    v.MAX_TARGET_STEP_X = 1000
    v.removeBackround()
    
    plotTest4B(preset=None, colourBar=False)
    
    #### FIG2 B
    
    plot.subplot(2,2,2)
    
    plot.title('B', loc='left')
    
    aarange = [35,44,57,82,146,180]
    aarange = range(10,200,2)
    
    plots      = []
    markers    = [ "o",   "s",  "v",  "D"  ]
    labels     = ['0.88 nm', '1.76 nm', '3.52 nm', '7.04 nm']
    colourVals = [ 0.0,   -1.4, -2.8, -4.1 ] 
    
    for kuhn, marker, label, colourVal in zip (
         [ 0.88,  1.76, 3.52, 7.04 ]
        ,markers
        ,labels
        ,colourVals
        ):
        for aa in aarange:
            v.initVars(v.PRESET_001)
            v.KUHN = kuhn
            
            v.MAX_TARGET_STEP_X = 1000
            v.removeBackround()
            
            v.setSegment(aa=aa)
            
            
            scatter = plotTest4B(preset=None, colourBar=False, marker=marker, label=label, colourVal=colourVal, cmName='Spectral')
        plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels)
    
    
    #### FIG2 C
    plot.subplot(2,2,3)
    
    plot.title('C', loc='left')
    
    v.initVars(v.PRESET_055)
    plotBackground()
    plotTest4B(preset=None, colourBar=False, marker='o')
    
    ######
    
    plotCB(left=0.42,bottom=0.11,width=0.01, height=0.1)
    
    #### FIG2 D
    plot.subplot(2,2,4)
    
    plot.title('D', loc='left')
    
    v.initVars(v.PRESET_013)
    v.removeBackround()
    
    plots = []
    labels = ['fibril', 'sheet']
    markers = ['o','^']
    colourVals = [0.0,-4.2]
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[0], label=labels[0], colourVal=colourVals[0] )
    plots.append(scatter)
    
    v.initVars(v.PRESET_056)
    v.removeBackround()
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[1], label=labels[1], colourVal=colourVals[1] )
    plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )

    
    
    ###########
    

    #######
    
    v.makeINAME('fig2')
    
    print(v.IMAGE_FULL_JPG)
    
    fig.set_size_inches(16,9)
    
    plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print('Image saved:',v.IMAGE_FULL_JPG)
    
    plot.show()
    plot.close()
    
    


def generateFig2_2():
    fig = plot.figure() 
    
    plot.subplots_adjust(bottom=0.1, right=0.88, top=0.9, hspace=0.3)
    
    #### FIG2 A
    plot.subplot(2,2,1)
    
    plot.title('A', loc='left')
    
    v.initVars(v.PRESET_001)
    v.MAX_TARGET_STEP_X = 1000
    v.removeBackround()
    
    #plotTest4B(preset=None, ax=axs[0,0], colourBar=False)
    plotTest4B(preset=None, colourBar=False)
    
    #### FIG2 B
    
    plot.subplot(2,2,2)
    
    plot.title('B', loc='left')
    
    aarange = [35,44,57,82,146,180]
    aarange = range(10,200,2)
    
    plots      = []
    markers    = [ "o",   "s",  "v",  "D"  ]
    labels     = ['0.88 nm', '1.76 nm', '3.52 nm', '7.04 nm']
    colourVals = [ 0.0,   -1.4, -2.8, -4.1 ] 
    
    for kuhn, marker, label, colourVal in zip (
         [ 0.88,  1.76, 3.52, 7.04 ]
        ,markers
        ,labels
        ,colourVals
        ):
        for aa in aarange:
            v.initVars(v.PRESET_001)
            v.KUHN = kuhn
            
            v.MAX_TARGET_STEP_X = 1000
            v.removeBackround()
            
            v.setSegment(aa=aa)
            
            
            scatter = plotTest4B(preset=None, colourBar=False, marker=marker, label=label, colourVal=colourVal)
        plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels)
    
    #### FIG2 C
    plot.subplot(2,2,3)
    plot.title('C', loc='left')
    
    plots      = [ ]
    markers    = [ "o", "s",  "^",  "D" ]
    labels     = [ "{:05.3f} nm".format(v.CELLULOSE_STEP_X), "{:05.3f} nm".format(v.CELL_STEP3_X), "{:05.3f} nm".format(v.CELL_STEP5_X), "{:05.3f} nm".format(v.CELL_STEP7_X) ] 
    colourVals = [ 0.0, -1.4, -2.8, -4,1 ]
    
    for preset, marker, label, colourVal in zip(
         [ v.PRESET_058, v.PRESET_059, v.PRESET_060, v.PRESET_061 ]
        ,markers
        ,labels
        ,colourVals
        ) :
        v.initVars(preset)
        v.removeBackround()
    
        scatter = plotTest4B(preset=None,  colourBar=False, marker=marker, label=label, colourVal=colourVal)
        plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='RdYlBu', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots,labels)
    
    
    #### FIG2 D
    plot.subplot(2,2,4)
    
    plot.title('D', loc='left')
    
    v.initVars(v.PRESET_013)
    v.removeBackround()
    
    plots = []
    labels = ['fibril', 'sheet']
    markers = ['o','^']
    colourVals = [0.0,-4.2]
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[0], label=labels[0], colourVal=colourVals[0] )
    plots.append(scatter)
    
    v.initVars(v.PRESET_056)
    v.removeBackround()
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[1], label=labels[1], colourVal=colourVals[1] )
    plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    ###########
    
#     plotCB()

    #######
    
    v.makeINAME('fig2_2')
    
    print(v.IMAGE_FULL_JPG)
    
    fig.set_size_inches(16,9)
    
    
    plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print('Image saved:',v.IMAGE_FULL_JPG)
    
    plot.show()
    plot.close()

def generateFig2_3():
    fig = plot.figure()
    
    plot.subplots_adjust(bottom=0.1, right=0.88, top=0.9, hspace=0.3)
    
    #### FIG2 A
    plot.subplot(2,2,1)
    
    plot.title('A', loc='left')
    
    v.initVars(v.PRESET_001)
    v.MAX_TARGET_STEP_X = 1000
    v.removeBackround()
    
    plotTest4B(preset=None, colourBar=False)
    
    #### FIG2 B
    plot.subplot(2,2,2)
    plot.title('B', loc='left')
    
    aarange = [35,44,57,82,146,180]
    aarange = range(10,200,2)
    
    plots      = []
    markers    = [ "o",   "s",  "v",  "D"  ]
    labels     = ['0.88 nm', '1.76 nm', '3.52 nm', '7.04 nm']
    colourVals = [ 0.0,   -1.4, -2.8, -4.1 ] 
    
    for kuhn, marker, label, colourVal in zip (
         [ 0.88,  1.76, 3.52, 7.04 ]
        ,markers
        ,labels
        ,colourVals
        ):
        for aa in aarange:
            v.initVars(v.PRESET_001)
            v.KUHN = kuhn
            
            v.MAX_TARGET_STEP_X = 1000
            v.removeBackround()
            
            v.setSegment(aa=aa)
                        
            scatter = plotTest4B(preset=None, colourBar=False, marker=marker, label=label, colourVal=colourVal, cmName='Spectral')
        plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels)
    
    #### FIG2 C
    plot.subplot(2,2,3)
    plot.title('C', loc='left')
    
    v.initVars(v.PRESET_013)
    v.removeBackround()
       
    plots = []
    labels = ['fibril', 'sheet']
    markers = ['o','^']
    colourVals = [0.0,-4.2]
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[0], label=labels[0], colourVal=colourVals[0] )
    plots.append(scatter)
    
    v.initVars(v.PRESET_056)
    v.removeBackround()
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[1], label=labels[1], colourVal=colourVals[1] )
    plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    ###########
    
#     plotCB()
    
    v.makeINAME('fig2_3')
    
    print(v.IMAGE_FULL_JPG)
    
    fig.set_size_inches(16,9)
    
    plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print('Image saved:',v.IMAGE_FULL_JPG)
    
    plot.show()
    plot.close()
    
def generateFig2_4():
    fig = plot.figure()
    
    plot.subplots_adjust(bottom=0.1, right=0.88, top=0.9, hspace=0.3)
    
    #### FIG2 A
    plot.subplot(2,2,1)
    
    plot.title('A', loc='left')
    
    v.initVars(v.PRESET_001)
    v.MAX_TARGET_STEP_X = 1000
    v.removeBackround()
    
    plotTest4B(preset=None, colourBar=False)
    
    #### FIG2 B
    plot.subplot(2,2,2)
    plot.title('B', loc='left')
    
    aarange = [35,44,57,82,146,180]
    aarange = range(10,200,2)
    
    plots      = []
    markers    = [ "o",   "s",  "v",  "D"  ]
    labels     = ['0.88 nm', '1.76 nm', '3.52 nm', '7.04 nm']
    colourVals = [ 0.0,   -1.4, -2.8, -4.1 ] 
    
    for kuhn, marker, label, colourVal in zip (
         [ 0.88,  1.76, 3.52, 7.04 ]
        ,markers
        ,labels
        ,colourVals
        ):
        for aa in aarange:
            v.initVars(v.PRESET_001)
            v.KUHN = kuhn
            
            v.MAX_TARGET_STEP_X = 1000
            v.removeBackround()
            
            v.setSegment(aa=aa)
            
            
            scatter = plotTest4B(preset=None, colourBar=False, marker=marker, label=label, colourVal=colourVal, cmName='Spectral')
        plots.append(scatter)

    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)    
#     plotLegend(plots,labels)
    
    #### FIG2 C
    plot.subplot(2,2,3)
    plot.title('C', loc='left')
    
    plots      = [ ]
    markers    = [ "o", "s",  "^",  "D" ]
    labels     = [ "{:05.3f} nm".format(v.CELLULOSE_STEP_X), "{:05.3f} nm".format(v.CELL_STEP3_X), "{:05.3f} nm".format(v.CELL_STEP5_X), "{:05.3f} nm".format(v.CELL_STEP7_X) ] 
    colourVals = [ 0.0, -1.4, -2.8, -4,1 ]
    
    for preset, marker, label, colourVal in zip(
         [v.PRESET_058, v.PRESET_059, v.PRESET_060, v.PRESET_061 ]
        ,markers
        , labels
        ,colourVals
        ) :
        v.initVars(preset)
        v.removeBackround()
        
        print("ColourVal: {}".format(colourVal))
        scatter = plotTest4B(preset=None,  colourBar=False, marker=marker, label=label, colourVal=colourVal, cmName='Spectral')
        plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    
    
    #### FIG2 D
    plot.subplot(2,2,4)
    plot.title('D', loc='left')
    
    v.initVars(v.PRESET_013)
    v.removeBackround()
    
    plots = []
    labels = ['fibril', 'sheet']
    markers = ['o','^']
    colourVals = [0.0,-4.2]
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[0], label=labels[0], colourVal=colourVals[0] )
    plots.append(scatter)
    
    v.initVars(v.PRESET_056)
    v.removeBackround()
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[1], label=labels[1], colourVal=colourVals[1] )
    plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    ###########
    
    v.makeINAME('fig2_4')
    
    print(v.IMAGE_FULL_JPG)
    
    fig.set_size_inches(16,9)
    
    plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print('Image saved:',v.IMAGE_FULL_JPG)
    
    plot.show()
    plot.close()
    
def generateFig2_5():
    fig = plot.figure()
    bottom = 0.15
    
    plot.subplots_adjust(bottom=bottom, right=0.9, top=0.85, hspace=0.15)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    
    image = plot.imread('C:/Dropbox/processivity_images/png/8cel_crop.png')
    
    #### FIG2 A
    
    plot.subplot(gs[0])
    
    plot.title('A', loc='left')
    plot.xticks([])
    plot.yticks([])
    plot.imshow(image, interpolation='bicubic')
    
    #### FIG2 B
    
    plot.subplot(gs[1])
    
    plot.title('B', loc='left')
    
    v.initVars(v.PRESET_055)
    plotBackground()
    
    v.DRAW_HLINE = False
    
    plotTest4B(preset=None, colourBar=False, marker='o')
    
    
    #### FIG2 C
    
    #### FIG2 D

    
    ###########
    
    plotCB(bottom=bottom+0.01, left=0.86,width=0.01,height=0.2, ticks=[0.0,-1.0,-2.0,-3.0,-4.0], ticklabels=[' 0.0','-1.0','-2.0','-3.0','-4.0'])
    
    v.makeINAME('fig2_5')
    
    print(v.IMAGE_FULL_JPG)
    
    fig.set_size_inches(16,9)
    
    plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print('Image saved:',v.IMAGE_FULL_JPG)
    
    plot.show()
    plot.close()
    
def generateFig2_6():
    fig = plot.figure()
    
    plot.subplots_adjust(bottom=0.1, right=0.88, top=0.9, hspace=0.3)
    
    #### FIG2 A
    plot.subplot(1,3,1)
    
    plot.title('A', loc='left')
    
    v.initVars(v.PRESET_001)
    v.MAX_TARGET_STEP_X = 1000
    v.YLIM_MIN = 1e-6
    v.removeBackround()
    
    plotTest4B(preset=None, colourBar=False)
    
    #### FIG2 B
    plot.subplot(1,3,2)
    plot.title('B', loc='left')
    
    aarange = [35,44,57,82,146,180]
    aarange = range(10,200,2)
    
    plots      = []
    markers    = [ "o",   "s",  "v",  "D"  ]
    labels     = ['0.88 nm', '1.76 nm', '3.52 nm', '7.04 nm']
    colourVals = [ 0.0,   -1.4, -2.8, -4.1 ] 
    
    for kuhn, marker, label, colourVal in zip (
         [ 0.88,  1.76, 3.52, 7.04 ]
        ,markers
        ,labels
        ,colourVals
        ):
        for aa in aarange:
            v.initVars(v.PRESET_001)
            v.KUHN = kuhn
            
            v.MAX_TARGET_STEP_X = 1000
            v.removeBackround()
            
            v.setSegment(aa=aa)
            
            
            scatter = plotTest4B(preset=None, colourBar=False, marker=marker, label=label, colourVal=colourVal, cmName='Spectral')
            plot.subplot(1,3,2)
        plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    #### FIG2 C
    plot.subplot(1,3,3)
    plot.title('C', loc='left')
    
    v.initVars(v.PRESET_013)
    v.YLIM_MIN = 1e-6
    v.removeBackround()
    
    plots = []
    labels = ['fibril', 'sheet']
    markers = ['o','^']
    colourVals = [0.0,-4.2]
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[0], label=labels[0], colourVal=colourVals[0] )
    plots.append(scatter)
    
    v.initVars(v.PRESET_056)
    v.removeBackround()
    scatter = plotTest4B(preset=None, colourBar=False, marker=markers[1], label=labels[1], colourVal=colourVals[1] )
    plots.append(scatter)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    #### FIG2 D

    
    ###########
    
    v.makeINAME('fig2_6')
    
    print(v.IMAGE_FULL_JPG)
    
    fig.set_size_inches(16,9)
    
    plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print('Image saved:',v.IMAGE_FULL_JPG)
    
    plot.show()
    plot.close()

def testfig():
    fig = plot.figure()
    plot.subplots_adjust(bottom=0.1, right=0.88, top=0.9, hspace=0.3)
    
    axx = []
    ax1 = plot.subplot(2,2,1)
    plot.title('A', loc='left')
    plot.xticks([])
    plot.yticks([])
    axx.append(ax1)
    
    
    ax2 = plot.subplot(2,2,2)
    plot.title('B', loc='left')
    axx.append(ax2)
    
    v.initVars(v.PRESET_057)
    v.removeBackround()
    
    plotTest4B(preset=None, colourBar=False)
    
    ax3 = plot.subplot(2,2,3)
    plot.title('C', loc='left')
    axx.append(ax3)
    
    plots      = [ ]
    markers    = [ "o", "s",  "^",  "D" ]
    labels     = [ "{:05.3f} nm".format(v.CELLULOSE_STEP_X), "{:05.3f} nm".format(v.CELL_STEP3_X), "{:05.3f} nm".format(v.CELL_STEP5_X), "{:05.3f} nm".format(v.CELL_STEP7_X) ] 
    colourVals = [ 0.0, -1.4, -2.8, -4,1 ]
    
    for preset, marker, label, colourVal in zip(
         [v.PRESET_058, v.PRESET_059, v.PRESET_060, v.PRESET_061 ]
        ,markers
        , labels
        ,colourVals
        ) :
        v.initVars(preset)
        v.removeBackround()
        
        print("ColourVal: {}".format(colourVal))
        scatter = plotTest4B(preset=None,  colourBar=False, marker=marker, label=label, colourVal=colourVal, cmName='Spectral')
        plots.append(scatter)
        
#         plot.subplot(2,2,3)
    
    genPlotLegend(markers=markers,legends=labels,colourVals=colourVals, vmin=-4.2, vmax=0.0, cmName='Spectral', loc='lower right', bbox_to_anchor=None, apad=None)
#     plotLegend(plots, labels, bbox_to_anchor=None, loc="lower right", apad = None )
    
    ax4 = plot.subplot(2,2,4)
    plot.title('D', loc='left')
    axx.append(ax4)
    
    v.initVars(v.PRESET_001)
    v.removeBackround()
    
    ax = plot.gca()
    print("AX",ax)
    
    plotTest4B(preset=None, colourBar=False, axs = axx, ax=ax)
    
    ###########
    
    plotCB(left=0.42,bottom=0.11,width=0.01, height=0.1)

    #######
    
    fig.set_size_inches(16,9)
    
#     plot.tight_layout()
    
    plot.show()
    plot.close()
    

PLIST = [ 
          v.PRESET_001 
         ,v.PRESET_002
         ,v.PRESET_003
         ,v.PRESET_004
         ,v.PRESET_005
         ,v.PRESET_006
         ,v.PRESET_007
         ,v.PRESET_008
         ,v.PRESET_009
         
         ,v.PRESET_010
         ,v.PRESET_011
         ,v.PRESET_012
         ,v.PRESET_013
         ,v.PRESET_014
         ,v.PRESET_015
         ,v.PRESET_016
         ,v.PRESET_017
         ,v.PRESET_018
         ,v.PRESET_019
         
         ,v.PRESET_020
         ,v.PRESET_021
         ,v.PRESET_022
         ,v.PRESET_023
         ,v.PRESET_024
         ,v.PRESET_025
         ,v.PRESET_026
         ,v.PRESET_027
         ,v.PRESET_028
         ,v.PRESET_029
         
         ,v.PRESET_030
         ,v.PRESET_031
         ,v.PRESET_032
         ,v.PRESET_033
         ,v.PRESET_034
         ,v.PRESET_035
         ,v.PRESET_036
         ,v.PRESET_037
         ,v.PRESET_038
         ,v.PRESET_039
         
         ,v.PRESET_040
         ,v.PRESET_041
         ,v.PRESET_042
         ,v.PRESET_043
         ,v.PRESET_044
         ,v.PRESET_045
         ,v.PRESET_046
         ,v.PRESET_047
         ,v.PRESET_048
         ,v.PRESET_049
         
         ,v.PRESET_050
         ,v.PRESET_051
         ,v.PRESET_052
         ,v.PRESET_053
         ,v.PRESET_054
         
         ,v.PRESET_055
         ,v.PRESET_056
         ,v.PRESET_057
         ,v.PRESET_058
         ,v.PRESET_059
         ,v.PRESET_060
         ,v.PRESET_061
         
        
        ]

# PLIST = [ v.PRESET_033 ]        

if __name__ == '__main__':3
    
    configureMPL()
    
#     testfig()
    generateFig2()
    generateFig2_2()
    generateFig2_3()
    generateFig2_4()
    generateFig2_5()
    generateFig2_6()
    
    
    
    plot.close()
    
    for preset in PLIST :
        #plotTest0B()
        #plotTest5B()        
        # plotTest3B()
        fig, axs = plot.subplots(nrows=1,ncols=1)
        fig.set_size_inches(16,9)
        
        v.initVars(preset)
        
        plotBackground()
        
        plotTest4B(preset)
        
        #plot.legend() 
        #plot.grid(True)
        #plot.yscale('log') 
        plot.savefig(v.IMAGE_FULL_NAME, format='svg', dpi=1200, pad_inches=0.1)
        plot.savefig(v.IMAGE_FULL_JPG, format='png', dpi=120, pad_inches=0.1)
        plot.close()
    