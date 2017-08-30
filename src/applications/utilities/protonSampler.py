'''
This python script will (hopefully) generate proton disordered configurations of ice-Ih, given an input of the (.xyz) coordinates of the oxygen positions in the lattice.
'''

__author__ = "Patrick Louden (plouden@nd.edu)"
__copyright__ = "Copyright (c) 2017 by the University of Notre Dame"
__license__ = "OpenMD"

import sys
import argparse
import textwrap
import numpy as np
import random
import math as math
from fractions import Fraction, gcd
from argparse import RawDescriptionHelpFormatter


# define and zero variables used in the program here.
totOxygens = 0
totHydrogens = 0
Hmat = []

oxygenList = []



def usage():
    print __doc__

class oxygenAtom:
    def __init__(self):
        self.pos_ = [0.0, 0.0, 0.0]
        self.donors_ = []
        self.acceptors_ = []
        self.isSurface_ = 'false'
        self.neighbors_ = []

    def setPos(self, pos):
        self.pos_ = pos
        
    def getPos(self):
        return self.pos_

class proton:
    def __init__(self):
        self.pos_ = [0.0, 0.0, 0.0]

    def setPos(self, pos):
        self.pos_ = pos

    def getPos(self):
        return self.pos_

    
          
def roundMe(x):
    if (x >= 0.0):
        return math.floor(x + 0.5)
    else:
        return math.ceil(x - 0.5)

def wrapMe(disp, boxl):
    scaled = 0.0
    scaled = disp / boxl
    scaled = scaled - roundMe(scaled)
    disp = scaled * boxl
    return disp

def normalize(L1):
    L2 = []
    myLength = math.sqrt(dot(L1, L1))
    for i in range(len(L1)):
        L2.append(L1[i] / myLength)
    return L2

def dot(L1, L2):
    myDot = 0.0
    for i in range(len(L1)):
        myDot = myDot + L1[i]*L2[i]
    return myDot

        
def readInputFile(inputFileName):
    inputFile = open(inputFileName,"r")
    
    line = inputFile.readline()
    totOxygens = int(line.split()[0])
    print "total number of oxygens to read in = ", totOxygens

    line = inputFile.readline()
    Hxx = float(line.split()[1].strip(';'))
    Hxy = float(line.split()[2].strip(';'))
    Hxz = float(line.split()[3].strip(';'))
    Hyx = float(line.split()[4].strip(';'))
    Hyy = float(line.split()[5].strip(';'))
    Hyz = float(line.split()[6].strip(';'))
    Hzx = float(line.split()[7].strip(';'))
    Hzy = float(line.split()[8].strip(';'))
    Hzz = float(line.split()[9].strip(';'))
    Hmat.append([Hxx, Hxy, Hxz])
    Hmat.append([Hyx, Hyy, Hyz])
    Hmat.append([Hzx, Hzy, Hzz])
    print Hmat
    
    while 1:
        line = inputFile.readline()
        if not line:
            break
        
        newOxygen = oxygenAtom()
        posVec = []
        posX = float(line.split()[1])
        posY = float(line.split()[2])
        posZ = float(line.split()[3])
        posVec.append([posX, posY, posZ])
        newOxygen.setPos(posVec)

        oxygenList.append(newOxygen)

    print "total number of oxygens read in = ", len(oxygenList)
    return (oxygenList, Hmat)


                
    
def findNeighbors(oxygenList):
    print "finding neighbors now"
    for i in range(0, len(oxygenList)):
        oxygenA = oxygenList[i]
        for j in range(i+1, len(oxygenList)):
            oxygenB = oxygenList[j]

            if (oxygenA == oxygenB):
                print "Querying same oxygen!"

            # Want to wrap x & y displacements, but leave z as is
            # since z exposes to a surface
            xDisp = oxygenA.getPos()[0][0] - oxygenB.getPos()[0][0]
            xDisp = wrapMe(xDisp,Hmat[0][0])
            yDisp = oxygenA.getPos()[0][1] - oxygenB.getPos()[0][1]
            yDisp = wrapMe(yDisp,Hmat[1][1])
            zDisp = oxygenA.getPos()[0][2] - oxygenB.getPos()[0][2]
            #zDisp = wrapMe(zDisp,Hmat[2][2])
            
            dist = np.sqrt( pow(xDisp,2) + pow(yDisp,2) + pow(zDisp,2) )

            if (dist < 3.0):
                oxygenA.neighbors_.append(j)
                oxygenB.neighbors_.append(i)

    nSurfaceAtoms = 0
    for i in range(0, len(oxygenList)):
        oxygenA = oxygenList[i]
        if ( len(oxygenA.neighbors_) < 4):
            nSurfaceAtoms = nSurfaceAtoms + 1
            oxygenA.isSurface_ = 'true'
    print "nSurfaceAtoms =", nSurfaceAtoms

    #sanity check to make sure no-one is a self-neighbor
    for i in range(0, len(oxygenList)):
        oxygenA = oxygenList[i]
        if i in oxygenA.neighbors_ :
            print i, oxygenA.neighbors_
      

def randomizeProtons(oxygenList):

            
    # select random oxygen to start with and perform protonBucketBrigades.
    for i in range(0, 1):
        startingOxygenIndex = random.randint(0,len(oxygenList)-1)
        protonBucketBrigade(oxygenList,startingOxygenIndex)

   
    '''
    # check to see if any molecules have to many or still need protons
    bulkUnderCoord = []
    surfUnderCoord = []
    for i in range(0, len(oxygenList)):
        oxygen = oxygenList[i]
        
        if (len(oxygen.donors_) > 2):
            pass
            #print "oxygen no.", i, "has", len(oxygen.donors_),"donors"
        elif (len(oxygen.acceptors_) > 2):
            pass
             #print "oxygen no.", i, "has", len(oxygen.acceptors_),"acceptors"
             
        if (oxygen.isSurface_ == 'true'):
            numProtons = len(oxygen.donors_) + len(oxygen.acceptors_)
            if (numProtons < 3):
                surfUnderCoord.append(i)
        elif (oxygen.isSurface_ == 'false'):
            numProtons = len(oxygen.donors_) + len(oxygen.acceptors_)
            if (numProtons < 4):
                bulkUnderCoord.append(i)

    print "There are", len(bulkUnderCoord), "bulk oxygens needing protons"
    print "There are", len(surfUnderCoord), "surf. oxygens needing protons"

    print "Attempting more protonBucketBridages starting with the deficient oxygens"
    for j in range(0, 10):
        for i in range(0, len(bulkUnderCoord)):
            protonBucketBrigade(oxygenList, bulkUnderCoord[i])

    for j in range(0, 10):        
        for i in range(0, len(surfUnderCoord)):
            protonBucketBrigade(oxygenList, surfUnderCoord[i])



    # re-check to see if any molecules have to many or still need protons
    bulkUnderCoord = []
    surfUnderCoord = []
    for i in range(0, len(oxygenList)):
        oxygen = oxygenList[i]
        
        if (len(oxygen.donors_) > 2):
            pass
            #print "oxygen no.", i, "has", len(oxygen.donors_),"donors"
        elif (len(oxygen.acceptors_) > 2):
            pass
             #print "oxygen no.", i, "has", len(oxygen.acceptors_),"acceptors"
             
        if (oxygen.isSurface_ == 'true'):
            numProtons = len(oxygen.donors_) + len(oxygen.acceptors_)
            if (numProtons < 3):
                surfUnderCoord.append(i)
                print 'surfUC', i, oxygen.neighbors_, oxygen.donors_, oxygen.acceptors_
        elif (oxygen.isSurface_ == 'false'):
            numProtons = len(oxygen.donors_) + len(oxygen.acceptors_)
            if (numProtons < 4):
                bulkUnderCoord.append(i)
                print 'bulkUC', i, oxygen.neighbors_, oxygen.donors_, oxygen.acceptors_
                #print 'bulkUCNeighbor',oxygen.neighbors_[0], oxygenList[oxygen.neighbors_[0]].donors_, oxygenList[oxygen.neighbors_[0]].acceptors_

    print "There are", len(bulkUnderCoord), "bulk oxygens needing protons"
    print "There are", len(surfUnderCoord), "surf. oxygens needing protons"


    numDonors = 0
    numAcceptors = 0
    for i in range(0, len(oxygenList)):
        oxygen = oxygenList[i]
        if (oxygen.isSurface_ == 'true'):
            numDonors = numDonors + len(oxygen.donors_)
            numAcceptors = numAcceptors + len(oxygen.acceptors_)
    print "numDonors =", numDonors
    print "numAcceptors =", numAcceptors

    
    for i in range(0, len(oxygenList)):
        oxygen = oxygenList[i]
        if (oxygen.isSurface_ == 'false'):
            if (len(oxygen.donors_) < 2):
                print '\nDonorDeficient'
                for i in range(0, len(oxygen.neighbors_)):
                    nborOx = oxygenList[oxygen.neighbors_[i]]
                    print "neighbor",oxygen.neighbors_[i], nborOx.neighbors_, nborOx.donors_, nborOx.acceptors_
            if (len(oxygen.acceptors_) < 2):
                print '\nAcceptorDeficient'
                print "UCoxygen", i, oxygen.neighbors_, oxygen.donors_, oxygen.acceptors_
                for i in range(0, len(oxygen.neighbors_)):
                    nborOx = oxygenList[oxygen.neighbors_[i]]
                    print "neighbor",oxygen.neighbors_[i], nborOx.neighbors_, nborOx.donors_, nborOx.acceptors_
    '''
def protonBucketBrigade(oxygenList,startOxygenIndex):
    # assign starting oxygen
    startOxygen = oxygenList[startOxygenIndex]
    
    # Going to keep track of bonding with the protonLoop, and only if successful will we implement bonds.
    protonLoop = []
    protonLoop.append(startOxygenIndex)
    
    # First step using the startingOxygen
    
    # badNeighbors are oxygens already donating to/accepting from current
    badNbors = set(startOxygen.donors_).union(set(startOxygen.acceptors_))
    # goodNbors are my nbors not in badnbors
    goodNbors = set(startOxygen.neighbors_).difference(badNbors)

    
    if (len(goodNbors) == 0):
        return 

    # sample next oxygen from good neighbors
    if (startOxygen.isSurface_ == 'true'):
        numBonds = len(startOxygen.donors_) + len(startOxygen.acceptors_)
        if (numBonds < 3):
            if ( len(startOxygen.donors_) < 2):
                    nextOxygenIndex = random.sample(goodNbors,1)[0]
                    nextOxygen = oxygenList[nextOxygenIndex]
                    
                    if (nextOxygen.isSurface_ == 'true'):
                            nOnumBonds = len(nextOxygen.donors_) + len(nextOxygen.acceptors_)
                            if (nOnumBonds < 3):
                                    if (len(nextOxygen.acceptors_) < 2):
                                        protonLoop.append(nextOxygenIndex)
                                    elif (len(nextOxygen.acceptors_) == 2):
                                            return #nextOxygen has saturated acceptors
                            elif (nOnumBonds == 3):
                                    return #nextOxygen has total saturation of bonds
                            elif (nOnumBonds > 3):
                                print "Code is broken, nextOxygen is surface and has > 3 bonds."
                                    
                    elif (nextOxygen.isSurface_ == 'false'):
                        if (len(nextOxygen.acceptors_) < 2):
                            protonLoop.append(nextOxygenIndex)
                        elif (len(nextOxygen.acceptors_) == 2):
                            return 
                        elif (len(nextOxygen.acceptors_) > 2):
                            print "Code is broken, nextOxygen is bulk and has > 2 acceptors."

            # we are a surface oxygen, but we have 2 donors already.   
            elif ( len(startOxygen.donors_) == 2):
                    protonLoop = []
                    return
                
            elif ( len(startOxygen.donors_) > 2):
                print "Code is broken, startOxygen is surface and has > 2 donors"
                
        # we are a surface oxygen, but have saturation of bonds.        
        elif (numBonds == 3):
            return
        
        elif (numBonds > 3):
            print "numBonds > 3 (while)"    
        
    elif (startOxygen.isSurface_ == 'false'):
            if (len(startOxygen.donors_) < 2): #startOxygen < 2 donors, therefore can make a bond
                nextOxygenIndex = random.sample(goodNbors,1)[0]
                nextOxygen = oxygenList[nextOxygenIndex]

                if (nextOxygen.isSurface_ == 'true'):
                    nOnumBonds = len(nextOxygen.donors_) + len(nextOxygen.acceptors_)
                    if (nOnumBonds < 3):
                        if ( len(nextOxygen.acceptors_) < 2): #requirement met, make a bond
                            protonLoop.append(nextOxygenIndex)
                        elif ( len(nextOxygen.acceptors_) == 2):
                            return # nextOxygen has 2 acceptors, no room for one more.
                        elif (nOnumBonds == 3):
                            return #nextOxygen has 3 bonds, no room for one more.
                        elif (nOnumBonds > 3):
                            print "Code is broken, nextOxygen is surface and has > 3 bonds."
                                    
                elif (nextOxygen.isSurface_ == 'false'):
                    if (len(nextOxygen.acceptors_) < 2): #requirement met, make a bond
                        protonLoop.append(nextOxygenIndex)
                    elif (len(nextOxygen.acceptors_) == 2):
                        return # nextOxygen has 2 acceptors, no room for one more.
                    elif (len(nextOxygen.acceptors_) > 2):
                        print "Code is broken, nextOxygen is bulk and has > 2 acceptors."

                
            #startOxygen is not surface and has 2 donors, no more room for additional donors.    
            elif (len(startOxygen.donors_) == 2):
                return

            elif (len(startOxygen.donors_) > 2):
                print "Code broke, startOxygen is not surface and has > 2 donors."


    # update nextOxygen to currentOxygen, begin looping.
    previousOxygen = startOxygen
    previousOxygenIndex = startOxygenIndex
    currentOxygen = nextOxygen
    currentOxygenIndex = nextOxygenIndex

    
    # Now we will loop until we reach back to the starting Oxygen, guaranteeing a zero-net-dipole chain!
    while (currentOxygenIndex != startOxygenIndex):
        
        badNbors = set(currentOxygen.donors_).union(set(currentOxygen.acceptors_))
        goodNbors = set(currentOxygen.neighbors_).difference(badNbors)


        # if currentOxygen is surface, must test previousOxygen for surface or bulk.
        if (currentOxygen.isSurface_ == 'true'):
            #can't test the raw number of donors_ could have a (DAA) set of bonds and donors < 2 would give false positive.
            numBonds = len(currentOxygen.donors_) + len(currentOxygen.acceptors_)
            if (numBonds < 3): # with < 3 bonds, currentOxygen has room to donate to another oxygen.
                nextOxygenIndex = random.sample(goodNbors,1)[0]
                nextOxygen = oxygenList[nextOxygenIndex]
                
                # if previous and current both surface, nextOxygen must be a bulk.
                if (previousOxygen.isSurface_ == 'true'):
                    if (nextOxygen.isSurface_ == 'true'):
                        return
                    elif (nextOxygen.isSurface_ == 'false'):
                        if (len(nextOxygen.acceptors_) < 2):
                            protonLoop.append(nextOxygenIndex)
                        elif (len(nextOxygen.acceptors_) == 2):
                            return # nextOxygen has 2 acceptors, no room for one more.
                        elif (len(nextOxygen.acceptors_) > 2):
                            print "Code is broken (while), nextOxygen is bulk and has > 2 acceptors."

                            
                # if previousOxygen is not surface, then nextOxygen can be a surface or bulk oxygen
                elif (previousOxygen.isSurface_ == 'false'):
                    
                    if (nextOxygen.isSurface_ == 'true'):
                        nOnumBonds = len(nextOxygen.donors_) + len(nextOxygen.acceptors_)
                        if (nOnumBonds < 3):
                            if (len(nextOxygen.acceptors_) < 2):
                                protonLoop.append(nextOxygenIndex)
                            elif (len(nextOxygen.acceptors_) == 2):
                                return # nextOxygen has 2 acceptors, no room for one more.
                            elif (len(nextOxygen.acceptors_) > 2):
                                print "Code broken (while), nextOxygen has > 2 acceptors"
                        elif (nOnumBonds == 3):
                            return #nextOxygen has 3 bonds, no room for one more.
                        elif (nOnumBonds > 3):
                            print "Code is broken (while), nextOxygen is surface and has > 3 bonds."
                            
                    elif (nextOxygen.isSurface_ == 'false'):
                        if (len(nextOxygen.acceptors_) < 2):
                            protonLoop.append(nextOxygenIndex)
                        elif (len(nextOxygen.acceptors_) == 2):
                            return # nextOxygen has 2 acceptors, no room for one more.
                        elif (len(nextOxygen.acceptors_) > 2):
                            print "Code is broken (while), nextOxygen is bulk and has > 2 acceptors."
                            
            elif (numBonds == 3):
                return #currentOxygen has a saturation of bonds.
            elif (numBonds > 3):
                print "Code is broken (while), currentOxygen is surface and has > 3 bonds."

        #if currentOxygen is not surface, nextOxygen can be bulk or surface.        
        elif (currentOxygen.isSurface_ == 'false'):
            nextOxygenIndex = random.sample(goodNbors,1)[0]
            nextOxygen = oxygenList[nextOxygenIndex]

            if (nextOxygen.isSurface_ == 'true'):
                nOnumBonds = len(nextOxygen.donors_) + len(nextOxygen.acceptors_)
                if (nOnumBonds < 3):
                    if (len(nextOxygen.acceptors_) < 2):
                        protonLoop.append(nextOxygenIndex)
                    elif (len(nextOxygen.acceptors_) == 2):
                        return # nextOxygen has 2 acceptors, no room for one more.
                    elif (len(nextOxygen.acceptors_) > 2):
                        print "Code broken (while), nextOxygen has > 2 acceptors"
                elif (nOnumBonds == 3):
                    return #nextOxygen has 3 bonds, no room for one more.
                elif (nOnumBonds > 3):
                    print "Code is broken (while), nextOxygen is surface and has > 3 bonds."
                
            elif (nextOxygen.isSurface_ == 'false'):
                if (len(nextOxygen.acceptors_) < 2):
                    protonLoop.append(nextOxygenIndex)
                elif (len(nextOxygen.acceptors_) == 2):
                    return # nextOxygen has 2 acceptors, no room for one more.
                elif (len(nextOxygen.acceptors_) > 2):
                    print "Code is broken (while), nextOxygen is bulk and has > 2 acceptors."


        # update nextOxygen to currentOxygen, begin looping.
        previousOxygen = currentOxygen
        previousOxygenIndex = nextOxygenIndex
        currentOxygen = nextOxygen
        currentOxygenIndex = nextOxygenIndex
      
    #post while loop: test for length of protonLoop, will be set to zero upon break above.
    if (len(protonLoop) > 0):
        assignBonds(oxygenList,protonLoop)





        
def assignBonds(oxygenList,protonLoop):
    print "protonLoop =", protonLoop
    # in this function, we will assign the donor/acceptor to the oxygens
    for i in range(0, len(protonLoop)-1):
        oxygenList[protonLoop[i]].donors_.append(protonLoop[i+1])
        oxygenList[protonLoop[i+1]].acceptors_.append(protonLoop[i])


        
def addProtonsToDonors(oxygenList):
    protonList = []
    for i in range(0, len(oxygenList)):
        oxygenA = oxygenList[i]
        for j in range(0, len(oxygenA.donors_)):
            oxygenB = oxygenList[oxygenA.donors_[j]]

            # place protons 1 unit vector away from the donating oxygen
            # towards the accepting oxygen
            xDisp = oxygenA.getPos()[0][0] - oxygenB.getPos()[0][0]
            yDisp = oxygenA.getPos()[0][1] - oxygenB.getPos()[0][1]
            zDisp = oxygenA.getPos()[0][2] - oxygenB.getPos()[0][2]

            rVec = [xDisp,yDisp,zDisp]
            rUnit = normalize(rVec)
           
            xPos = oxygenA.getPos()[0][0] + rUnit[0]
            yPos = oxygenA.getPos()[0][1] + rUnit[1]
            zPos = oxygenA.getPos()[0][2] + rUnit[2]
            protonVec = [xPos, yPos, zPos]
            
            newProton = proton()
            newProton.setPos(protonVec)

            protonList.append(newProton)

        # assign a -1 for any surface atom missing a donor or acceptor,
        # will reference the -1 when we compute quaternions.
        if (oxygenA.isSurface_ == 'true'):
            if (len(oxygenA.donors_) < 2):
                oxygenA.donors_.append(-1)
            elif (len(oxygenA.acceptors_) < 2):
                oxygenA.acceptors_.append(-1)
            #print oxygenA.donors_, oxygenA.acceptors_
            
    return protonList



def computeQuats(oxygenList):

    #for i in range(len(indices)):
    #    DonatingTo.append(list())
    #    AcceptingFrom.append(list())
    #    OldDonor.append(list())
    # print "Computing Quaternions"

    # Set up initial rotation matrix
    ux = [0.0, 0.0, 0.0]
    uy = [0.0, 0.0, 0.0]
    uz = [0.0, 0.0, 0.0]
    RotMat = [ux, uy, uz]
    totalDipole = [0.0, 0.0, 0.0]
    #for i in range(len(indices)):
        # print "doing quats for molecule %d" % i
        # print "Dlen = %d " % len(DonatingTo[i])
        # print DonatingTo[i]
    for i in range(0, len(oxygenList)):
        oxygenA = oxygenList[i]
        print oxygenA.donors_
        print oxygenA.acceptors_
        
        #myPos = positions[i]

        #acceptor1 = DonatingTo[i][0]
        #acceptor2 = DonatingTo[i][1]
        acceptor1 = oxygenA.donors_[0]
        acceptor2 = oxygenA.donors_[1]

        
        '''
        if (acceptor1 == -1 and acceptor2 == -1):
            donor1 = AcceptingFrom[i][0]            
            donor2 = AcceptingFrom[i][1]
            npos = positions[donor1]
            apos1 = [npos[0] - myPos[0], npos[1] - myPos[1], npos[2] - myPos[2]]
            apos1 = wrapVector(apos1)
            npos = positions[donor2]
            apos2 = [npos[0] - myPos[0], npos[1] - myPos[1], npos[2] - myPos[2]]
            apos2 = wrapVector(apos2)
            tempX = [0.0, 0.0, 0.0]
            tempY = [0.0, 0.0, 0.0]
            tempZ = [0.0, 0.0, 0.0]
            tempZ[0] = 0.5*(apos1[0] + apos2[0])
            tempZ[1] = 0.5*(apos1[1] + apos2[1])
            tempZ[2] = 0.5*(apos1[2] + apos2[2])
            tempX[0] = (apos2[0] - apos1[0])
            tempX[1] = (apos2[1] - apos1[1])
            tempX[2] = (apos2[2] - apos1[2])
            tempZ = normalize(tempZ)
            tempX = normalize(tempX)
            tempY = cross(tempX, tempZ)
            dpos1[0] = tempZ[0] - 0.5*tempY[0]
            dpos1[1] = tempZ[1] - 0.5*tempY[1]
            dpos1[2] = tempZ[2] - 0.5*tempY[2]
            dpos2[0] = tempZ[0] + 0.5*tempY[0]
            dpos2[1] = tempZ[1] + 0.5*tempY[1]
            dpos2[2] = tempZ[2] + 0.5*tempY[2]
            
        else:
            if (acceptor1 == -1):
                tempVec = [0.0, 0.0, 0.0]
                for j in range(3):
                    thisNeighbor = neighbors[i][j]
                    npos = positions[thisNeighbor]
                    npos1 = [npos[0] - myPos[0], npos[1] - myPos[1], npos[2] - myPos[2]]
                    npos1 = wrapVector(npos1)
                    tempVec[0] = tempVec[0]  + npos1[0]
                    tempVec[1] = tempVec[1]  + npos1[1]
                    tempVec[2] = tempVec[2]  + npos1[2]                
                dpos1 = [-tempVec[0]/3.0, -tempVec[1]/3.0, -tempVec[2]/3.0]
                dpos1 = normalize(dpos1)
            else:
                a1pos = positions[acceptor1]
                dpos1 = [a1pos[0] - myPos[0], a1pos[1] - myPos[1], a1pos[2] - myPos[2]]
                dpos1 = wrapVector(dpos1)
                dpos1 = normalize(dpos1)
        

            if (acceptor2 == -1):
                tempVec = [0.0, 0.0, 0.0]
                for j in range(3):
                    thisNeighbor = neighbors[i][j]
                    npos = positions[thisNeighbor]
                    npos1 = [npos[0] - myPos[0], npos[1] - myPos[1], npos[2] - myPos[2]]
                    npos1 = wrapVector(npos1)
                    tempVec[0] = tempVec[0]  + npos1[0]
                    tempVec[1] = tempVec[1]  + npos1[1]
                    tempVec[2] = tempVec[2]  + npos1[2]                
                dpos2 = [-tempVec[0]/3.0, -tempVec[1]/3.0, -tempVec[2]/3.0]
                dpos2 = normalize(dpos2)
            else:
                a2pos = positions[acceptor2]
                dpos2 = [a2pos[0] - myPos[0], a2pos[1] - myPos[1], a2pos[2] - myPos[2]]
                dpos2 = wrapVector(dpos2)
                dpos2 = normalize(dpos2)        

        for j in range(3):
            uz[j] = (dpos1[j] + dpos2[j])/2.0
        uz = normalize(uz)
        for j in range(3):
            uy[j] = dpos2[j] - dpos1[j]
        uy = normalize(uy)
        ux = cross(uy, uz)
        ux = normalize(ux)

        q = [0.0, 0.0, 0.0, 0.0]

        # RotMat to Quat code is out of OpenMD's SquareMatrix3.hpp code:

        RotMat[0] = ux
        RotMat[1] = uy
        RotMat[2] = uz

        t = RotMat[0][0] + RotMat[1][1] + RotMat[2][2] + 1.0
        
        if( t > 1e-6 ):
            s = 0.5 / math.sqrt( t )
            q[0] = 0.25 / s
            q[1] = (RotMat[1][2] - RotMat[2][1]) * s
            q[2] = (RotMat[2][0] - RotMat[0][2]) * s
            q[3] = (RotMat[0][1] - RotMat[1][0]) * s
        else:
            ad1 = RotMat[0][0]
            ad2 = RotMat[1][1]
            ad3 = RotMat[2][2]

            if( ad1 >= ad2 and ad1 >= ad3 ):
                s = 0.5 / math.sqrt( 1.0 + RotMat[0][0] - RotMat[1][1] - RotMat[2][2] )
                q[0] = (RotMat[1][2] - RotMat[2][1]) * s
                q[1] = 0.25 / s
                q[2] = (RotMat[0][1] + RotMat[1][0]) * s
                q[3] = (RotMat[0][2] + RotMat[2][0]) * s
            elif ( ad2 >= ad1 and ad2 >= ad3 ):
                s = 0.5 / math.sqrt( 1.0 + RotMat[1][1] - RotMat[0][0] - RotMat[2][2] )
                q[0] = (RotMat[2][0] - RotMat[0][2] ) * s
                q[1] = (RotMat[0][1] + RotMat[1][0]) * s
                q[2] = 0.25 / s
                q[3] = (RotMat[1][2] + RotMat[2][1]) * s
            else:
                s = 0.5 / math.sqrt( 1.0 + RotMat[2][2] - RotMat[0][0] - RotMat[1][1] )
                q[0] = (RotMat[0][1] - RotMat[1][0]) * s
                q[1] = (RotMat[0][2] + RotMat[2][0]) * s
                q[2] = (RotMat[1][2] + RotMat[2][1]) * s
                q[3] = 0.25 / s

        quaternions[i] = q
        totalDipole = [totalDipole[0] + uz[0], totalDipole[1] + uz[1], totalDipole[2] + uz[2]]
    totalDipole = [-totalDipole[0], -totalDipole[1], -totalDipole[2]]
    print totalDipole
    Dipole = math.sqrt(dot(totalDipole, totalDipole))
    print 'Total Dipole Moment = %d' % Dipole
    if (Dipole > 40 or Dipole < -40):
        print "Bad Dipole, starting over"
        for i in range(len(indices)):
            del OldDonor[:]
            del AcceptingFrom[:]
            del DonatingTo[:]
            # del badChoices[:]
        randomizeProtons()
        computeQuats()
    else:
        print "All Done!"

    '''

        

def writeXYZfile(oxygenList, protonList,Hmat):
    outFile = open("test.xyz","w")

    numObjects = len(oxygenList) + len(protonList)
    outFile.write(str(numObjects) + "\n")
    outFile.write("\t%f%s\t%f\t%f\t%f%s\t%f\t%f\t%f%s\t%f\t%f\t%f\n" % 
                     (0.00000,";",
                      Hmat[0][0], Hmat[0][1], Hmat[0][2],";", 
                      Hmat[1][0], Hmat[1][1], Hmat[1][2],";", 
                      Hmat[2][0], Hmat[2][1], Hmat[2][2]))
    for i in range(0, len(oxygenList)):
        oxygenA = oxygenList[i]
        outFile.write("%s\t%f\t%f\t%f\n" % ('O', oxygenA.getPos()[0][0], oxygenA.getPos()[0][1], oxygenA.getPos()[0][2]))

    for i in range(0, len(protonList)):
        protonA = protonList[i]
        outFile.write("%s\t%f\t%f\t%f\n" % ('H', protonA.getPos()[0], protonA.getPos()[1], protonA.getPos()[2]))


    

def main(argv):
    parser = argparse.ArgumentParser(
        description='OpenMD module that generates proton disordered ice-Ih lattices given an (.xyz) coordinate file for the Oxygen lattice.',
        formatter_class=RawDescriptionHelpFormatter,
        epilog="Example: protonSampler -i oxygenLattice.xyz -o protonDisorderedIce.omd")
    parser.add_argument("-i","--inputFile=", action="store",
                            dest="inputFileName", help="use specified input file")
    parser.add_argument("-o","--outputfile=", action="store", dest="outputFileName",help="use specified output (.omd) file")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()

    if (not args.inputFileName):
        parser.error("No inputFile was specified")

    if (not args.outputFileName):
        parser.print_help()
        parser.error("No outputFile was specified")


    
    (oxygenList, Hmat) = readInputFile(args.inputFileName)
    findNeighbors(oxygenList)
    randomizeProtons(oxygenList)
    protonList = addProtonsToDonors(oxygenList)
    #computeQuats(oxygenList)
    
    writeXYZfile(oxygenList,protonList,Hmat)
    

if __name__ == "__main__":
    #if len(sys.argv) == 1:
        #parser.print_help()
        #usage() # need to change to call argeparse stuffs
        #sys.exit()
    main(sys.argv[1:])
