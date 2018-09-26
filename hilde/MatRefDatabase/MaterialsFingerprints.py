import sqlite3
import numpy as np
import yaml
from time import clock
# import matplotlib.pyplot as plt

# This class will be split into a band/DOS fingerprint child classes with a unified input structure
# Right now it is done this way since we need to decide the format of lower level databases to see
# what is available at which level
class MaterialsFingerprint(object):
    kpoints = {}
    spectraFiles = []
    spectraYAML = ""
    isB = False
    isElec = False
    isDStream = False

    def __init__(self, isElec, isB, nbins=-1, dE=-1.0, fp={},kpoints={}, spectraFiles=[], spectraYAML="", minE=None, maxE=None):
        self.kpoints = kpoints
        self.spectraFiles = spectraFiles
        self.spectraYAML = spectraYAML
        self.isB = isB
        self.isElec = isElec
        self.minE = minE if(minE !=None) else  1e20 if(isB) else np.min( np.genfromtxt(self.spectraFiles[0])[:,0] )
        self.maxE = maxE if(maxE !=None) else -1e20 if(isB) else np.max( np.genfromtxt(self.spectraFiles[0])[:,0] )
        self.nbins = nbins
        self.nbins = nbins
        self.dE    = dE
        self.fingerprint = fp

    def __conform__(self, protocol):
        if protocol is sqlite3.PrepareProtocol:
            frmt = "%i;%r;%r" % (self.nbins, isB, isElec)
            for pt in self.fingerprint:
                frmt += ";%s" % (pt)
                for ff in self.fingerprint[pt]:
                    frmt += ";%f;%f" % (ff[0], ff[1])
            return frmt
        return ""

    def getEner(self):
        if(self.nbins != -1):
            return np.linspace( self.minE, self.maxE, self.nbins+1 )
        else:
            return np.arange( self.minE, self.maxE+self.dE/10.0, self.dE )


    def scalarProduct(self, otherFP):
        if(len(self.fingerprint) != len(otherFP.fingerprint) ):
            print("The two finger prints do not represent the same q points, the scalar product is zero")
            return 0.0
        scalar = 0.0
        for key in self.fingerprint:
            if key not in otherFP.fingerprint:
                print("The two finger prints do not represent the same q points, the scalar product is zero")
                return 0.0
            elif(self.isB and np.sum(self.fingerprint[key][:,1]) != np.sum(otherFP.fingerprint[key][:,1]) ):
                print("The two fingerprints should have the same number bands, the scalar product is 0")
                return 0.0
            elif( np.any(self.fingerprint[key][:,0] != otherFP.fingerprint[key][:,0] ) ):
                print("The energy bins need to be exactly the same, the scalar product is 0.0")
                return 0.0
            scalar += np.dot(self.fingerprint[key][:,1], otherFP.fingerprint[key][:,1])/np.linalg.norm(self.fingerprint[key][:,1])**2
        return scalar/len(self.fingerprint)

class DOSFingerprint(MaterialsFingerprint):
    def __init__(self, isElec, isB, nbins=None, dE=None, fp={},kpoints={}, spectraFiles=[], spectraYAML="", minE=None, maxE=None):
        self.kpoints = kpoints
        self.spectraFiles = spectraFiles
        self.spectraYAML = spectraYAML
        self.isElec = isElec
        self.minE = np.min( np.genfromtxt(self.spectraFiles[0])[:,0] ) if minE is None else minE
        self.maxE = np.max( np.genfromtxt(self.spectraFiles[0])[:,0] ) if maxE is None else maxE
        if nbins is None:
            self.nbins = 256 if dE is None else (self.maxE-self.minE)/dE
        else:
            self.nbins = nbins
        if dE is None:
            self.dE = (self.maxE - self.minE)/(256.0) if dE is None else (self.maxE-self.minE)/nbins
        else:
            self.dE = dE
        self.ener = self.getEner()

        if(fp == {}):
            dos = np.genfromtxt(self.spectraFiles[0])
            self.fingerprint = { "DOS" : np.zeros(shape=(len(self.ener)-1, 2)) }
            self.fingerprint["DOS"][:,0] = self.ener[:-1] + self.dE/2.0
            for ff in dos:
                if(ff[0] >= self.minE and ff[0] <= self.maxE):
                    self.fingerprint["DOS"][np.where(self.fingerprint["DOS"][:,0]-self.dE/(2.0-1e-10) <= ff[0])[0][-1], 1] += ff[1]
        else:
            self.fingerprint = fp

class BandStructureFingerprint(MaterialsFingerprint):
    def __init__(self, isElec, isB, nbins=None, dE=None, fp={},kpoints={}, spectraFiles=[], spectraYAML="", minE=None, maxE=None):
        start = clock()
        self.kpoints = kpoints
        self.spectraFiles = spectraFiles
        self.spectraYAML = spectraYAML
        self.isB = isB
        self.isElec = isElec
        self.bands = {}
        if(self.isElec):
            for pt in self.kpoints:
                for sFile in self.spectraFiles:
                    firstLine = list(filter(lambda x: x != '', open(sFile).readline().rstrip().split(" ")))
                    lastLine  = list(filter(lambda x: x != '', open(sFile).readlines()[-1].rstrip().split(" ")))
                    if( np.all(np.array(firstLine[1:4], dtype='float_') == self.kpoints[pt] )):
                        self.bands[pt] = np.array(firstLine[5::2], dtype='float_')
                    elif( np.all(np.array(lastLine[1:4], dtype='float_') == self.kpoints[pt] )):
                        self.bands[pt] = np.array(lastLine[5::2], dtype='float_')
        else:
            bsSpect = yaml.load(open(self.spectraYAML, 'r'))
            bsLimited = []
            for bandpt in bsSpect["phonon"]:
                if("label" in bandpt):
                    bsLimited.append(bandpt)
            for pt in self.kpoints:
                for bb in bsLimited:
                    if( np.all(bb['q-position'] == self.kpoints[pt]) ):
                        self.bands[pt] = [ff['frequency'] for ff in bb['band']]
        self.minE = self.findMinE() if minE is None else minE
        self.maxE = self.findMaxE() if maxE is None else maxE
        if nbins is None:
            self.nbins = 256 if dE is None else (self.maxE-self.minE)/dE
        else:
            self.nbins = nbins
        if dE is None:
            self.dE = (self.maxE - self.minE)/(256.0) if dE is None else (self.maxE-self.minE)/nbins
        else:
            self.dE = dE
        self.ener = self.getEner()

        if( fp == {} ):
            self.fingerprint = {}
            for pt in self.bands:
                self.fingerprint[pt] = np.zeros(shape=(self.nbins,2))
                self.fingerprint[pt][:,0] = self.ener[:-1] + (self.dE)/2.0
                for ff in self.bands[pt]:
                    if( self.minE <= ff and ff <= self.maxE ):
                        self.fingerprint[pt][np.where(self.fingerprint[pt][:,0]-self.dE/(2.0-1e-10) <= ff)[0][-1], 1] += 1.0
        else:
            self.fingerprint = fp

    def findMinE(self):
        minE = 1e20
        for pt in self.bands:
            if(np.min(self.bands[pt]) < minE):
                minE = np.min(self.bands[pt])
        return minE

    def findMaxE(self):
        maxE = -1e20
        for pt in self.bands:
            if(np.max(self.bands[pt]) > maxE):
                maxE = np.max(self.bands[pt])
        return maxE