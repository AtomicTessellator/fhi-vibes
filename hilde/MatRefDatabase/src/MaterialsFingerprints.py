import sqlite3
import numpy as np
import yaml
# import matplotlib.pyplot as plt

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
        self.dE    = dE
        if( nbins == -1 and dE == -1.0 ):
            self.nbins = 32 if isB else 256
            self.fingerprint = ( self.constructBEFingerprint() if isElec else self.constructBPFingerprint() ) if isB else self.constructDFingerprint()
        elif nbins != -1:
            self.fingerprint = ( self.constructBEFingerprint() if isElec else self.constructBPFingerprint() ) if isB else self.constructDFingerprint()
        elif dE != -1.0:
            self.fingerprint = ( self.constructBEFingerprint() if isElec else self.constructBPFingerprint() ) if isB else self.constructDFingerprint()
        else:
            self.fingerprint = fp
            self.nbins = nbins

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

    def constructDFingerprint(self):
        dos = np.genfromtxt(self.spectraFiles[0])
        ener = self.getEner( )
        fingerprint = { "DOS" : np.zeros(shape=(len(ener)-1, 2)) }
        dEner = (ener[1]-ener[0])
        fingerprint["DOS"][:,0] = ener[:-1] + dEner/2.0
        for ff in dos:
            if(ff[0] >= self.minE and ff[0] <= self.maxE):
                fingerprint["DOS"][np.where(fingerprint["DOS"][:,0]-dEner/(2.0-1e-10) <= ff[0])[0][-1], 1] += ff[1]
        return fingerprint

    def constructBEFingerprint(self):
        bands = {}
        minE = self.minE
        maxE = self.maxE
        for pt in self.kpoints:
            qq = self.kpoints[pt]
            for sFile in self.spectraFiles:
                firstLine = list(filter(lambda x: x != '', open(sFile).readline().rstrip().split(" ")))
                lastLine  = list(filter(lambda x: x != '', open(sFile).readlines()[-1].rstrip().split(" ")))
                if( np.all(np.array(firstLine[1:4], dtype='float_') == qq )):
                    bands[pt] = np.array(firstLine[5::2], dtype='float_')
                    if( (self.minE ==  1.0e20) and (np.min(bands[pt]) < minE) ):
                        minE = np.min(bands[pt])
                    if( (self.maxE == -1.0e20) and (np.max(bands[pt]) > maxE) ):
                        maxE = np.max(bands[pt])
                elif( np.all(np.array(lastLine[1:4], dtype='float_') == qq )):
                    bands[pt] = np.array(lastLine[5::2], dtype='float_')
                    if( (self.minE ==  1.0e20) and (np.min(bands[pt]) < minE) ):
                        minE = np.min(bands[pt])
                    if( (self.maxE == -1.0e20) and (np.max(bands[pt]) > maxE) ):
                        maxE = np.max(bands[pt])
        if(self.minE ==  1.0e20):
            self.minE = minE
        if(self.maxE == -1.0e20):
            self.maxE = maxE
        return self.makeFingerPrintFromBand(bands)

    def constructBPFingerprint(self):
        bsSpect = yaml.load(open(self.spectraYAML, 'r'))
        bsLimited = []
        for bandpt in bsSpect["phonon"]:
            if("label" in bandpt):
                bsLimited.append(bandpt)
        bands = {}
        minE = self.minE
        maxE = self.maxE

        for pt in self.kpoints:
            qq = self.kpoints[pt]
            for bb in bsLimited:
                if( np.all(bb['q-position'] == qq) ):
                    band = []
                    for ff in bb['band']:
                        band.append(ff['frequency'])
                    bands[pt] = np.array(band)
            if( (self.minE ==  1e20) and np.min(bands[pt]) < minE):
                minE = np.min(bands[pt])
            if( (self.maxE == -1e20) and np.max(bands[pt]) > maxE):
                maxE = np.max(bands[pt])

        if(self.minE ==  1.0e20):
            self.minE = minE
        if(self.maxE == -1.0e20):
            self.maxE = maxE
        return self.makeFingerPrintFromBand(bands)

    def makeFingerPrintFromBand(self, bands):
        fingerprint = {}
        ener = self.getEner()
        for pt in bands:
            fingerprint[pt] = np.zeros(shape=(self.nbins,2))
            dEner = ener[1]-ener[0]
            fingerprint[pt][:,0] = ener[:-1] + (dEner)/2.0
            for ff in bands[pt]:
                if( self.minE <= ff and ff <= self.maxE ):
                    fingerprint[pt][np.where(fingerprint[pt][:,0]-dEner/(2.0-1e-10) <= ff)[0][-1], 1] += 1.0
        return fingerprint

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