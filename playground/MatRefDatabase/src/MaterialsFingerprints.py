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

    def __init__(self, isElec, isB, nbins=0, fp={},kpoints={}, spectraFiles=[], spectraYAML="", minE=1e20, maxE=-1e20):
        self.kpoints = kpoints
        self.spectraFiles = spectraFiles
        self.spectraYAML = spectraYAML
        self.isB = isB
        self.isElec = isElec
        self.minE = minE
        self.maxE = maxE
        if( nbins == 0 ):
            self.nbins = 32 if isB else 256
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

    def constructDFingerprint(self):
        fingerprint = { "DOS" : np.zeros(shape=(self.nbins, 2)) }
        dos = np.genfromtxt(self.spectraFiles[0])
        ener = np.linspace( np.min( dos[:,0] ), np.max( dos[:,0] ), self.nbins+1 )
        dEner = (ener[1]-ener[0])
        fingerprint["DOS"][:,0] = ener[:-1] + dEner/2.0
        for ff in dos:
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
        ener = np.linspace(self.minE, self.maxE, self.nbins+1)
        for pt in bands:
            fingerprint[pt] = np.zeros(shape=(self.nbins,2))
            dEner = ener[1]-ener[0]
            fingerprint[pt][:,0] = ener[:-1] + (dEner)/2.0
            for ff in bands[pt]:
                if( self.minE <= ff and ff <= self.maxE ):
                    fingerprint[pt][np.where(fingerprint[pt][:,0]-dEner/(2.0-1e-10) <= ff)[0][-1], 1] += 1.0
        return fingerprint

def convert_fingerprint(s):
    sSplit = s.split(b";")
    fingerprint = {}
    nbins = int(sSplit[0])
    isB = bool(sSplit[1])
    isE = bool(sSplit[2])
    ii = 3;
    while ii < len(sSplit):
        pt = sSplit[ii].decode("ascii")
        ii += 1
        fingerprint[pt] = np.ndarray(shape=(nbins,2), dtype=float)
        for jj in range(nbins):
            fingerprint[pt][jj,:] = list(map(float,sSplit[ii:ii+2]))
            ii += 2
    return MaterialsFingerprint(isE, isB, nbins=nbins, fp=fingerprint)

def adapt_fingerprint(fp):
    frmt = ("%i;%r;%r" % (fp.nbins, fp.isB, fp.isElec))
    for pt in fp.fingerprint:
        frmt += (";%s" % (pt))
        for ff in fp.fingerprint[pt]:
            frmt += ";%f;%f" % (ff[0], ff[1])
    return frmt.encode()
