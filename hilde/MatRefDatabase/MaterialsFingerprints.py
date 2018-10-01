import sqlite3
import numpy as np
import yaml
from phonopy import Phonopy

# Functions to define the energy bins
def get_ener( useFrequencies, frequencies, minE, maxE, nbins):
    """
    Get the energy bins used for making a fingerprint
    Args:
        useFrequencies: if True use the band/DOS frequencies given as the bin boundaries
        minE          : minimum energy mode to be included
        maxE          : maximum energy mode to be included
    nbins         : number of bins for the histogram
    Returns:
        energy bin labels (energy of the band or bin center point)
        energy bin boundaries
    """
    if not useFrequencies:
        enerBounds = np.linspace( minE, maxE, nbins+1 )
        return enerBounds[:-1] + (enerBounds[1]-enerBounds[0])/2.0, enerBounds
    else:
        return np.array(frequencies), np.append(frequencies, [frequencies[-1] + frequencies[-1]/10])

def find_min_E(bands):
    """
    Calculates the minimum energy mode in a band structure
    Args:
        bands: dictionary describing the phonon/electronic modes at all high symmetry points
    Returns:
        The global minimum mode energy energy
    """
    return np.min( np.array([ bands[pt] for pt in bands ]).flatten() )

def find_max_E(bands):
    """
    Calculates the maximum energy mode in a band structure
    Args:
        bands: dictionary describing the phonon/electronic modes at all high symmetry points
    Returns:
        The global maximum mode energy energy
    """
    return np.max( np.array([ bands[pt] for pt in bands ]).flatten() )

# Given a band structure or dos get the finger print
def get_fingerprint_bs(bands, minE, maxE, nbins):
    """
    Creates a dictionary of the band structure fingerprint at all high symmetry points, where the high symmetry point is the key
    Args:
        bands: A dict storing the phonon/electron mode energies at the high symmetry points, which are the keys for the dict
        minE : The minimum mode energy to include in the fingerprint
        maxE : The maximum mode energy to include in the fingerprint
        nbins: Number of bins to be used in the fingerprint
    Returns:
        fingerprint: a dictionary of the bs fingerprint at all high symmetry points
    """
    fingerprint = {}
    for pt in bands:
        ener, enerBounds = get_ener( len(bands[pt]) < nbins, bands[pt], minE, maxE, nbins)
        fingerprint[pt] = np.zeros(shape=(len(ener), 2))
        fingerprint[pt][:,0] = ener
        fingerprint[pt][:,1] = np.histogram(bands[pt], enerBounds)[0]
    return fingerprint

def get_fingerprint_dos(dos, minE, maxE, nbins):
    """
    Creates a dictionary of the density of states fingerprint
    Args:
        dos  : an numpy array of the density of states
        minE : The minimum mode energy to include in the fingerprint
        maxE : The maximum mode energy to include in the fingerprint
        nbins: Number of bins to be used in the fingerprint
    Returns:
        fingerprint: a dictionary of the dos fingerprint where the only key is "DOS"
    """
    fingerprint = {}
    if(dos.shape[0] < nbins):
        fingerprint["DOS"] = dos[np.where((dos[:,0] >= minE) & (dos[:,0] <= maxE))]

    else:
        ener, enerBounds = get_ener(False, dos[:,0], minE, maxE, nbins)
        fingerprint = { "DOS" : np.zeros(shape=(len(ener), 2) ) }
        fingerprint["DOS"][:,0] = ener
        fingerprint["DOS"][:,1] = np.histogram(dos, enerBounds)[0]
    # print(fingerprint["DOS"])
    return fingerprint

# Function to calculate the modes at the high symmetry points
def get_elec_bands(spectraFiles, k_points):
    """
    Generates a dict describing the electronic band structure at from a list of files
    Args:
        spectraFiles: A list of files with the bands are defined
        k_points    : a list of high symmetry points
    Returns:
        bands: a dict describing the electronic band structure with the keys being the high symmetry pints
    """
    bands = {}
    for pt in k_points:
        for sFile in spectraFiles:
            firstLine = list(filter(lambda x: x != '', open(sFile).readline().rstrip().split(" ")))
            lastLine  = list(filter(lambda x: x != '', open(sFile).readlines()[-1].rstrip().split(" ")))
            if( np.all(np.array(firstLine[1:4], dtype='float_') == k_points[pt] )):
                bands[pt] = np.array(firstLine[5::2], dtype='float_')
            elif( np.all(np.array(lastLine[1:4], dtype='float_') == k_points[pt] )):
                bands[pt] = np.array(lastLine[5::2], dtype='float_')
    return bands

def get_phonon_bands_phonopy(phonon, q_points):
    """
    Generates a dict describing the phonon band structure at from a phonopy object
    Args:
        phonon  : The phonopy object with the band structure calculated
        q_points: a list of high symmetry points
    Returns:
        bands: a dict describing the phonon band structure with the keys being the high symmetry pints
    """
    bands = {}
    for pt in q_points:
        bands[pt] = phonon.get_frequencies(q_points[pt])
    return bands

def get_phonon_bands_yaml(yamlFile, q_points):
    """
    Generates a dict describing the phonon band structure at from a yaml file
    Args:
        yamlFile: The phonopy generated yaml file describing the band structure
        q_points: a list of high symmetry points
    Returns:
        bands: a dict describing the phonon band structure with the keys being the high symmetry pints
    """
    bands = {}
    bsSpect = yaml.load(open(yamlFile, 'r'))
    bsLimited = []
    for bandpt in bsSpect["phonon"]:
        if("label" in bandpt):
            bsLimited.append(bandpt)
    bands = {}
    for pt in q_points:
        for bb in bsLimited:
            if( np.all(bb['q-position'] == q_points[pt]) ):
                bands[pt] = [ff['frequency'] for ff in bb['band']]
    return bands

# Functions to get the fingerprint from various input values
def get_phonon_bs_fingerprint_phononpy(phonon, q_points, minE=None, maxE=None, nbins=32):
    """
    Generates the phonon band structure fingerprint for a bands structure stored in a phonopy object
    Args:
        phonon  : The phonopy generated yaml file describing the band structure
        q_points: a list of high symmetry points
        minE    : The minimum mode energy to include in the fingerprint
        maxE    : The maximum mode energy to include in the fingerprint
        nbins   : Number of bins to be used in the fingerprint
    Returns:
        The phonon band structure fingerprint
    """
    bands = get_phonon_bands_phonopy(phonon, q_points)
    return get_fingerprint_bs(bands, find_min_E(bands) if minE is None else minE, find_max_E(bands) if maxE is None else maxE, nbins)

def get_phonon_bs_fingerprint_yaml(yamlFile, q_points, minE=None, maxE=None, nbins=32):
    """
    Generates the phonon band structure fingerprint for a bands structure stored in a phonopy generated yaml file
    Args:
        yankFile: The phonopy generated yaml file describing the band structure
        q_points: a list of high symmetry points
        minE    : The minimum mode energy to include in the fingerprint
        maxE    : The maximum mode energy to include in the fingerprint
        nbins   : Number of bins to be used in the fingerprint
    Returns:
        The phonon band structure fingerprint
    """
    bands = get_phonon_bands_yaml(yamlFile, bands)
    return get_fingerprint_bs(bands, find_min_E(bands) if minE is None else minE, find_max_E(bands) if maxE is None else maxE, nbins)

def get_elec_bs_fingerprint(spectraFiles, k_points, minE=None, maxE=None, nbins=32):
    """
    Generates the electronic band structure fingerprint for a band structure stored in a set of files
    Args:
        spectraFiles: A list of files with the bands are defined
        k_points    : a list of high symmetry points
        minE        : The minimum mode energy to include in the fingerprint
        maxE        : The maximum mode energy to include in the fingerprint
        nbins       : Number of bins to be used in the fingerprint
    Returns:
        The electronic band structure fingerprint
    """
    bands = get_elec_bands(spectraFiles, k_points)
    return get_fingerprint_bs(bands, find_min_E(bands) if minE is None else minE, find_max_E(bands) if maxE is None else maxE, nbins)

def get_dos_fingerprint(dosFile, minE=None, maxE=None, nbins=256):
    """
    Generates the density of states fingerprint from a file describing the density of states
    Args:
        dosFile: A file storing the density of states in 2 columns: (energy, electron density)
        minE   : The minimum mode energy to include in the fingerprint
        maxE   : The maximum mode energy to include in the fingerprint
        nbins  : Number of bins to be used in the fingerprint
    Returns:
        The density of states fingerprint
    """
    dos = np.genfromtxt(dosFile)
    return get_fingerprint_dos(dos, np.min(dos[:,0]) if minE is None else minE, np.max(dos[:,0]) if maxE is None else maxE, nbins)

def get_phonon_dos_fingerprint_phononpy(phonon, minE=None, maxE=None, nbins=256):
    dos = np.array(phonon.get_total_DOS()).transpose()
    return get_fingerprint_dos(dos, np.min(dos[:,0]) if minE is None else minE, np.max(dos[:,0]) if maxE is None else maxE, nbins)


# Class definitions for incorporation into databases
class MaterialsFingerprint(object):
    kpoints = {}
    spectraFiles = []
    spectraYAML = ""
    isB = False
    isElec = False
    isDStream = False

    def __init__(self, isElec, isB, nbins=None, dE=None, minE=None, maxE=None, fp = {}):
        self.isB = isB
        self.isElec = isElec
        self.minE = minE
        self.maxE = maxE
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

    def scalar_product(self, otherFP):
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
    def __init__(self, isElec, isB, nbins=None, dE=None, minE=None, maxE=None, fp={}, kpoints={}, spectraFiles=[]):
        self.kpoints = kpoints
        self.spectraFiles = spectraFiles
        self.isElec = isElec
        # Determine Energy range
        dos = np.genfromtxt(self.spectraFiles[0])
        self.minE = np.min( dos[:,0] ) if minE is None else minE
        self.maxE = np.max( dos[:,0] ) if maxE is None else maxE
        self.nbins = nbins
        if self.nbins is None:
            self.nbins = 256 if dE is None else (self.maxE-self.minE)/dE
        self.dE = dE
        if self.dE is None:
            self.dE = (self.maxE - self.minE)/(256.0) if dE is None else (self.maxE-self.minE)/nbins
        #make the fingerprint
        if(fp == {}):
            self.fingerprint = getDOSFingerprint(dos, self.minE, self.maxE, self.nbins)
        else:
            self.fingerprint = fp

class BandStructureFingerprint(MaterialsFingerprint):
    def __init__(self, isElec, isB, nbins=None, dE=None, minE=None, maxE=None, fp={}, kpoints={}, spectraFiles=[], spectraYAML="", phonon=None):
        start = clock()
        self.kpoints = kpoints
        self.spectraFiles = spectraFiles
        self.spectraYAML = spectraYAML
        self.isB = isB
        self.isElec = isElec
        bands = {}
        # Take in band structure data
        if(self.isElec):
            bands = get_elec_bands(self.spectraFiles, self.kpoints)
        else:
            if(phonon == None):
                bands = get_phonon_bands_yaml(self.yamlFile, self.kpoints)
            else:
                bands = get_phonon_bands_phonon(self.yamlFile, self.kpoints)

        # Find energy bins
        self.minE = find_min_E(bands) if minE is None else minE
        self.maxE = find_max_E(bands) if maxE is None else maxE
        self.nbins = nbins
        if self.nbins is None:
            self.nbins = 32 if dE is None else (self.maxE-self.minE)/dE
        self.dE = dE
        if self.dE is None:
            self.dE = (self.maxE - self.minE)/(256.0) if dE is None else (self.maxE-self.minE)/nbins

        # Make the fingerprint
        if( fp == {} ):
            self.fingerprint = get_fingerprint_bs(bands, self.minE, self.maxE, self.nbins)
        else:
            self.fingerprint = fp
