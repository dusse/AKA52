import numpy as np

class Initializer:

    
    def __init__(self):
        #   0: run from scratch
        #   1: restart, set input file
        self.runType = 0
        self.inputFile = "./harris_out/har_restart.h5"
        
        
        # BOX
        self.boxSize    = [25.6, 12.8, 0.064] # in d0 - ion inertia length
        self.boxSizePxl = [400 , 200 , 1]

        self.bcType = [1,0,1] #1 - periodic 0 - ideal
        
        self.partclBcType = [1,2,1] #2 - reflect 1 - periodic 0 - outflow
        
        self.mpiCores  = [1,2,1]
        
        # time
        self.ts = 0.0002
        self.maxtsnum = 100
        self.outputStride = 50
        
        self.middle = [0.5*self.boxSize[0], 0.5*self.boxSize[1]]
        
        # output. need to create it before
        self.outputDir = "harris_out/"
        self.fileTemplate = "harris_"
    
        # Particles
        self.ppc = 20
        self.ppcMinDens = 0.2
        
        self.numOfSpecies = 2
        self.masses  = [1, 1]
        self.charges = [1, 1]
        self.dens = [0.2, 1.0]
        TeTi = 0.2
        pmag = 0.5
        Tion = pmag / (1 + TeTi)
        Tele = pmag * TeTi / (1 + TeTi)
        
        print 'Tion =  ', Tion
        print 'Tele =  ', Tele
    
        self.vel1 = np.sqrt(Tion)
        self.vel2 = np.sqrt(Tion)
        
        self.sheetHalfWidth = 0.5
        
        self.tempEle = Tele
        
        #magnetic field magnitude
        self.Bfield = [0.0, 0.0, 0.0]
        
        #ohm's law, small scale dissipation
        self.hyperviscosity = 0.001

        # pressure tensor
        self.emass = 0.04
        self.tau   = 1.0 # ~ time
        self.relaxFactor = 1/self.tau
        self.smoothStride = 2*self.maxtsnum # no smoothing

        #misc
        self.ZERO = 0.0
    
    
    #   0: run from scratch
    #   1: restart, set input file
    def getRunType(self):
        return self.runType
    
    def getInputFile(self):
        return self.inputFile
    
    #   spatial: left lower corner is (0,0,0)
    #   spatial: total box length in normalized units
    def getXright(self):
        return self.boxSize[0]
    
    def getYright(self):
        return self.boxSize[1]
    
    def getZright(self):
        return self.boxSize[2]
    
    
    # if number of pixels is less than 2, there is no direvative
    #   total box length in pixels Lx
    def getXresolution(self):
        return self.boxSizePxl[0]
    #   total box length in pixels Ly
    def getYresolution(self):
        return self.boxSizePxl[1]
    #   total box length in pixels Lz
    def getZresolution(self):
        return self.boxSizePxl[2]
    
    def getXmpiDomainNum(self):
        return self.mpiCores[0]
    
    def getYmpiDomainNum(self):
        return self.mpiCores[1]
    
    def getZmpiDomainNum(self):
        return self.mpiCores[2]
    
    
    #   BC
    def getFieldBCTypeX(self):
        return self.bcType[0]
    
    def getFieldBCTypeY(self):
        return self.bcType[1]
    
    def getFieldBCTypeZ(self):
        return self.bcType[2]
    
    def getParticleBCTypeX(self):
        return self.partclBcType[0]
    
    def getParticleBCTypeY(self):
        return self.partclBcType[1]
    
    def getParticleBCTypeZ(self):
        return self.partclBcType[2]
    
    
    #   time
    def getTimestep(self):
        return self.ts
    
    def getMaxTimestepsNum(self):
        return self.maxtsnum
    
    #   output
    def getOutputDir(self):
        return self.outputDir
    
    def getOutputFilenameTemplate(self):
        return self.fileTemplate
    
    def getOutputTimestep(self):
        return self.outputStride



    def getElectronPressure(self, x, y, z):
        y0 = (y-self.middle[1])/self.sheetHalfWidth
        coshSquare = np.cosh(y0)*np.cosh(y0)
        mainpop = self.tempEle*self.dens[1]/coshSquare
        bgpop = self.tempEle*self.dens[0]
        return mainpop+bgpop
    

    #   physics: particles
#   first set number of used species
    def getNumOfSpecies(self):
        return self.numOfSpecies
    
    def getParticlesPerCellNumber(self):
        return self.ppc
    
    
    def getMinimumDens2ResolvePPC(self):
        return self.ppcMinDens
    
    #           species 1
    def getMass4species1(self):
        return self.masses[0]
    
    def getCharge4species1(self):
        return self.charges[0]
    
    def getDensity4species1(self, x, y, z):
        return self.dens[0]


#           species 1: modulus of velocity for Maxwell distribution
    def getVelocityX4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityY4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityZ4species1(self, x, y, z):
        return self.vel1
    
    #           species 2
    def getMass4species2(self):
        return self.masses[1]
    
    def getCharge4species2(self):
        return self.charges[1]
    
    
    def getDensity4species2(self, x, y, z):
        y0 = (y-self.middle[1])/self.sheetHalfWidth
        return 1./(np.cosh(y0)*np.cosh(y0))

#           species 2: modulus of velocity for Maxwell distribution
    def getVelocityX4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityY4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityZ4species2(self, x, y, z):
        return self.vel2
    
    

    
    #   physics: magnetic field
    def getBfieldX(self, x, y, z):
        sigma = 2.0
        weight = 0.1
        x0 = (x-self.middle[0])/self.sheetHalfWidth
        y0 = (y-self.middle[1])/self.sheetHalfWidth
        pinch = - weight * y0 * np.exp(-(x0*x0 + y0*y0) / (sigma*sigma));
        
        return np.tanh(y0)+pinch
    
    def getBfieldY(self, x, y, z):
        sigma = 2.0
        weight = 0.1
        x0 = (x-self.middle[0])/self.sheetHalfWidth
        y0 = (y-self.middle[1])/self.sheetHalfWidth
        pinch = weight * x0 * np.exp(-(x0*x0 + y0*y0) / (sigma*sigma));
        
        return pinch
    
    
    def getBfieldZ(self, x, y, z):
        return 0.0
    
    #   physics: ohm's law hyper-viscosity
    def getHyperviscosity(self):
        return self.hyperviscosity
    
    
    #   physics: pressure evolution
    def getElectronMass(self):
        return self.emass

    #   physics: pressure evolution : isotropization
    def getRelaxFactor(self):
        return self.relaxFactor
    
    #   physics: pressure evolution : smoothing
    def getElectronPressureSmoothingStride(self):
        return self.smoothStride

#   physics: laser
    def getNumberOfLaserSpots(self):
        return 0











