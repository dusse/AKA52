class Initializer:

    
    def __init__(self):
        #   0: run from scratch
        #   1: restart, set input file
        self.runType = 0
        self.inputFile = "./single_las/single_laser3D_restart.h5"
        
        # BOX
        self.dim=3
        
        self.boxSize    = [16.0, 16.0, 4.0] # in d0 - ion inertia length
        self.boxSizePxl = [80 , 80 , 20]

        self.bcType = [1,1,0] #1 - periodic 0 - ideal
        
        self.partclBcType = [1,1,2] #2 - reflect 1 - periodic 0 - outflow
        
        self.mpiCores  = [2,1,1]
        
        # time
        self.ts = 0.001
        self.maxtsnum = 80001
        self.outputStride = 50
        
        # output. need to create it before
        self.outputDir = "single_las/"
        self.fileTemplate = "single_laser3D_"
    
        # Particles
        self.ppc = 10
        self.ppcMinDens = 0.1
        
        self.numOfSpecies = 2
        self.masses  = [1, 1]
        self.charges = [1, 1]
        self.dens = [0.4, 2.0]
        self.vel1 = 0.022
        self.vel2 = 0.01
       
        self.Pele0 = 0.0001 # constant pressure everywhere 
        self.surface_pos = 1.0
        
        # Laser spots
        self.spotsNum = 1
        self.spotPos  = [8.0,8.0,1.0]
        
        self.radius = 6.0

        self.type2Heat = 2
        self.temp2load = 0.0001
        self.dens2sustain = 2.0
        self.temp2sustain = 10.0
    
        #magnetic field magnitude
        self.Bfield = [1.0, 0.0, 0.0]
        
        #ohm's law
        self.hyperviscosity = 0.001

        # pressure tensor
        self.emass = 0.1
        self.tau   = 0.01 # ~ time
        self.relaxFactor = 1/self.tau
        self.smoothStride = 4

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
        return self.Pele0

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
        checkVar = 0
        if(self.dim == 1):
            checkVar = x
        if(self.dim == 2):
            checkVar = y
        if(self.dim == 3):
            checkVar = z
        
        if(checkVar>=self.surface_pos):
            return self.dens[0]
        else:
            return self.ZERO


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
        checkVar = 0
        if(self.dim == 1):
            checkVar = x
        if(self.dim == 2):
            checkVar = y
        if(self.dim == 3):
            checkVar = z
        
        if(checkVar<self.surface_pos):
            return self.dens[1]
        else:
            return self.ZERO

#           species 2: modulus of velocity for Maxwell distribution
    def getVelocityX4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityY4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityZ4species2(self, x, y, z):
        return self.vel2
    
    
    #   physics: laser
    def getNumberOfLaserSpots(self):
        return self.spotsNum
    
    # type defined above
    def getParticleType2Heat(self):
        return self.type2Heat
    
    def getParticleTemp2Load(self):
        return self.temp2load
    
    def getParticleDensity2sustain(self):
        return self.dens2sustain
    
    def getParticleTemperature2sustain(self):
        return self.temp2sustain
    
    def getLaserSpotCenterXposition1(self):
        return self.spotPos[0]
    
    def getLaserSpotCenterYposition1(self):
        return self.spotPos[1]
    
    def getLaserSpotCenterZposition1(self):
        return self.spotPos[2]
    
    def getLaserSpotRadius1(self):
        return self.radius
    
    def getLaserSpotCenterXposition2(self):
        return self.spot2Pos[0]
    
    def getLaserSpotCenterYposition2(self):
        return self.spot2Pos[1]
    
    def getLaserSpotCenterZposition2(self):
        return self.spot2Pos[2]
    
    def getLaserSpotRadius2(self):
        return self.radius
    
    
    #   physics: magnetic field
    def getBfieldX(self, x, y, z):
        return self.Bfield[0]
    
    def getBfieldY(self, x, y, z):
        return self.Bfield[1]
    
    def getBfieldZ(self, x, y, z):
        return self.Bfield[2]
    
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











