
class Initializer:

    
    def __init__(self):
        #   0: run from scratch
        #   1: restart, set input file
        self.runType = 0
        self.inputFile = "./output/_restart.h5"
        
        # BOX
        self.dim=3
        
        self.boxSize    = [30.0, 30.0, 30.0] # in d0 - ion inertia length
        self.boxSizePxl = [30, 30, 30]

        self.pxSize = self.boxSize[2]/self.boxSizePxl[2]
        
        self.bcType = [0,0,0] #1 - periodic 0 - ideal
        
        self.partclBcType = [0,0,0] #2 - reflect 1 - periodic 0 - outflow
        
        self.dampingBoundaryWidth = [[5,5], [5,5], [5,5]]
         
        self.mpiCores  = [1,2,1]
        
        # time
        self.ts = 0.01
        self.maxtsnum = 501
        self.outputStride = 50

        # output. need to create it before
        self.outputDir = "output/"
        self.fileTemplate = "free_3D_"
    
        # Particles
        self.ppc4load = 20
        self.ppc = [1, 1]
        self.ppcMinDens = 0.005

        self.numOfSpecies = 2
        self.masses  = [1, 1]
        self.charges = [1, 1]
        self.dens = [self.ppcMinDens, 2.0]
        self.vel1 = 0.0
        self.vel2 = 0.0
       
        self.Pele0 = 0.0000001
        self.surface_pos = 8.0
        self.trgtwidth = 3.0
        # Laser spots
        self.spotsNum = 1
        self.spotPos = [0.5*self.boxSize[0], 0.5*self.boxSize[1], self.surface_pos]
                        
       
        self.radius1 = 10.0
        self.radius2 = 10.0

        self.type2Load = 2
        self.thermVion = 0.5

        self.dens2sustain = self.dens[1]
        self.temp2sustain = 5.0
    
        #magnetic field magnitude
        self.Bfield = [0.0, 0.0, 0.0]
        
        #ohm's law
        self.resistivity = 0.0

        # pressure tensor
        self.emass = 0.1
        self.tau   = 0.01 # ~ time
        self.relaxFactor = 1/self.tau
        self.smoothStride = 1.5*self.maxtsnum

        #misc
        self.ZERO = 0.0
        self.PI = 3.14159265358979323846
        self.FROZEN = 1
    
    
    
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

    def getDampingBoundaryWidthXleft(self):
        return self.dampingBoundaryWidth[0][0]

    def getDampingBoundaryWidthXright(self):
        return self.dampingBoundaryWidth[0][1]

    def getDampingBoundaryWidthYleft(self):
        return self.dampingBoundaryWidth[1][0]

    def getDampingBoundaryWidthYright(self):
        return self.dampingBoundaryWidth[1][1]

    def getDampingBoundaryWidthZleft(self):
        return self.dampingBoundaryWidth[2][0]

    def getDampingBoundaryWidthZright(self):
        return self.dampingBoundaryWidth[2][1]
    
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
    
    def getMinimumDens2ResolvePPC(self):
        return self.ppcMinDens
    
    #           species 1
    def getPPC4species1(self):
        return self.ppc[0]


    def getMass4species1(self):
        return self.masses[0]
    
    def getCharge4species1(self):
        return self.charges[0]

    def getIfParticleTypeIsFrozen4species1(self):
        return self.FROZEN
    
    def getDensity4species1(self, x, y, z):
        Ly2 = 0.5*self.boxSize[1]
        Lz = self.boxSize[2]
        if (z <= (self.surface_pos-self.trgtwidth) or z >= self.surface_pos):
            return self.dens[0]
        else:
            return self.ZERO

    def getIfParticleTypeIsFrozen4species1(self):
        return self.FROZEN

    # species 1: modulus of velocity for Maxwell distribution
    def getVelocityX4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityY4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityZ4species1(self, x, y, z):
        return self.vel1


    # species 1: modulus of fluid velocity
    def getFluidVelocityX4species1(self, x, y, z):
        return 0.0
    
    def getFluidVelocityY4species1(self, x, y, z):
        return 0.0
    
    def getFluidVelocityZ4species1(self, x, y, z):
        return 0.0
    
    #           species 2
    def getPPC4species2(self):
        return self.ppc[1]
    
    def getMass4species2(self):
        return self.masses[1]
    
    def getCharge4species2(self):
        return self.charges[1]
        
    def getDensity4species2(self, x, y, z):
        Ly2 = 0.5*self.boxSize[1]
        Lz = self.boxSize[2]
        if (z < self.surface_pos and z > (self.surface_pos - self.trgtwidth) ):
            return self.dens[1]
        else:
            return self.ZERO

    def getIfParticleTypeIsFrozen4species2(self):
        return self.FROZEN


    # species 2: modulus of velocity for Maxwell distribution
    def getVelocityX4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityY4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityZ4species2(self, x, y, z):
        return self.vel2
    
    # species 2: modulus of fluid velocity
    def getFluidVelocityX4species2(self, x, y, z):
        return 0.0
    
    def getFluidVelocityY4species2(self, x, y, z):
        return 0.0
    
    def getFluidVelocityZ4species2(self, x, y, z):
        return 0.0

    #   physics: laser
    def getNumberOfLaserSpots(self):
        return self.spotsNum
    
    #   physics: laser pulse duration in timesteps
    def getLaserPulseDuration(self):
        return self.maxtsnum
    
    # particles type (1,2,...) is needed to copy mass and charge for injected particles
    def getParticleType2Load(self):
        return self.type2Load

    # can specify ppc for the injected fraction
    def getPPC4loadedParticles(self):
        return self.ppc4load
    
    # injected particles: modulus of fluid velocity
    def getFluidVelocityX4InjectedParticles(self, x, y, z):
        return 0.0
    
    def getFluidVelocityY4InjectedParticles(self, x, y, z):
        return 0.0
    
    def getFluidVelocityZ4InjectedParticles(self, x, y, z):
        return 0.0

    # injected particles: modulus of thermal velocity
    def getVelocityX4InjectedParticles(self, x, y, z):
        Ly2 = 0.5*self.boxSize[1]
        Lz = self.boxSize[2]
        if (z > self.surface_pos):
            return self.ZERO
        
        tx = x - self.spotPos[0]
        ty = y - self.spotPos[1]
        tz = z - self.spotPos[2]+self.pxSize
            
        A = self.radius1
        B = self.radius2
        C = 2.0*self.pxSize
        rad = pow(pow(tx/A, 2)+pow(ty/B, 2)+pow(tz/C, 2), 0.5)
        polVal = self.polynomByRochSmets(rad)
            
        if (rad < 1.0):
            return polVal*self.thermVion

        return self.ZERO
    
    def getVelocityY4InjectedParticles(self, x, y, z):
        return getVelocityX4InjectedParticles(x, y, z)
    
    def getVelocityZ4InjectedParticles(self, x, y, z):
        return getVelocityX4InjectedParticles(x, y, z)

    # set density profile to sustain by ablation operator
    def getTargetIonDensity2sustain(self, x, y, z):
        Ly2 = 0.5*self.boxSize[1]
        Lz = self.boxSize[2]
        if (z < self.surface_pos and z > (self.surface_pos - self.trgtwidth) ):
            return self.dens2sustain
        else:
            return self.ZERO


    # smooth polynom by roch smets 2014 PoP
    def polynomByRochSmets(self, x): # x = |x|
        return -6.0*x**5+15.0*x**4-10.0*x**3+1.0
    
    
     # set electron pressure profile to sustain by ablation operator
    def getElectronPressure2sustain(self, x, y, z):
        Ly2 = 0.5*self.boxSize[1]
        Lz = self.boxSize[2]
        if (z > self.surface_pos):
            return self.ZERO

        tx = x - self.spotPos[0]
        ty = y - self.spotPos[1]
        tz = z - self.spotPos[2]+self.pxSize
            
        A = self.radius1
        B = self.radius2
        C = 2.0*self.pxSize
        rad = pow(pow(tx/A, 2)+pow(ty/B, 2)+pow(tz/C, 2), 0.5)
        polVal = self.polynomByRochSmets(rad)
        
            
        if (rad < 1.0):
            return polVal*self.dens2sustain*self.temp2sustain

        return self.ZERO


   
    #   physics: magnetic field
    def getBfieldX(self, x, y, z):
        return self.Bfield[0]
    
    def getBfieldY(self, x, y, z):
        return self.Bfield[1]
    
    def getBfieldZ(self, x, y, z):
        return self.Bfield[2]
    
    #   physics: ohm's law resistivity
    def getResistivity(self):
        return self.resistivity    
    
    #   physics: pressure evolution
    def getElectronMass(self):
        return self.emass

    #   physics: pressure evolution : isotropization
    def getRelaxFactor(self):
        return self.relaxFactor

    #   physics: pressure evolution : smoothing
    def getElectronPressureSmoothingStride(self):
        return self.smoothStride

    #   physics: pressure evolution : 1 - isothermal (optional), 0 - evolution equation
    def getIfWeUseIsothermalClosure(self):
        return 0
    #   physics: isothemal closure : electron temperature
    def getElectronTemperature4IsothermalClosure(self):
        return 0.0










