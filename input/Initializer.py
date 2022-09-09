
class Initializer:

    
    def __init__(self):
        #   0: run from scratch
        #   1: restart, set input file
        self.runType = 0
        self.inputFile = "./output/_restart.h5"
        
        # BOX
        dx = 0.2
        self.boxSize    = [60.0, 20.0, dx] # in d0 - ion inertia length
        self.boxSizePxl = [self.boxSize[0]/dx, self.boxSize[1]/dx, 1]

        self.pxSize = self.boxSize[0]/self.boxSizePxl[0]
        
        self.bcType = [0,1,1] #1 - periodic 0 - ideal
        
        self.partclBcType = [0,1,1] #2 - reflect 1 - periodic 0 - outflow
        
        self.dampingBoundaryWidth = [[1,2], [0,0], [0,0]]
         
        self.mpiCores  = [8,4,1]
        
        # time
        self.ts = 0.005
        self.maxtsnum = 8001
        self.outputStride = 800

        # output. need to create it before
        self.outputDir = "/gpfs/fs1/home/a/anticipa/weipeng/scratch/Data_aka52/stream/output/"
        self.fileTemplate = "test_"
    
        # Particles
        self.ppc4load = 30
        self.ppc = [10, 10]
        self.ppcMinDens = 1.0

        self.numOfSpecies = 2
        self.masses  = [1, 20]   # [mp]
        self.charges = [1, 10]   # [qe]
        self.dens = [0.1, 1.0] # [nc]
        self.vel1 = 0.1         # [c?] -> ~ 5 MeV
        self.vfl1 = [0.0,0.0,0.0]
        self.vel2 = 0.0
        self.vfl2 = [0.0,0.0,0.0]

        self.DFtype4species1 = 1           # 1 - rectangular distribution function, since v1.6
        self.DFtype4InjectedParticles = 1  # same above
       
        self.Pele0 = 0.0000001  # [nc*eV?]
        self.trgtwidth = 2.0    # [di]
        # Laser spots
        self.spotsNum = 1                

        self.type2Load = 2      # [?]
        self.thermVion = 0.1    # [c?]
	
        self.vfl2Load = [1.,0.0,0.0] # [?]	

        self.dens2sustain = self.dens[1]
    
        #magnetic field magnitude
        self.Bfield = [0.00625, 0.0, 0.0] # [so that wce = wpe? 1e18 cc -> 320 T; 1e19 cc -> 1000 T]
                                          # [B = 20 T & ne = 1e20 cc -> B_norm = 0.00625]   
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

    def getElectronPressureXX(self, x, y, z):
        return self.Pele0

    def getElectronPressureYY(self, x, y, z):
        return self.Pele0
    
    def getElectronPressureZZ(self, x, y, z):
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

    def getDensity4species1(self, x, y, z):
        return self.dens[0]

    def getIfParticleTypeIsFrozen4species1(self):
        return self.ZERO

    # species 1: modulus of velocity for Maxwell distribution
    def getVelocityX4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityY4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityZ4species1(self, x, y, z):
        return self.vel1


    # species 1: modulus of fluid velocity
    def getFluidVelocityX4species1(self, x, y, z):
        return self.vfl1[0]
    
    def getFluidVelocityY4species1(self, x, y, z):
        return self.vfl1[1]
    
    def getFluidVelocityZ4species1(self, x, y, z):
        return self.vfl1[2]
    
    #           species 2
    def getPPC4species2(self):
        return self.ppc[1]
    
    def getMass4species2(self):
        return self.masses[1]
    
    def getCharge4species2(self):
        return self.charges[1]
        
    def getDensity4species2(self, x, y, z):
        if x <  self.trgtwidth:
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
        return self.vfl2[0]
    
    def getFluidVelocityY4species2(self, x, y, z):
        return self.vfl2[1]
    
    def getFluidVelocityZ4species2(self, x, y, z):
        return self.vfl2[2]

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
        return self.vfl2Load[0]
    
    def getFluidVelocityY4InjectedParticles(self, x, y, z):
        return self.vfl2Load[1]
    
    def getFluidVelocityZ4InjectedParticles(self, x, y, z):
        return self.vfl2Load[2]

    # injected particles: modulus of thermal velocity
    def getVelocityX4InjectedParticles(self, x, y, z):
        return self.thermVion
    
    def getVelocityY4InjectedParticles(self, x, y, z):
        return self.getVelocityX4InjectedParticles(x, y, z)
    
    def getVelocityZ4InjectedParticles(self, x, y, z):
        return self.getVelocityX4InjectedParticles(x, y, z)

    # set density profile to sustain by ablation operator
    def getTargetIonDensity2sustain(self, x, y, z):
        if (x < self.trgtwidth):
            return self.dens2sustain
        else:
            return self.ZERO

    #v1.6 change distribution function of particle specie to rectangular
    def getDFtype4species1(self, x, y, z):
        return 1
        # return self.DFtype4species1  # to be better formatted.

    def getDFtype4InjectedParticles(self, x, y, z):
        return 1
        # return self.DFtype4InjectedParticles  


    # smooth polynom by roch smets 2014 PoP
    def polynomByRochSmets(self, x): # x = |x|
        return -6.0*x**5+15.0*x**4-10.0*x**3+1.0
    
    
     # set electron pressure profile to sustain by ablation operator
    def getElectronPressure2sustain(self, x, y, z):
        if (x < self.trgtwidth):
            return self.Pele0
        else:
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
        return 1

    #   physics: isothemal closure : electron temperature
    def getElectronTemperature4IsothermalClosure(self):
        return self.Pele0


    #   v1.6 physics: ion-ion collisions:
    def getDefaultColoumbLogarithm(self):
        return 10

    def getIonIonCollisionFrequencyFactor(self):
        return 1    









