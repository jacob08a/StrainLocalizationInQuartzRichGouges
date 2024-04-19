readParamsFromTable(
    SN=5.e6,
    coeff=0.0001,
    RATE_NS1=0.1,
    RATE_NS2=0.1,   #original 0.1
    RATE_shear=1,
    TSSC=0.3,
    shearStrain=5,
    wallFRIC=0.0, # degree
    compFRIC=27, # degree, during compaction (controls porosity)
    shearFRIC=27, # degree, during shear
    DAMPSHEAR=0.,
    roll_stiff=0.0,
    roll_fric=0.0,
    Ystep=100,
    Zstep=50,
    OUT='frictionlessWall_Hertz_single_thickness15d_ss5_compFRIC27(mu05)_shearFRIC27(mu05)_TSSC03_withaccele_221216'
    )

from yade.params.table import *
from yade import pack,plot,export,ymport
from scipy import stats
import numpy as np
import math

sp1=pack.SpherePack()
sp2=pack.SpherePack()
sp3=pack.SpherePack()

O.periodic=True

# dimensions of sample (fixed by particle size such as L/D~15)
DIAMETER_packing=0.000250
RADIUS_packing=0.000125

cellSize_packing_502525 = np.load('cellSize_single_125r_sameSize.npy')
cellSize_boundary_502525 = np.load('cellSize_single_boundary_125r_sameSize.npy')

length=cellSize_boundary_502525[0]
height=cellSize_boundary_502525[1]/3
width=cellSize_boundary_502525[2]

# friction angles (see above)

# boundary conditions
PI=1.e5
conf = 0.
if conf<PI:
    PC = PI
else :
    PC = conf
SN=SN # normal stress
RATE_NS1=RATE_NS1 # velocity of top plate during compaction
RATE_NS2=RATE_NS2 # velocity of top plate during shear
RATE_shear=RATE_shear # shear velocity

# simulation control
ITER=2e5
OUT=OUT

#####Mair and Hazzard (2007)'s microproperties
O.cell.hSize=Matrix3(length,0,0,0,3*height,0,0,0,width)

#O.materials.append(FrictMat(density=2500,young=5.5e10,poisson=0.25,frictionAngle=wallFRIC,label='boxMat'))
#O.materials.append(FrictMat(density=2500,young=5.5e10,poisson=0.25,frictionAngle=boundaryFRIC,label='boundaryMat'))
#O.materials.append(FrictMat(density=2500,young=5.5e10,poisson=0.25,frictionAngle=spFRIC,label='sphereMat'))
O.materials.append(FrictMat(density=2500,young=5.5e10,poisson=0.25,frictionAngle=radians(wallFRIC),label='boxMat'))
O.materials.append(FrictMat(density=2500,young=5.5e10,poisson=0.25,frictionAngle=radians(compFRIC),label='boundaryMat'))
O.materials.append(FrictMat(density=2500,young=5.5e10,poisson=0.25,frictionAngle=radians(compFRIC),label='sphereMat'))

upBox = utils.box(center=(length/2,1.8875*height+1.1*DIAMETER_packing,1.5*width),orientation=Quaternion(1,0,0,0),extents=(5*length,DIAMETER_packing/2,5*width),fixed=1,wire=False,color=(1,0,0),material='boxMat')
lowBox = utils.box(center=(length/2,1.1125*height-1.1*DIAMETER_packing,1.5*width),orientation=Quaternion(1,0,0,0),extents=(5*length,DIAMETER_packing/2,5*width),fixed=1,wire=False,color=(1,0,0),material='boxMat')
frontBox = utils.box(center=(length/2,1.8875*height+1.1*DIAMETER_packing,1.5*width),orientation=Quaternion(1,0,0,0),extents=(5*length,DIAMETER_packing/2,5*width),fixed=1,wire=False,color=(1,0,0),material='boxMat')
backBox = utils.box(center=(length/2,1.8875*height+1.1*DIAMETER_packing,1.5*width),orientation=Quaternion(1,0,0,0),extents=(5*length,DIAMETER_packing/2,5*width),fixed=1,wire=False,color=(1,0,0),material='boxMat')

O.bodies.append([upBox,lowBox])

### Import isotropic-compressed packing and boundary
packing = ymport.text("Packing_single_125r_50-25-25.isoCompression",color=(1,1,0),material='sphereMat')
boundary_top = ymport.text("Packing_single_boundary_125r_50-25-25.isoCompression",color=(148/255,0,211/255),material='boundaryMat')
boundary_bottom = ymport.text("Packing_single_boundary_125r_50-25-25.isoCompression",color=(148/255,0,211/255),material='boundaryMat')

### Filter only desired packing
def sphereWanted(b):
    x,y,z = b.state.pos
    return x >= 0 and x <= length and y >= 1.1625*height and y <= 1.8375*height and z >= 0 and z <= width

def boundaryWanted(b):
    x,y,z = b.state.pos
    return x >= 0 and x <= cellSize_packing_502525[0] and z >= 0 and z <= cellSize_packing_502525[2]

packing = [b for b in packing if sphereWanted(b)]
boundary_top = [b for b in boundary_top if boundaryWanted(b)]
boundary_bottom = [b for b in boundary_bottom if boundaryWanted(b)]

### Shift position
cellSize_deviation_x = abs(cellSize_boundary_502525[0]/2 - cellSize_packing_502525[0]/2)
cellSize_deviation_y = abs(cellSize_boundary_502525[1]/2 - cellSize_packing_502525[1]/2)
cellSize_deviation_z = abs(cellSize_boundary_502525[2]/2 - cellSize_packing_502525[2]/2)

for b in packing:
    b.state.pos += ((cellSize_deviation_x),0,(cellSize_deviation_z))
for b in boundary_top:
    b.state.pos += (-(cellSize_deviation_x),-(cellSize_boundary_502525[1]/2-1.8875*height),-(cellSize_deviation_z))
for b in boundary_bottom:
    b.state.pos += ((cellSize_deviation_x),-(cellSize_boundary_502525[1]/2-1.1125*height),(cellSize_deviation_z))

O.bodies.append(packing)
O.bodies.append(boundary_top)
O.bodies.append(boundary_bottom)

list_packingThickness = []
for b in packing:
    list_packingThickness.append(b.state.pos[1])

print ('generated layer thickness = ', max(list_packingThickness)-min(list_packingThickness))

effCellVol=(O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2]
volRatio=(O.cell.hSize[0,0]*O.cell.hSize[1,1]*O.cell.hSize[2,2])/effCellVol

print ('dataName=',OUT)
print ('volRatio=',volRatio)

O.engines=[
 ForceResetter()
 ,InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb()],verletDist=-0.1,allowBiggerThanPeriod=True)
 ,InteractionLoop(
  [Ig2_Sphere_Sphere_ScGeom6D(),Ig2_Box_Sphere_ScGeom6D()],
  [Ip2_FrictMat_FrictMat_MindlinPhys(krot=roll_stiff,eta=roll_fric)],
  [Law2_ScGeom_MindlinPhys_Mindlin(includeMoment=True)]
 )
 ,GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=1,timestepSafetyCoefficient=TSSC,defaultDt=-1)
 ,PyRunner(command='fixVelocity(RATE_shear)',iterPeriod=1,label='fixVel',dead=True)
 ,PeriTriaxController(dynCell=True,mass=10,maxUnbalanced=1e-3,relStressTol=1e-4,stressMask=7,goal=(-PI/volRatio,-PI/volRatio,-PI/volRatio),globUpdate=1,maxStrainRate=(1,1,1),doneHook='triaxDone()',label='triax')
 ,NewtonIntegrator(gravity=(0,0,0),damping=0.3,label='newton')
 ,PyRunner(command='dataRecorder_preshear()',iterPeriod=1,label='recData_preshear',dead=True)
 ,PyRunner(command='dataRecorder()',iterPeriod=1,label='recData',dead=True)
 ,PyRunner(iterPeriod=100,command="vtkExport()",label='vtkExp',dead=True)
 #,VTKRecorder(fileName='3d-vtk-',recorders=['all'],iterPeriod=1000,label='recVTK',dead=True)
 #,PyRunner(command='fixVelocity(RATE_shear)',iterPeriod=1,label='fixVel',dead=True)
 ]
 
TW=TesselationWrapper()
TW.setState()
TW.computeVolumes()
s=bodyStressTensors()

vtk_spheres = export.VTKExporter(OUT)
vtk_interactions = export.VTKExporter(OUT)
        
def vtkExport():
    list_intrs = []
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            vol_sph = 4.*pi/3.*b.shape.radius**3
            b.porosity = 1- (vol_sph/TW.volume(b.id))
            b.momentum = b.state.mass * b.state.vel
            b.linVel = b.state.vel
            b.angVel = b.state.angVel
            b.force = O.forces.f(b.id)
            b.torque = O.forces.t(b.id)
            b.color = b.shape.color
        for i in b.intrs():
            O.interactions.erase(i.id1,0)
            O.interactions.erase(i.id1,1)
            O.interactions.erase(0,i.id2)
            O.interactions.erase(1,i.id2)
            if i.isReal:
                list_intrs.append([i.id1,i.id2])
    vtk_spheres.exportSpheres(what=dict(color="b.color",porosity="b.porosity",momentum="b.momentum",linearVel="b.linVel",angVel="b.angVel",force="b.force",torque="b.torque"))
    vtk_interactions.exportInteractions(ids=list_intrs,what=dict(contactNormalForce="i.phys.normalForce", contactShearForce="i.phys.shearForce"))

def dataRecorder_preshear():
 h=vol=vol_s=nb_s=nbSph=Rmean=Rmax=0.
 Rmin=1e6
 h=O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
 vol=h*O.cell.hSize[0,0]*O.cell.hSize[2,2]
 contactStress=getStress(vol)
 for o in O.bodies:
  if isinstance(o.shape,Sphere):
   nbSph+=1
   Rmean+=o.shape.radius
   if o.shape.radius>Rmax: Rmax=o.shape.radius
   if o.shape.radius<Rmin: Rmin=o.shape.radius
   vol_s += 4.*pi/3.*(o.shape.radius)**3
 Rmean=Rmean/nbSph
 n = 1-vol_s/vol
 plot.addData(
  iter=O.iter
  ,time=O.time
  ,dt=O.dt
  ,shear_strain=(O.bodies[0].state.pos[0])/h
  ,shear_strain2=boundary_top[0].state.pos[0]/h
  ,TW_x=O.bodies[0].state.pos[0]
  ,TW_y=O.bodies[0].state.pos[1]
  ,BW_y=O.bodies[1].state.pos[1]
  ,thickness=h
  ,porosity=n
  )
 
def dataRecorder():
 h=vol=vol_s=nb_s=nbSph=Rmean=Rmax=0.
 Rmin=1e6
 h=O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
 vol=h*O.cell.hSize[0,0]*O.cell.hSize[2,2]
 contactStress=getStress(vol)
 for o in O.bodies:
  if isinstance(o.shape,Sphere):
   nbSph+=1
   Rmean+=o.shape.radius
   if o.shape.radius>Rmax: Rmax=o.shape.radius
   if o.shape.radius<Rmin: Rmin=o.shape.radius
   vol_s += 4.*pi/3.*(o.shape.radius)**3
 Rmean=Rmean/nbSph
 n = 1-vol_s/vol
 totalShearForce_x=0
 totalShearForce_y=0  
 totalShearForce_z=0 
 for b in boundary_top:  
  totalShearForce_x+=O.forces.f(b.id)[0]
  totalShearForce_y+=O.forces.f(b.id)[1]
  totalShearForce_z+=O.forces.f(b.id)[2]
 plot.addData(
  iter=O.iter
  ,time=O.time
  ,dt=O.dt
  ,stress_topsp_x=abs(O.forces.f(0)[0]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
  ,stress_topsp_y=abs(O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
  ,stress_topsp_z=abs(O.forces.f(0)[2]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
  ,shearStress=abs(totalShearForce_x/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
  ,Friction_coefficient=abs(totalShearForce_x/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))/abs(O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
  ,shear_strain=(O.bodies[0].state.pos[0])/h
  ,shear_strain2=boundary_top[0].state.pos[0]/h
  ,TW_x=O.bodies[0].state.pos[0]
  ,TW_y=O.bodies[0].state.pos[1]
  ,BW_y=O.bodies[1].state.pos[1]
  ,thickness=h
  ,porosity=n
  )
 
phase = 0
def triaxDone():
 global phase
 volRatio=(O.cell.hSize[0,0]*O.cell.hSize[1,1]*O.cell.hSize[2,2])/((O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2])
 if phase == 0:
  h=O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
  vol=h*O.cell.hSize[0,0]*O.cell.hSize[2,2]
  contactStress=getStress(vol)
  vol_s=Rmean=Rmax=nbSph=0
  Rmin=1e6
  for o in O.bodies:
    if isinstance(o.shape,Sphere):
      nbSph+=1
      Rmean+=o.shape.radius
      if o.shape.radius>Rmax: Rmax=o.shape.radius
      if o.shape.radius<Rmin: Rmin=o.shape.radius
      vol_s += 4.*pi/3.*(o.shape.radius)**3
  Rmean=Rmean/nbSph
  n = 1-vol_s/vol
  print('DONE! iter=',O.iter,'| sample generated: nb spheres',nbSph,', Rmean=',Rmean,', Rratio=',Rmax/Rmin,', porosity=',n)
  print('Changing contact properties now')
  #tt=TriaxialCompressionEngine()
  #tt.setContactProperties(FRIC)
  utils.setContactFriction(radians(shearFRIC))
  print("APPLYING CONFINING PRESSURE : sx ,sy and sz will go to PC =",PC)
  triax.goal=(-PC/volRatio,-PC/volRatio,-PC/volRatio)
  phase+=1
 elif phase ==1:
  print ('DONE ! iter =',O.iter ,'| isotropic confinement done : stresses =',volRatio * triax.stress )
  triax.dead=True
  O.pause()

#### Initialization
print ('SAMPLE PREPARATION!')

recData_preshear.dead = False # uncomment if you want to record what is happening during preparation of sample (isotropic compaction)
#O.run (1000000 ,1)
#recVTK.dead = False
#O.step()
#recVTK.dead = True
O.run(10000000000,1)
#plot.saveDataTxt(OUT+"_isoConfined_"+str(O.iter)+'.txt')

print ('Normal stress (platen) = ',O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
print ('Normal stress (contacts) = ',getStress((O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2])[1,1])

#### Applying normal stress
print ('NORMAL LOADING! iter=',O.iter)

stage=0
stiff=fnPlaten=currentSN=0.
def servo():
 global stage,stiff,fnPlaten,currentSN,TSSC
 if stage==0:
  currentSN=O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2])
  unbF=unbalancedForce()
  #print ('SN=',SN,'| current SN = ',currentSN,' | unbF=',unbF )
  boundaryVel=copysign(min(RATE_NS1,abs(0.5*(currentSN-SN))),currentSN-SN)
  #print ('boundaryVel=',boundaryVel)
  O.bodies[0].state.vel[1]=boundaryVel
  if ( (abs(currentSN-SN)/SN)<0.001 and unbF<0.001 ):
   stage+=1
   fnPlaten=O.forces.f(0)[1]
   print ('Normal stress =',currentSN,' | unbF=',unbF)
   ## the following computes the stiffness of the plate (used for stress control of the top plate)
   for i in O.interactions.withBody(O.bodies[0].id):
    stiff+=i.phys.kn
   print ('DONE! iter=',O.iter)
   O.pause()
 if stage==1:
  fnDesired=SN*(O.cell.hSize[0,0]*O.cell.hSize[2,2])
  #boundaryVel_2=copysign(abs(0.5*(O.forces.f(0)[1]-fnDesired)/stiff/O.dt),O.forces.f(0)[1]-fnDesired)
  boundaryVel=copysign(min(RATE_NS2,abs(0.333*(O.forces.f(0)[1]-fnDesired)/stiff/O.dt/TSSC)),O.forces.f(0)[1]-fnDesired)
  #print('check: ', RATE_NS2 == min(RATE_NS2,abs(0.333*(O.forces.f(0)[1]-fnDesired)/stiff/O.dt)))
  #print ('boundaryVel=',boundaryVel)
  O.bodies[0].state.vel[1]=boundaryVel
  #plot.addData(
  # wall_velocity=boundaryVel
  # )
  
O.engines = O.engines[:5]+[PyRunner(command='servo()',iterPeriod=1,label='servo')]+O.engines[5:]

O.trackEnergy = True

O.run(10000000000,1)
plot.saveDataTxt(OUT+"_normalConfined_"+str(O.iter)+'.txt')

print ('STABILIZING! iter=',O.iter)

# coloring particles to see vertical stripes in material
list_marker = [0,2]
dxi=8*(2*RADIUS_packing) # can be a function of cell length dxi=O.cell.hSize[0,0]/5.
n=int(length/dxi)
xmin=1e6
for i in list_marker:
 for o in O.bodies:
  if o.id>1:
   if o.state.pos[0]<xmin: xmin=o.state.pos[0]
   if (o.state.pos[0]>=(i+0.3)*dxi) and (o.state.pos[0]<((i+0.8)*dxi)):
    o.shape.color=(1,0,0)

print ('Normal stress (platen) = ',O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
print ('Normal stress (contacts) = ',getStress((O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2])[1,1])

################################################################ Start shearing
#### Initial X position of all spheres
initialposition_x = []
for i in range(2,len(O.bodies)):
    initialposition_x.append(O.bodies[i].state.pos[0])
    
#### Initial position of all spheres
initialposition = []
for i in range(2,len(O.bodies)):
    initialposition.append(O.bodies[i].state.pos)

#### fixed bottom bottomLayer_id
for b in boundary_bottom:
    b.state.blockedDOFs='xyzXYZ'
    b.state.dynamics = False
    b.state.vel = Vector3(0,0,0)
    
refIteration = O.iter
print ('SHEARING! iter=',O.iter)

refTime = O.time
print ("Simulation_time = ", O.time)

print("Layer_thickness = ", O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])

initial_particle_pos = boundary_top[0].state.pos[0] ### The initial position of one of the topBoundaryLayer particles

recData_preshear.dead = True
recData.dead=False
#vtkExp.dead=False
#recVTK.dead=False
newton.damping=DAMPSHEAR
fixVel.dead = False

def fixVelocity(RATE_shear):
    shearVel = 0
    O.bodies[0].state.blockedDOFs='x'
    if shearVel < RATE_shear :
        shearVel +=( RATE_shear /100.)
    O.bodies[0].state.vel = (shearVel,0,0)
    for b in boundary_top:
        b.state.blockedDOFs='x'
        b.state.vel=(shearVel,0,0)
    slip = boundary_top[0].state.pos[0] - initial_particle_pos
    h=O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
    ss = slip/h
    if ss > shearStrain:
        O.pause()

O.run(10000000000,1)
plot.saveDataTxt(OUT+'.txt')

### Contact force network visualization
#from yade import qt

#qt.View()

#rr=yade.qt.Renderer()
#rr.shape=False
#rr.intrPhys=True

################################################################################ Post processing
def weight_function_displacement(z,z1,coeff):
    return(math.exp(-(((z-z1)**2)/(2*(coeff**2))))) 

def weight_function_linearVel(v,v1,coeff):
    return(math.exp(-(((v-v1)**2)/(2*(coeff**2))))) 

def weight_function_force(y,y1,coeff):
    return(math.exp(-(((y-y1)**2)/(2*(coeff**2)))))

def weight_function_force_2D(y2d,y2d1,coeff):
    return(math.exp(-(((y2d-y2d1)**2)/(2*(coeff**2)))))

def weight_function_contactNormalForce(n,n1,coeff):
    return(math.exp(-(((n-n1)**2)/(2*(coeff**2))))) 

def weight_function_contactShearForce(s,s1,coeff):
    return(math.exp(-(((s-s1)**2)/(2*(coeff**2))))) 

def weight_function_torque(t,t1,coeff):
    return(math.exp(-(((t-t1)**2)/(2*(coeff**2))))) 

def weight_function_momentum_x(mx,mx1,coeff):
    return(math.exp(-(((mx-mx1)**2)/(2*(coeff**2))))) 

def weight_function_momentum_y(my,my1,coeff):
    return(math.exp(-(((my-my1)**2)/(2*(coeff**2))))) 

def weight_function_momentum_z(mz,mz1,coeff):
    return(math.exp(-(((mz-mz1)**2)/(2*(coeff**2))))) 

def weight_function_microStress(a,a1,coeff):
    return(math.exp(-(((a-a1)**2)/(2*(coeff**2)))))

def weight_function_porosity(p,p1,coeff):
    return(math.exp(-(((p-p1)**2)/(2*(coeff**2)))))

def weight_function_angVel_x(avx,avx1,coeff):
    return(math.exp(-(((avx-avx1)**2)/(2*(coeff**2)))))

def weight_function_angVel_y(avy,avy1,coeff):
    return(math.exp(-(((avy-avy1)**2)/(2*(coeff**2)))))

def weight_function_angVel_z(avz,avz1,coeff):
    return(math.exp(-(((avz-avz1)**2)/(2*(coeff**2)))))

def weight_function_angMomentum_x(amx,amx1,coeff):
    return(math.exp(-(((amx-amx1)**2)/(2*(coeff**2)))))

def weight_function_angMomentum_y(amy,amy1,coeff):
    return(math.exp(-(((amy-amy1)**2)/(2*(coeff**2)))))

def weight_function_angMomentum_z(amz,amz1,coeff):
    return(math.exp(-(((amz-amz1)**2)/(2*(coeff**2)))))

#### Final X position of all spheres
finalposition_x = []           # Displacement of each particle in X direction
for i in range(2,len(O.bodies)):
    finalposition_x.append(O.bodies[i].state.pos[0])
    
#### Final position of all spheres
finalposition = []           # Displacement of each particle in X direction
for i in range(2,len(O.bodies)):
    finalposition.append(O.bodies[i].state.pos)
 
#### Slip displacement in X direction
displacement_x = []
for i in range(2,len(O.bodies)):
    displacement_x.append(O.bodies[i].state.pos[0] - initialposition_x[i-2])
    
def getDisplacement(z1):
  weight=0 
  weightedDisplacement=0
  totalWeightedDisplacement=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_displacement(O.bodies[i].state.pos[1], z1, coeff)     # Weight function coefficient
    weightedDisplacement = displacement_x[i-2] * weight         # Total weighted particle displacement at certain Y position
    totalWeightedDisplacement = totalWeightedDisplacement + weightedDisplacement
    totalWeight = totalWeight + weight                                       # All weighted particle displacement that is involved in at certain Y position
  return (totalWeightedDisplacement/totalWeight)  

def getLinearVel(v1):
  weight=0 
  weightedLinearVel=0
  totalWeightedLinearVel=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_linearVel(O.bodies[i].state.pos[1], v1, coeff)    
    weightedLinearVel = O.bodies[i].state.vel[0] * weight       
    totalWeightedLinearVel = totalWeightedLinearVel + weightedLinearVel
    totalWeight = totalWeight + weight                                       
  return (totalWeightedLinearVel/totalWeight)  

def getForce(y1):
  weight=0
  weightedForce=0
  totalWeightedForce=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_force(O.bodies[i].state.pos[1], y1, coeff)
    weightedForce = O.forces.f(i)[0] * weight
    totalWeightedForce = totalWeightedForce + weightedForce
    totalWeight = totalWeight + weight
  return (totalWeightedForce/totalWeight)

def getForce_2D(y2d1):
  weight=0
  weightedForce_2d=0
  totalWeightedForce_2d=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_force_2D(O.bodies[i].state.pos[2], y2d1, coeff)
    weightedForce_2d = O.forces.f(i)[0] * weight
    totalWeightedForce_2d = totalWeightedForce_2d + weightedForce_2d
    totalWeight = totalWeight + weight
  return (totalWeightedForce_2d/totalWeight)

def getTorque(t1):
  weight=0
  weightedTorque=0
  totalWeightedTorque=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_torque(O.bodies[i].state.pos[1], t1, coeff)
    weightedTorque = O.forces.t(i)[0] * weight
    totalWeightedTorque = totalWeightedTorque + weightedTorque
    totalWeight = totalWeight + weight
  return (totalWeightedTorque/totalWeight)

def getContactNormalForce(n1):
  weight=0 
  weighedContactNormalForce=0
  totalContactNormalForce=0
  totalWeight=0
  for i in O.interactions:
      if i.isReal:
          weight = weight_function_contactNormalForce(i.geom.contactPoint[1], n1, coeff)     # Weight function coefficient
          weighedContactNormalForce = (i.phys.normalForce[0]) * weight       
          totalContactNormalForce = totalContactNormalForce + weighedContactNormalForce
          totalWeight = totalWeight + weight                                       #
  return (totalContactNormalForce/totalWeight)  

def getContactShearForce(s1):
  weight=0 
  weighedContactShearForce=0
  totalContactShearForce=0
  totalWeight=0
  for i in O.interactions:
      if i.isReal:
          weight = weight_function_contactShearForce(i.geom.contactPoint[1], s1, coeff)     # Weight function coefficient
          weighedContactShearForce = (i.phys.shearForce[0]) * weight       
          totalContactShearForce = totalContactShearForce + weighedContactShearForce
          totalWeight = totalWeight + weight                                       #
  return (totalContactShearForce/totalWeight) 

def getMomentum_x(mx1):
  weight=0 
  weightedMomentum_x=0
  totalWeightedMomentum_x=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_momentum_x(O.bodies[i].state.pos[1], mx1, coeff)     
    weightedMomentum_x = O.bodies[i].state.mass * O.bodies[i].state.vel[0] * weight         
    totalWeightedMomentum_x = totalWeightedMomentum_x + weightedMomentum_x
    totalWeight = totalWeight + weight                                       
  return (totalWeightedMomentum_x/totalWeight) 

def getMomentum_y(my1):
  weight=0 
  weightedMomentum_y=0
  totalWeightedMomentum_y=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_momentum_y(O.bodies[i].state.pos[1], my1, coeff)     
    weightedMomentum_y = O.bodies[i].state.mass * O.bodies[i].state.vel[1] * weight         
    totalWeightedMomentum_y = totalWeightedMomentum_y + weightedMomentum_y
    totalWeight = totalWeight + weight                                       
  return (totalWeightedMomentum_y/totalWeight) 

def getMomentum_z(mz1):
  weight=0 
  weightedMomentum_z=0
  totalWeightedMomentum_z=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_momentum_z(O.bodies[i].state.pos[1], mz1, coeff)     
    weightedMomentum_z = O.bodies[i].state.mass * O.bodies[i].state.vel[2] * weight         
    totalWeightedMomentum_z = totalWeightedMomentum_z + weightedMomentum_z
    totalWeight = totalWeight + weight                                       
  return (totalWeightedMomentum_z/totalWeight) 


def getMicroStress(a1):
  weight=0 
  weightedMicroStress=0
  totalWeightedMicroStress=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_microStress(O.bodies[i].state.pos[1], a1, coeff)
    mystress = s[i]*4.*pi/3.*O.bodies[i].shape.radius**3/TW.volume(i)
    weightedMicroStress = mystress[0][0] * weight
    totalWeightedMicroStress = totalWeightedMicroStress + weightedMicroStress
    totalWeight = totalWeight + weight
  return (totalWeightedMicroStress/totalWeight)


TW=TesselationWrapper()
TW.setState()
TW.computeVolumes()
s=bodyStressTensors()

def getPorosity(p1):
  weight=0 
  weightedPorosity=0
  totalWeightedPorosity=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_porosity(O.bodies[i].state.pos[1], p1, coeff)
    vol_sph = 4.*pi/3.*O.bodies[i].shape.radius**3
    porosity = 1 - (vol_sph/TW.volume(i))
    weightedPorosity = porosity * weight
    totalWeightedPorosity = totalWeightedPorosity + weightedPorosity
    totalWeight = totalWeight + weight
  return (totalWeightedPorosity/totalWeight)

def getAngVel_x(avx1):
  weight=0 
  weightedAngVel_x=0
  totalWeightedAngVel_x=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_angVel_x(O.bodies[i].state.pos[1], avx1, coeff)
    weightedAngVel_x = O.bodies[i].state.angVel[0] * weight       
    totalWeightedAngVel_x = totalWeightedAngVel_x + weightedAngVel_x
    totalWeight = totalWeight + weight                                       
  return (totalWeightedAngVel_x/totalWeight) 

def getAngVel_y(avy1):
  weight=0 
  weightedAngVel_y=0
  totalWeightedAngVel_y=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_angVel_y(O.bodies[i].state.pos[1], avy1, coeff)
    weightedAngVel_y = O.bodies[i].state.angVel[1] * weight       
    totalWeightedAngVel_y = totalWeightedAngVel_y + weightedAngVel_y
    totalWeight = totalWeight + weight                                       
  return (totalWeightedAngVel_y/totalWeight)

def getAngVel_z(avz1):
  weight=0 
  weightedAngVel_z=0
  totalWeightedAngVel_z=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_angVel_z(O.bodies[i].state.pos[1], avz1, coeff)
    weightedAngVel_z = O.bodies[i].state.angVel[2] * weight       
    totalWeightedAngVel_z = totalWeightedAngVel_z + weightedAngVel_z
    totalWeight = totalWeight + weight                                       
  return (totalWeightedAngVel_z/totalWeight)

def getAngMomentum_x(amx1):
  weight=0 
  weightedAngMomentum_x=0
  totalWeightedAngMomentum_x=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_angMomentum_x(O.bodies[i].state.pos[1], amx1, coeff)
    weightedAngMomentum_x = O.bodies[i].state.angVel[0] * O.bodies[i].state.mass * O.bodies[i].shape.radius**2 * weight       
    totalWeightedAngMomentum_x = totalWeightedAngMomentum_x + weightedAngMomentum_x
    totalWeight = totalWeight + weight                                       
  return (totalWeightedAngMomentum_x/totalWeight)

def getAngMomentum_y(amy1):
  weight=0 
  weightedAngMomentum_y=0
  totalWeightedAngMomentum_y=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_angMomentum_y(O.bodies[i].state.pos[1], amy1, coeff)
    weightedAngMomentum_y = O.bodies[i].state.angVel[1] * O.bodies[i].state.mass * O.bodies[i].shape.radius**2 * weight       
    totalWeightedAngMomentum_y = totalWeightedAngMomentum_y + weightedAngMomentum_y
    totalWeight = totalWeight + weight                                       
  return (totalWeightedAngMomentum_y/totalWeight)

def getAngMomentum_z(amz1):
  weight=0 
  weightedAngMomentum_z=0
  totalWeightedAngMomentum_z=0
  totalWeight=0
  for i in range(2,len(O.bodies)):
    weight = weight_function_angMomentum_z(O.bodies[i].state.pos[1], amz1, coeff)
    weightedAngMomentum_z = O.bodies[i].state.angVel[2] * O.bodies[i].state.mass * O.bodies[i].shape.radius**2 * weight       
    totalWeightedAngMomentum_z = totalWeightedAngMomentum_z + weightedAngMomentum_z
    totalWeight = totalWeight + weight                                       
  return (totalWeightedAngMomentum_z/totalWeight)

dl = (O.cell.hSize[2][2]) / Zstep  #Thickness of particle length divided over 100 steps


list_force_2D = []
for i in range(Zstep):
    y2d1 = dl*i
    list_force_2D.append([getForce_2D(y2d1),y2d1])

dh = (O.bodies[0].state.pos[1] - O.bodies[1].state.pos[1]) / Ystep  #Thickness of particle bed divided over 100 steps

list_displacement = []
for i in range(Ystep):
    z1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    #z1 = O.bodies[i].state.pos[1]
    list_displacement.append([getDisplacement(z1), z1])
    #list_displacement.append([getDisplacement(z1), O.bodies[i].state.pos[1]])
    
list_linearVel = []
for i in range(Ystep):
    v1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    #z1 = O.bodies[i].state.pos[1]
    list_linearVel.append([getLinearVel(v1), v1])
    #list_displacement.append([getDisplacement(z1), O.bodies[i].state.pos[1]])
    
list_contactNormalForce = []
for i in range(Ystep):
    n1 = dh*i + O.bodies[1].state.pos[1]
    list_contactNormalForce.append([getContactNormalForce(n1), n1])    
    
list_contactShearForce = []
for i in range(Ystep):
    s1 = dh*i + O.bodies[1].state.pos[1]
    list_contactShearForce.append([getContactShearForce(s1), s1])    
    
list_force = []
for i in range(Ystep):
    y1 = dh*i + O.bodies[1].state.pos[1]
    list_force.append([getForce(y1),y1])
    
list_torque = []
for i in range(Ystep):
    t1 = dh*i + O.bodies[1].state.pos[1]
    list_torque.append([getTorque(t1),t1])
    
list_momentum_x = []
for i in range(Ystep):
    mx1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_momentum_x.append([getMomentum_x(mx1), mx1])

list_momentum_y = []
for i in range(Ystep):
    my1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_momentum_y.append([getMomentum_y(my1), my1])
    
list_momentum_z = []
for i in range(Ystep):
    mz1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_momentum_z.append([getMomentum_z(mz1), mz1])
    
list_microstress = []
for i in range(Ystep):
    a1 = dh*i + O.bodies[1].state.pos[1] 
    list_microstress.append([getMicroStress(a1), a1])
    
list_porosity = []
for i in range(Ystep):
    p1 = dh*i + O.bodies[1].state.pos[1] 
    list_porosity.append([getPorosity(p1), p1])
    
list_angVel_x = []
for i in range(Ystep):
    avx1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_angVel_x.append([getAngVel_x(avx1), avx1])
    
list_angVel_y = []
for i in range(Ystep):
    avy1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_angVel_y.append([getAngVel_y(avy1), avy1])

list_angVel_z = []
for i in range(Ystep):
    avz1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_angVel_z.append([getAngVel_z(avz1), avz1])

list_angMomentum_x = []
for i in range(Ystep):
    amx1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_angMomentum_x.append([getAngMomentum_x(amx1), amx1])
    
list_angMomentum_y = []
for i in range(Ystep):
    amy1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_angMomentum_y.append([getAngMomentum_y(amy1), amy1])

list_angMomentum_z = []
for i in range(Ystep):
    amz1 = dh*i + O.bodies[1].state.pos[1] #Calculate the absolute z1 value 
    list_angMomentum_z.append([getAngMomentum_z(amz1), amz1])
    
np.savetxt(OUT+'-contactNormalForce',list_contactNormalForce, delimiter=',') 
np.savetxt(OUT+'-contactShearForce',list_contactShearForce, delimiter=',')    
np.savetxt(OUT+'-Torque',list_torque, delimiter=',')
np.savetxt(OUT+'-Pos_disp',list_displacement, delimiter=',')
np.savetxt(OUT+'-linearVel',list_linearVel, delimiter=',')
np.savetxt(OUT+'-Force',list_force, delimiter=',')
np.savetxt(OUT+'-Force_2D',list_force_2D, delimiter=',')
np.savetxt(OUT+'-Momentum_x',list_momentum_x, delimiter=',')
np.savetxt(OUT+'-Momentum_y',list_momentum_y, delimiter=',')
np.savetxt(OUT+'-Momentum_z',list_momentum_z, delimiter=',')
np.savetxt(OUT+'-microStress',list_microstress, delimiter=',')
np.savetxt(OUT+'-porosity',list_porosity, delimiter=',')
np.savetxt(OUT+'-angVel_x',list_angVel_x, delimiter=',')
np.savetxt(OUT+'-angVel_y',list_angVel_y, delimiter=',')
np.savetxt(OUT+'-angVel_z',list_angVel_z, delimiter=',')
np.savetxt(OUT+'-angMomentum_x',list_angMomentum_x, delimiter=',')
np.savetxt(OUT+'-angMomentum_y',list_angMomentum_y, delimiter=',')
np.savetxt(OUT+'-angMomentum_z',list_angMomentum_z, delimiter=',')

############################################################
### plot figure
plot.plots={'shear_strain':('Friction_coefficient'),'shear_strain2':('Friction_coefficient')}
plot.plot()

'''
### Micro-strain analysis
fixVel.dead = True
O.bodies[0].state.vel = (RATE_shear,0,0)
for b in boundary_top:
    b.state.blockedDOFs='x'
    b.state.vel=(RATE_shear,0,0)
    
TW=TesselationWrapper()
TW.triangulate()        
TW.computeVolumes()     
TW.volume(10)      

def particles_that_contain_in_the_box_below():
    ret = []
    for b in packing:
        if not isinstance(b.shape,Sphere): continue # skip non-spherical, not clear if yes or not
        ret.append(b)
    return ret

particles_0 = particles_that_contain_in_the_box_below()
with open("packing0_"+OUT+".txt","w") as f:
   for b in particles_0:
      x,y,z = b.state.pos
      f.write("{} {} {} {}\n".format(x,y,z,b.shape.radius))

#TW.setState(0)          
O.run(1000,True) 

particles_1 = particles_that_contain_in_the_box_below()
with open("packing1_"+OUT+".txt","w") as f:
   for p in particles_1:
      x,y,z = p.state.pos
      f.write("{} {} {} {}\n".format(x,y,z,p.shape.radius))
   
'''


#TW.defToVtkFromPositions( "packing0_"+OUT+".txt","packing1_"+OUT+".txt" , 'microStrain_'+OUT+'.vtk')
