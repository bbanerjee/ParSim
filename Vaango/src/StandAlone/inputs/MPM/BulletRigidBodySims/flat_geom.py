# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=338.666656494141, 
    height=208.144454956055)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
Mdb()

#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.2)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=3)
s.setPrimaryObject(option=STANDALONE)

# Set up points
pt1 = (-0.04, 0.01)
pt2 = (-0.04, -0.02)
pt3 = (0.03, -0.02)
pt4 = (0.03, 0.01)
pt5 = (0.025, 0.01)
pt6 = (0.025, -0.015)
pt7 = (-0.035, -0.015)
pt8 = (-0.035, 0.01)

# Create spots at each point
s.Spot(pt1)
s.Spot(pt2)
s.Spot(pt3)
s.Spot(pt4)
s.Spot(pt5)
s.Spot(pt6)
s.Spot(pt7)
s.Spot(pt8)

# Create lines
s.Line(point1=pt1, point2 = pt2)
s.Line(point1=pt2, point2 = pt3)
s.Line(point1=pt3, point2 = pt4)
s.Line(point1=pt4, point2 = pt5)
s.Line(point1=pt5, point2 = pt6)
s.Line(point1=pt6, point2 = pt7)
s.Line(point1=pt7, point2 = pt8)
s.Line(point1=pt8, point2 = pt1)

# Extrude
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=0.1)
s.unsetPrimaryObject()

# Change viewport
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.198013, 
    farPlane=0.344668, width=0.112336, height=0.0723188, cameraPosition=(
    0.138306, 0.134014, 0.233855), cameraUpVector=(-0.560715, 0.617318, 
    -0.551831), cameraTarget=(-0.00377472, -0.00987317, 0.0536479))
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)

# Create independent instance
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=OFF)

# Set mesh contraols
partInstances =(a.instances['Part-1-1'], )
a.seedPartInstance(regions=partInstances, size=0.003, deviationFactor=0.1, 
    minSizeFactor=0.1)

c1 = a.instances['Part-1-1'].cells
#cells1 = c1.getSequenceFromMask(mask=('[#1 ]', ), )
pt13D = (pt1[0], pt1[1], 0.0)
cells1 = c1.findAt(pt13D,)
pickedRegions =(cells1, )
a.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)

# Generate mesh
a.generateMesh(regions=partInstances)

# Write input file
mdb.Job(name='flat_hull', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1)
mdb.jobs['flat_hull'].writeInput(consistencyChecking=OFF)

# Get the face area
f = p.faces
face = f.findAt(pt13D,)
print p.getArea((face,))
