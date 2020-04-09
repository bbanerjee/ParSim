# Visit 2.10.3 log file
ScriptVersion = "2.10.3"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
RestoreSession("/home/banerjee/ParSim/Matiti/runs/RigidHullSims/vflathull_local_upd.session", 0)
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.0394661, -0.942269, 0.332522)
View3DAtts.focus = (0, 0, 0.625)
View3DAtts.viewUp = (-0.0135797, 0.333257, 0.942738)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 1.11159
View3DAtts.nearPlane = -2.22317
View3DAtts.farPlane = 2.22317
View3DAtts.imagePan = (0.0563231, -0.00175693)
View3DAtts.imageZoom = 1.21
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 0.625)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state
Source("/home/banerjee/VisIt/build/bin/makemovie.py")
ToggleCameraViewMode()
#Set rendering attributes
RenderingAtts = RenderingAttributes()
RenderingAtts.antialiasing = 0
RenderingAtts.multiresolutionMode = 0
RenderingAtts.multiresolutionCellSize = 0.002
RenderingAtts.geometryRepresentation = RenderingAtts.Surfaces  # Surfaces, Wireframe, Points
RenderingAtts.displayListMode = RenderingAtts.Auto  # Never, Always, Auto
RenderingAtts.stereoRendering = 0
RenderingAtts.stereoType = RenderingAtts.CrystalEyes  # RedBlue, Interlaced, CrystalEyes, RedGreen
RenderingAtts.notifyForEachRender = 0
RenderingAtts.scalableActivationMode = RenderingAtts.Auto  # Never, Always, Auto
RenderingAtts.scalableAutoThreshold = 2000000
RenderingAtts.specularFlag = 0
RenderingAtts.specularCoeff = 0.6
RenderingAtts.specularPower = 10
RenderingAtts.specularColor = (255, 255, 255, 255)
RenderingAtts.doShadowing = 0
RenderingAtts.shadowStrength = 0.5
RenderingAtts.doDepthCueing = 0
RenderingAtts.depthCueingAutomatic = 1
RenderingAtts.startCuePoint = (-10, 0, 0)
RenderingAtts.endCuePoint = (10, 0, 0)
RenderingAtts.compressionActivationMode = RenderingAtts.Never  # Never, Always, Auto
RenderingAtts.colorTexturingFlag = 1
RenderingAtts.compactDomainsActivationMode = RenderingAtts.Never  # Never, Always, Auto
RenderingAtts.compactDomainsAutoThreshold = 256
SetRenderingAttributes(RenderingAtts)
#Set save window attributes
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "/home/banerjee/ParSim/Matiti/runs/RigidHullSims"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1109
SaveWindowAtts.height = 610
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
#Save the image
for time_cycle in range(1,21):
  SetTimeSliderState(time_cycle)
  SaveWindowAtts.fileName = "RigidCentrifuge_Full_FlatHull_NoWalls"+'%.4d'%time_cycle
  SetSaveWindowAttributes(SaveWindowAtts)
  SaveWindow()
