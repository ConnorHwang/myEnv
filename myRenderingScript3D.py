#### Hanul Connor Hwang
#### Paraview rendering script

import sys
import os
from paraview.simple import *

casePath = os.getcwd()
caseName = casePath.replace('/', 'X', casePath.count('/')-1)
caseName = caseName[(caseName.find('/')+1):]
casePath = casePath + '/'

#targetVariable = 'epsilon'
#targetVariable = 'alpha.water'
#targetVariable = 'nut'
#targetVariable = 'U'

myFileLoad = casePath + 'wave.foam'
mySavePath = './images'
if not os.path.isdir(mySavePath):
  print(' >>> Making images directory...')
  os.makedirs(mySavePath)

#if ((len(sys.argv) < 2) or (len(sys.argv) > 2)):
if not (len(sys.argv) == 3):
  print(' >>> WARNING: no input arguments <frame> <step>')
  frames = 999999
  step = 1
else:
  frames = int(sys.argv[1])
  step = int(sys.argv[2])
print(' >>> First input argument: %d, and second: %d' %(frames,step))

#try:
#  (sys.argv[1] or sys.argv[2])
#except noInputArg:
#  print(' >>> WARNING: no input arguments <frame> <step>')
#else: 
#  frames = int(sys.argv[1])
#  step = float(sys.argv[2])
#  print(' >>> First input argument: %d, and second: %d' %(frames,step))

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
ihfoam = OpenFOAMReader(FileName=myFileLoad)
ihfoam.MeshRegions = ['internalMesh']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1221, 759]

# Properties modified on renderView1
renderView1.UseTexturedBackground = 0

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 0

# Properties modified on ihfoam
ihfoam.SkipZeroTime = 0
ihfoam.CaseType = 'Decomposed Case'
ihfoam.ScalarSize = '32-bit (SP)'

# get animation scene
animationScene1 = GetAnimationScene()
#AnimationScene1.EndTime = (frames-1)*step
#AnimationScene1.PlayMode = 'Snap To timeSteps'

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Get the length of the time list
timeKeeper1 = GetTimeKeeper()
timeSeries = timeKeeper1.TimestepValues
#reader = GetActiveSource()
#timeSeries = reader.TimestepValues
#timeSeries = renderView1.TimestepValues
print(' >>> Number of time frames: ', len(timeSeries))
if (frames >= len(timeSeries)):
  print(' >>> Requested time frame is out of bounds!')
  print('     Setting the available last time frame')
  frames = len(timeSeries)
  print(' >>> End of the frame: %d with step interval %d' %(frames,step))
print(' >>> Time series from t = ', timeSeries[0], ', to t = ', timeSeries[-1])

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# show data in view
ihfoamDisplay = Show(ihfoam, renderView1)
# trace defaults for the display properties.
ihfoamDisplay.Representation = 'Surface'
ihfoamDisplay.ColorArrayName = ['POINTS', 'p']
ihfoamDisplay.LookupTable = pLUT
ihfoamDisplay.OSPRayScaleArray = 'p'
ihfoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ihfoamDisplay.SelectOrientationVectors = 'U'
ihfoamDisplay.ScaleFactor = 3.4626688
ihfoamDisplay.SelectScaleArray = 'p'
ihfoamDisplay.GlyphType = 'Arrow'
ihfoamDisplay.GlyphTableIndexArray = 'p'
ihfoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
ihfoamDisplay.PolarAxes = 'PolarAxesRepresentation'
ihfoamDisplay.ScalarOpacityFunction = pPWF
ihfoamDisplay.ScalarOpacityUnitDistance = 0.6244405985833273
ihfoamDisplay.GaussianRadius = 1.3449999809265138
ihfoamDisplay.SetScaleArray = ['POINTS', 'p']
ihfoamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ihfoamDisplay.OpacityArray = ['POINTS', 'p']
ihfoamDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
ihfoamDisplay.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
#renderView1.ResetCamera()
# create a new 'Contour'
contour1 = Contour(Input=ihfoam)
contour1.ContourBy = ['POINTS', 'p']
contour1.Isosurfaces = [4861.752126693726]
contour1.PointMergeMethod = 'Uniform Binning'

# Properties modified on ihfoam
ihfoam.CellArrays = ['U', 'alpha.water', 'k', 'nut', 'omega', 'p', 'p_rgh']

# Properties modified on contour1
contour1.ContourBy = ['POINTS', 'alpha.water']
contour1.Isosurfaces = [0.5]

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['CELLS', 'p']
contour1Display.LookupTable = pLUT
contour1Display.OSPRayScaleArray = 'p'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'U'
contour1Display.ScaleFactor = 3.082184063769546
contour1Display.SelectScaleArray = 'p'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'p'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'
contour1Display.GaussianRadius = 1.541092
contour1Display.SetScaleArray = ['POINTS', 'alpha.water_0']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'alpha.water_0']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(ihfoam, renderView1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# change representation type
#contour1Display.SetRepresentationType('Points')

# Hide the scalar bar for this color map if no visible data is colored by it.
#HideScalarBarIfNotNeeded(pLUT, renderView1)

ColorBy(contour1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
contour1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# change solid color
contour1Display.DiffuseColor = [0.0, 0.6666666666666666, 1.0]

#------------------------------------------------------------
# set active source
SetActiveSource(None)

# create a new 'STL Reader'
soars3DZstl = STLReader(FileNames=['/p/lustre2/hwang8/niche/soars3D/constant/triSurface/soars3DZ.stl'])

# get color transfer function/color map for 'STLSolidLabeling'
sTLSolidLabelingLUT = GetColorTransferFunction('STLSolidLabeling')

# show data in view
soars3DZstlDisplay = Show(soars3DZstl, renderView1)
# trace defaults for the display properties.
soars3DZstlDisplay.Representation = 'Surface'
soars3DZstlDisplay.ColorArrayName = ['CELLS', 'STLSolidLabeling']
soars3DZstlDisplay.LookupTable = sTLSolidLabelingLUT
soars3DZstlDisplay.OSPRayScaleArray = 'STLSolidLabeling'
soars3DZstlDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
soars3DZstlDisplay.SelectOrientationVectors = 'None'
soars3DZstlDisplay.ScaleFactor = 3.4626715057861475
soars3DZstlDisplay.SelectScaleArray = 'STLSolidLabeling'
soars3DZstlDisplay.GlyphType = 'Arrow'
soars3DZstlDisplay.GlyphTableIndexArray = 'STLSolidLabeling'
soars3DZstlDisplay.DataAxesGrid = 'GridAxesRepresentation'
soars3DZstlDisplay.PolarAxes = 'PolarAxesRepresentation'
soars3DZstlDisplay.GaussianRadius = 1.7313357528930737
soars3DZstlDisplay.SetScaleArray = [None, '']
soars3DZstlDisplay.ScaleTransferFunction = 'PiecewiseFunction'
soars3DZstlDisplay.OpacityArray = [None, '']
soars3DZstlDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
soars3DZstlDisplay.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(sTLSolidLabelingLUT, renderView1)

# Properties modified on soars3DZstlDisplay
soars3DZstlDisplay.Opacity = 0.1

# reset view to fit data
#renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# 3D view: current camera placement for renderView1
renderView1.CameraPosition = [41.600638971658526, -17.094727317379085, 13.312764661401507]
renderView1.CameraFocalPoint = [20.645708638015922, 2.867557096604591, 0.9328125914846953]
renderView1.CameraViewUp = [-0.29985247945826393, 0.25529017497115625, 0.919192807372821]
renderView1.CameraParallelScale = 17.46405062816704

## 2D view: current camera placement for renderView1
#renderView1.CameraPosition = [17.3, -45.97, 1.949]
#renderView1.CameraFocalPoint = [17.3, 0.0, 1.949]
#renderView1.CameraViewUp = [0.0, -0.008278146961434309, 0.9999657355544164]
#renderView1.CameraParallelScale = 17.423

#-----------------------------------------------------

annotateTimeFilter1 = AnnotateTimeFilter(Input=contour1)
SetActiveSource(annotateTimeFilter1)
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)
#annotateTimeFilter1Display.FontFamily = 'Times'
annotateTimeFilter1Display.FontSize = 8
annotateTimeFilter1Display.WindowLocation = 'AnyLocation'
annotateTimeFilter1Display.Position = [0.75,0.70]

#-----------------------------------------------------

# current camera placement for renderView1
renderView1.CameraPosition = [41.600638971658526, -17.094727317379085, 13.312764661401507]
renderView1.CameraFocalPoint = [20.645708638015922, 2.867557096604591, 0.9328125914846953]
renderView1.CameraViewUp = [-0.29985247945826393, 0.25529017497115625, 0.919192807372821]
renderView1.CameraParallelScale = 17.46405062816704

for i in range(0,int(frames),int(step)):
  print(' >>> (%d/%d) Saving screenshots...' %(i,len(timeSeries)))
  saveFileName = mySavePath + '/' + '3dvof' + str(i).zfill(4) + '.png'
  renderView1.ViewTime = timeSeries[i]
  contour1Display.RescaleTransferFunctionToDataRange(True, False)
  soars3DZstlDisplay.RescaleTransferFunctionToDataRange(True, False)
  print(' >>> Current time: ', timeSeries[i])
  Render()
  # save screenshot
  SaveScreenshot(saveFileName, renderView1, ImageResolution=[1221, 759])

print(' >>> Aborting... ')
