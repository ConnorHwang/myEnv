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
targetVariable = 'alpha.water'
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
ihfoamDisplay.ScaleFactor = 2.6899999618530277
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
ihfoam.CellArrays = ['U', 'alpha.water', 'alpha.water_0', 'epsilon', 'k', 'nut', 'p', 'p_rgh']

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
contour1Display.ScaleFactor = 2.6
contour1Display.SelectScaleArray = 'p'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'p'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'
contour1Display.GaussianRadius = 1.3
contour1Display.SetScaleArray = ['POINTS', 'alpha.water_0']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'alpha.water_0']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(ihfoam, renderView1)

# show color bar/color legend
#contour1Display.SetScalarBarVisibility(renderView1, True)
contour1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# change representation type
contour1Display.SetRepresentationType('Points')

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

ColorBy(contour1Display, None)

# change solid color
contour1Display.AmbientColor = [0.0, 0.0, 0.0]

#------------------------------------------------------------
# set active source
SetActiveSource(ihfoam)

# create a new 'Clip'
clip1 = Clip(Input=ihfoam)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'p']
clip1.Value = 4861.752126693726

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [13.449999809265137, 0.004999999888241291, 1.2000000476837158]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# show data in view
clip1Display = Show(clip1, renderView1)
# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'p']
clip1Display.LookupTable = pLUT
clip1Display.OSPRayScaleArray = 'p'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'U'
clip1Display.ScaleFactor = 2.6899999618530277
clip1Display.SelectScaleArray = 'p'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'p'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = pPWF
clip1Display.ScalarOpacityUnitDistance = 0.6244405664781962
clip1Display.GaussianRadius = 1.3449999809265138
clip1Display.SetScaleArray = ['POINTS', 'p']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'p']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(ihfoam, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, False)
#clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#### IF we want to change the variable
if targetVariable == 'alpha.water': targetVariable_ = 'alphawater'
else: targetVariable_ = targetVariable

# set scalar coloring
if targetVariable == 'U':
  ColorBy(clip1Display, ('CELLS', targetVariable, 'Magnitude'))
else:
  ColorBy(clip1Display, ('CELLS', targetVariable))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'nut'
myLUT = GetColorTransferFunction(targetVariable_)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
myLUT.ApplyPreset('Linear Blue (8_31f)', True)
#myLUT.ApplyPreset('Grayscale', True)
#myLUT.ApplyPreset('Rainbow Desaturated', True)
#myLUT.ApplyPreset('erdc_iceFire_L', True)
#myLUT.ApplyPreset('Cold and Hot', True)
#myLUT.ApplyPreset('Viridis (matplotlib)', True)
#myLUT.ApplyPreset('Cool to Warm (Extended)', True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [17.3, -45.97, 1.949]
renderView1.CameraFocalPoint = [17.3, 0.0, 1.949]
renderView1.CameraViewUp = [0.0, -0.008278146961434309, 0.9999657355544164]
renderView1.CameraParallelScale = 17.423

#### get color legend/bar for alphawaterLUT in view renderView1
myColorBar = GetScalarBar(myLUT, renderView1)

# change scalar bar placement
myColorBar.Orientation = 'Horizontal'
myColorBar.WindowLocation = 'AnyLocation'
myColorBar.Position = [0.3335954135954138, 0.6504874835309615]
myColorBar.ScalarBarLength = 0.3299999999999996
#myColorBar.
#-----------------------------------------------------

annotateTimeFilter1 = AnnotateTimeFilter(Input=clip1)
SetActiveSource(annotateTimeFilter1)
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)
#annotateTimeFilter1Display.FontFamily = 'Times'
annotateTimeFilter1Display.FontSize = 8
annotateTimeFilter1Display.WindowLocation = 'AnyLocation'
annotateTimeFilter1Display.Position = [0.75,0.65]

#-----------------------------------------------------

# current camera placement for renderView1
renderView1.CameraPosition = [17.3, -45.97, 1.949]
renderView1.CameraFocalPoint = [17.3, 0.0, 1.949]
renderView1.CameraViewUp = [0.0, -0.008278146961434309, 0.9999657355544164]
renderView1.CameraParallelScale = 17.423

for i in range(0,int(frames),int(step)):
  print(' >>> (%d/%d) Saving screenshots...' %(i,len(timeSeries)))
  saveFileName = mySavePath + '/' + targetVariable_ + str(i).zfill(4) + '.png'
  renderView1.ViewTime = timeSeries[i]
  clip1Display.RescaleTransferFunctionToDataRange(True, False)
  print(' >>> Current time: ', timeSeries[i])
  Render()
  # save screenshot
  SaveScreenshot(saveFileName, renderView1, ImageResolution=[1221, 759])

print(' >>> Aborting... ')
