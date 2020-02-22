    #### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import numpy 
import subprocess, re, os
import tarfile

vtu_files = []
results_tar = "../results/seed20.tar"

def extract_vtu(file, path=r'./', step=50):
    tar = tarfile.open(results_tar)
    n_files = len(tar.getmembers())
    files = []
    for i in numpy.arange(0,n_files-1, step):
        name = tar.getmembers()[i].name
        vtu_files.append(name)
    for i, filename in enumerate(vtu_files):
        tar.extract(filename, path=path)
        files.append('{:}/{:}'.format(path, filename))
    return files
def find_vtu(path=r'./', pattern='.+\.vtu'):
    findCommand ='find {:}  -regex "{:}" |sort'.format(path,pattern)
    files = subprocess.check_output(findCommand, shell=True)
    files = files.decode("utf-8")
    files=files.strip()
    files= re.split('\n',files)
    return files

def remove_tmp_files(tmp):
    for f in os.listdir(tmp):
        os.remove(f)
    os.rmdir(tmp)

vtu_files = extract_vtu(results_tar, path='output')
vtu_files = find_vtu(path='output')


# create a new 'XML Unstructured Grid Reader'
snap_0 = XMLUnstructuredGridReader(FileName=vtu_files)
#snap_0.PointArrayStatus = ['radius', 'vel']

snap_0.CellArrayStatus = ['bondVector', 'connected', 'bond_elo', 'bond_p', 'bond_q', 'bond_len']
snap_0.PointArrayStatus = ['v', 'a', 'diam', 's1', 's2', 's3', 'elmtSet', 'pbc']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1256, 518]

# show data in view
snap_0Display = Show(snap_0, renderView1)

# trace defaults for the display properties.
snap_0Display.Representation = 'Surface'
snap_0Display.ColorArrayName = [None, '']
snap_0Display.OSPRayScaleArray = 'a'
snap_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
snap_0Display.SelectOrientationVectors = 'None'
snap_0Display.ScaleFactor = 5.685805734884344
snap_0Display.SelectScaleArray = 'None'
snap_0Display.GlyphType = 'Arrow'
snap_0Display.GlyphTableIndexArray = 'None'
snap_0Display.GaussianRadius = 0.28429028674421714
snap_0Display.SetScaleArray = ['POINTS', 'a']
snap_0Display.ScaleTransferFunction = 'PiecewiseFunction'
snap_0Display.OpacityArray = ['POINTS', 'a']
snap_0Display.OpacityTransferFunction = 'PiecewiseFunction'
snap_0Display.DataAxesGrid = 'GridAxesRepresentation'
snap_0Display.SelectionCellLabelFontFile = ''
snap_0Display.SelectionPointLabelFontFile = ''
snap_0Display.PolarAxes = 'PolarAxesRepresentation'
snap_0Display.ScalarOpacityUnitDistance = 1.1783244563324338

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
snap_0Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.5, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
snap_0Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.5, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
snap_0Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.5, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
snap_0Display.DataAxesGrid.XTitleFontFile = ''
snap_0Display.DataAxesGrid.YTitleFontFile = ''
snap_0Display.DataAxesGrid.ZTitleFontFile = ''
snap_0Display.DataAxesGrid.XLabelFontFile = ''
snap_0Display.DataAxesGrid.YLabelFontFile = ''
snap_0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
snap_0Display.PolarAxes.PolarAxisTitleFontFile = ''
snap_0Display.PolarAxes.PolarAxisLabelFontFile = ''
snap_0Display.PolarAxes.LastRadialAxisTextFontFile = ''
snap_0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
snap_0Display.SetRepresentationType('Point Gaussian')

# set scalar coloring
ColorBy(snap_0Display, ('POINTS', 'elmtSet'))

# rescale color and/or opacity maps used to include current data range
snap_0Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
snap_0Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'elmtSet'
elmtSetLUT = GetColorTransferFunction('elmtSet')

# get opacity transfer function/opacity map for 'elmtSet'
elmtSetPWF = GetOpacityTransferFunction('elmtSet')

# Properties modified on snap_0Display
snap_0Display.GaussianRadius = 0.5

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [95.83240241020361, 110.21964309086397, 185.83177053046484]
renderView1.CameraFocalPoint = [28.42922007155722, 28.389059998095014, 28.429323986711093]
renderView1.CameraViewUp = [-0.16071432116105802, 0.9022151485963638, -0.40022335341279336]
renderView1.CameraParallelScale = 49.11766582493702

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).