# import argparse
# import re
# import os
import sys
from paraview.simple import *
import os
import re
import numpy as np


def get_steps(case_name):
    sdirs = []
    for root, subdirs, _ in os.walk(case_name):
        if root == case_name:
            sdirs = subdirs
            break
    sdirs = np.unique(np.array(sdirs))
    dirs = np.array([])
    for s in sdirs:
        if (re.match('^\d\d', s)):
            dirs = np.append(dirs, s)
    if len(dirs) > 0:
        dirs = np.sort(dirs.astype(int))
    else:
        dirs = 'None'
    return dirs

def StreamLines(case, surface, z_level=15, y_inlet=250):
    mid_plane = Transform(Input=surface)
    mid_plane.Transform.Translate = [0.0, 0.0, z_level]
    mid_plane.Transform.Rotate = [0.0, 0.0, 0.0]
    mid_plane.Transform.Scale = [1.0, 1.0, 1.0]
    RenameSource('stream_level', mid_plane)
    Hide3DWidgets(proxy=mid_plane.Transform)

    slice = Slice(Input=mid_plane)
    slice.SliceType = 'Plane'
    slice.SliceOffsetValues = [0.0]
    slice.SliceType.Origin = [0.0, y_inlet, 0.0]
    slice.SliceType.Normal = [0.0, 1.0, 0.0]
    RenameSource('stream_font_line', slice)
    Hide3DWidgets(proxy=slice.SliceType)

    seed_points = ProgrammableFilter(Input=slice)
    seed_points.Script = ""\
        "import numpy as np\n"\
        "import pandas as pd\n"\
        "\n"\
        "def flatten(input, output):\n"\
        "    #output.ShallowCopy(input)\n"\
        "    newPoints = vtk.vtkPoints()\n"\
        "    numPoints = input.GetNumberOfPoints()\n"\
        "    xx = np.zeros(numPoints)\n"\
        "    yy = np.zeros(numPoints)\n"\
        "    zz = np.zeros(numPoints)\n"\
        "\n"\
        "    for i in range(0, numPoints):\n"\
        "        coord = input.GetPoint(i)\n"\
        "        x, y, z = coord[:3]\n"\
        "        xx[i] = x\n"\
        "        yy[i] = y\n"\
        "        zz[i] = z\n"\
        "\n"\
        "    df = pd.DataFrame({'x':xx, 'y':yy, 'z':zz})\n"\
        "    df = df.sort_values(by=['x'])\n"\
        "    j = 0\n"\
        "    for i in range(0, numPoints):\n"\
        "        x = df.x[i]\n"\
        "        y = df.y[i]\n"\
        "        z = df.z[i]\n"\
        "        if x < 500 and np.abs(x % 8) < 1:\n"\
        "           newPoints.InsertPoint(j, x, y, z)\n"\
        "           j = j+1\n"\
        "\n"\
        "    output.SetPoints(newPoints)\n"\
        "\n"\
        "input = self.GetInputDataObject(0, 0);\n"\
        "output = self.GetOutputDataObject(0);\n"\
        "\n"\
        "if input.IsA('vtkMultiBlockDataSet'):\n"\
        "    output.CopyStructure(input)\n"\
        "    iter = input.NewIterator()\n"\
        "    iter.UnRegister(None)\n"\
        "    iter.InitTraversal()\n"\
        "    while not iter.IsDoneWithTraversal():\n"\
        "        curInput = iter.GetCurrentDataObject()\n"\
        "        curOutput = curInput.NewInstance()\n"\
        "        curOutput.UnRegister(None)\n"\
        "        output.SetDataSet(iter, curOutput)\n"\
        "        flatten(curInput, curOutput)\n"\
        "        iter.GoToNextItem();\n"\
        "else:\n"\
        "  flatten(input, output)\n"
    RenameSource('seed_points', seed_points)

    stream_lines = StreamTracerWithCustomSource(Input=case, SeedSource=seed_points)
    stream_lines.Vectors = ['POINTS', 'UN']
    stream_lines.MaximumStreamlineLength = 50000
    stream_lines.IntegrationDirection = 'BOTH'
    stream_lines.IntegratorType = 'Runge-Kutta 2'
    stream_lines.IntegrationStepUnit = 'Cell Length'
    stream_lines.InitialStepLength = 2.0
    stream_lines.MinimumStepLength = 0.1
    stream_lines.MaximumStepLength = 2.0
    RenameSource('stream_lines', stream_lines)

    return stream_lines

y_clip = 30.0
y_cross = np.linspace(0,500-y_clip-1, 4) + y_clip+1
x_cross = [1, 250, 500]
hurst = 0.8
lx = 500
ly = 500 - y_clip
lz = 113

cx = 0.5*lx
cy = 0.5*ly + y_clip
cz = 16.7
# z: [-34.4469, 67.8187]

case_name = 'crack_{:0.2f}_seed_002_d_40.00'.format(hurst)
print('case={:}'.format(case_name))
file = './{:}/case.foam'.format(case_name)
full_case = OpenFOAMReader(FileName=file)
full_case.MeshRegions=['internalMesh', 'Bottom']
full_case.CellArrays=['U', 'p']
RenameSource('full_case', full_case)

renderView1 = GetActiveViewOrCreate('RenderView')

# first clip the head of the sample....
clip = Clip(Input=full_case)
clip.ClipType.Normal = [0., -1.0, 0.]
clip.ClipType.Origin = [250.0, y_clip, 15]
RenameSource('clipped_head', clip)

#compute normalized vel and z value
case = Calculator(Input=clip)
case.ResultArrayName = 'UN'
case.Function = 'U/5.0'   # normalize velocity, should be Re = 100
RenameSource('normalize_U', case)
internal = Calculator(Input=case)
internal.ResultArrayName = 'z'
internal.Function = 'coordsZ'
RenameSource('internal', internal)

bottom = ExtractBlock(Input=internal, BlockIndices=[2])
RenameSource('bottom', bottom)
bottomDisplay = Show(bottom, renderView1)
bottomDisplay.Opacity = 0.75
ColorBy(bottomDisplay, 'z')

top = Transform(Input=bottom)
top.Transform.Translate = [0.0, 0.0, 40]
top.Transform.Rotate = [0.0, 0.0, 0.0]
top.Transform.Scale = [1.0, 1.0, 1.0]
RenameSource('top', top)
Hide3DWidgets(proxy=top.Transform)

stream = StreamLines(internal, bottom, z_level=cz, y_inlet=cy)
streamDisplay = Show(stream, renderView1)
ColorBy(streamDisplay, 'UN')

# cross sections in y-planes
for i, yi in enumerate(y_cross):
    s = Slice(Input=internal)
    s.SliceOffsetValues = [0.0]
    s.SliceType.Origin = [250.0, yi, 15.0]
    s.SliceType.Normal = [0.0, 1.0, 0.0]
    s.SliceType = 'Plane'
    RenameSource('y{:}'.format(yi), s)
    sDisplay = Show(s, renderView1)
    ColorBy(sDisplay, 'UN')
    Hide3DWidgets(proxy=s.SliceType)

    if i == 0:
        integral = IntegrateVariables(Input=s)
        integral.DivideCellDataByVolume = 1

    e = FeatureEdges(Input=s)
    e.BoundaryEdges = 1
    e.FeatureEdges = 0
    e.NonManifoldEdges = 0
    RenameSource('edges_y{:}'.format(yi), e)
    eDisplay = Show(e, renderView1)
    ColorBy(eDisplay, None)
    eDisplay.AmbientColor = [0.0, 0.0, 0.0]
    eDisplay.DiffuseColor = [0.0, 0.0, 0.0]

# could add a set of cross sections perpendicular to the flow
for xi in x_cross:
    s = Slice(Input=internal)
    s.SliceOffsetValues = [0.0]
    s.SliceType.Origin = [xi, 250, 15.0]
    s.SliceType.Normal = [1.0, 0.0, 0.0]
    s.SliceType = 'Plane'
    RenameSource('x{:}'.format(xi), s)
    sDisplay = Show(s, renderView1)
    ColorBy(sDisplay, 'UN')
    Hide3DWidgets(proxy=s.SliceType)


    e = FeatureEdges(Input=s)
    e.BoundaryEdges = 1
    e.FeatureEdges = 0
    e.NonManifoldEdges = 0
    RenameSource('edges_x{:}'.format(xi), e)
    eDisplay = Show(e, renderView1)
    ColorBy(eDisplay, None)
    eDisplay.AmbientColor = [0.0, 0.0, 0.0]
    eDisplay.DiffuseColor = [0.0, 0.0, 0.0]

uLUT = GetColorTransferFunction('U')
uPWF = GetOpacityTransferFunction('U')
uLUT.RescaleTransferFunction(0.0, 20)
uPWF.RescaleTransferFunction(0.0, 20)

unLUT = GetColorTransferFunction('UN')
unPWF = GetOpacityTransferFunction('UN')
unLUT.RescaleTransferFunction(0.0, 2.0)
unPWF.RescaleTransferFunction(0.0, 2.0)


zLUT = GetColorTransferFunction('z')
zPWF = GetOpacityTransferFunction('z')
zLUT.RescaleTransferFunction(-35.0, 35)
zLUT.ApplyPreset('Grayscale', True)
zPWF.RescaleTransferFunction(-35.0, 35)

box1 = Box()
box1.XLength = lx
box1.YLength = ly
box1.ZLength = lz
box1.Center = [cx, cy, cz]
box1Display = Show(box1, renderView1)
box1Display.Representation = 'Wireframe'
box1Display.LineWidth = 0.25

plane1 = Plane()
plane1.Origin = [0.5*lx-0.01, y_clip, cz - 0.5*lz]
plane1.Point1 = [0.5*lx-0.01, y_clip, cz + 0.5*lz]
plane1.Point2 = [0.5*lx-0.01, y_clip + ly, cz - 0.5*lz]
plane1.XResolution = 1
plane1.YResolution = 1
plane1Display = Show(plane1, renderView1)
plane1Display.Representation = 'Surface'
plane1Display.Opacity = 0.25
plane1Display.AmbientColor = [0.33, 0.67, 0.5]
plane1Display.DiffuseColor = [0.33, 0.67, 0.5]


plane2 = Plane()
plane2.Origin = [0.5*lx-0.01, y_clip, cz - 0.5*lz]
plane2.Point1 = [0.5*lx-0.01, y_clip, cz + 0.5*lz]
plane2.Point2 = [0.5*lx-0.01, y_clip + ly, cz - 0.5*lz]
plane2.XResolution = 1
plane2.YResolution = 1
plane2Display = Show(plane2, renderView1)
plane2Display.Representation = 'Wireframe'
plane2Display.LineWidth = 0.25


LoadPalette(paletteName='PrintBackground')

renderView1.Update()
