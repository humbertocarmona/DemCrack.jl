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

def StreamLines(case, surface, z_level=20, y_inlet=50, on_ratio=50):
    mid_plane = Transform(Input=surface)
    mid_plane.Transform.Translate = [0.0, 0.0, z_level]
    mid_plane.Transform.Rotate = [0.0, 0.0, 0.0]
    mid_plane.Transform.Scale = [1.0, 1.0, 1.0]
    RenameSource('stream_level', mid_plane)

    slice = Slice(Input=mid_plane)
    slice.SliceType = 'Plane'
    slice.SliceOffsetValues = [0.0]
    slice.SliceType.Origin = [0.0, y_inlet, 0.0]
    slice.SliceType.Normal = [0.0, 1.0, 0.0]
    RenameSource('stream_font_line', slice)


    # just sort slide in a different way?
    re_order = ProgrammableFilter(Input=slice)
    re_order.Script = ""\
        "from paraview.vtk.numpy_interface import dataset_adapter as DA\n" \
        "import numpy as np;\n\n" \
        "pdi = self.GetInputDataObject(0, 0);\n" \
        "pdo = self.GetOutputDataObject(0);\n" \
        "# pdo.CopyAttributes(pdi);\n\n" \
        "old_pts = inputs[0].Points;\n" \
        "new_pts = old_pts[np.lexsort(np.fliplr(old_pts).T)]\n" \
        "arr = DA.numpyTovtkDataArray(new_pts, 'newpts');\n" \
        "pdo.GetPoints().SetData(arr);"
    RenameSource('ordered_line', re_order)

    mask = MaskPoints(Input=re_order)
    mask.OnRatio = on_ratio
    mask.MaximumNumberofPoints = 2000
    RenameSource('stream_point_seeds', mask)

    stream_lines = StreamTracerWithCustomSource(Input=case, SeedSource=mask)
    stream_lines.Vectors = ['POINTS', 'UN']
    stream_lines.MaximumStreamlineLength = 50000
    stream_lines.IntegrationDirection = 'BOTH'
    stream_lines.IntegratorType = 'Runge-Kutta 2'
    stream_lines.IntegrationStepUnit = 'Cell Length'
    stream_lines.InitialStepLength = 2.0
    stream_lines.MinimumStepLength = 0.1
    stream_lines.MaximumStepLength = 2.0
    RenameSource('stream_lines', stream_lines)

    return stream_lines, mask

y_clip = 30.0
y_cross = np.linspace(0,500-y_clip-1,5) + y_clip+1
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
steps = get_steps(case_name)
print(steps)
step = steps[0]
print('case={:}, step={:}'.format(case_name, step))
# open the case  "step"
file = './{:}/case.foam'.format(case_name)
full_case = OpenFOAMReader(FileName=file)
full_case.MeshRegions=['internalMesh', 'Bottom']
full_case.CellArrays=['U', 'p']
full_case.UpdatePipeline(step)
RenameSource('full_case', full_case)

renderView1 = GetActiveViewOrCreate('RenderView')


# first clip the head of the sample....
clip = Clip(Input=full_case)
clip.ClipType.Normal = [0., -1.0, 0.]
clip.ClipType.Origin = [250.0, y_clip, 15]
clip.UpdatePipeline(step)
RenameSource('clipped_head', clip)

case = Calculator(Input=clip)
case.ResultArrayName = 'UN'
case.Function = 'U/15.0'   # normalize velocity, should be Re = 100
RenameSource('normalize_U', case)
internal = Calculator(Input=case)
internal.ResultArrayName = 'z'
internal.Function = 'coordsZ'
RenameSource('internal', internal)

bottom = ExtractBlock(Input=internal, BlockIndices=[2])
bottom.UpdatePipeline(step)
RenameSource('bottom', bottom)
bottomDisplay = Show(bottom, renderView1)
ColorBy(bottomDisplay, 'z')

top = Transform(Input=bottom)
top.Transform.Translate = [0.0, 0.0, 40]
top.Transform.Rotate = [0.0, 0.0, 0.0]
top.Transform.Scale = [1.0, 1.0, 1.0]
RenameSource('top', top)

sl1, m1 = StreamLines(internal, bottom, z_level=cz,
                      y_inlet=cy, on_ratio=12)
slDisplay = Show(sl1, renderView1)
ColorBy(slDisplay, 'UN')
# bottom = MergeBlocks(Input=bottom)
# internal = ExtractBlock(Input=full_case, BlockIndices=[1])
# internal.UpdatePipeline(step)
# internal = MergeBlocks(internal)

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
uLUT.RescaleTransferFunction(0.0, 2)
uPWF.RescaleTransferFunction(0.0, 2)


zLUT = GetColorTransferFunction('z')
zPWF = GetOpacityTransferFunction('z')
zLUT.RescaleTransferFunction(-35.0, 35)
zLUT.ApplyPreset('copper_Matlab', True)
zPWF.RescaleTransferFunction(-35.0, 35)


box1 = Box()


box1.XLength = lx
box1.YLength = ly
box1.ZLength = lz
box1.Center = [cx, cy, cz]

box1Display = Show(box1, renderView1)
box1Display.Representation = 'Wireframe'


plane1 = Plane()
plane1.Origin = [0.5*lx-0.01, y_clip, cz - 0.5*lz]
plane1.Point1 = [0.5*lx-0.01, y_clip, cz + 0.5*lz]
plane1.Point2 = [0.5*lx-0.01, y_clip + ly, cz - 0.5*lz]
plane1.XResolution = 1
plane1.YResolution = 1
plane1Display = Show(plane1, renderView1)
plane1Display.Representation = 'Surface'


plane2 = Plane()
plane2.Origin = [0.5*lx-0.01, y_clip, cz - 0.5*lz]
plane2.Point1 = [0.5*lx-0.01, y_clip, cz + 0.5*lz]
plane2.Point2 = [0.5*lx-0.01, y_clip + ly, cz - 0.5*lz]
plane2.XResolution = 1
plane2.YResolution = 1
plane2Display = Show(plane1, renderView1)
plane2Display.Representation = 'Wireframe'

renderView1.Update()
