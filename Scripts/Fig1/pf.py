#from paraview.vtk.numpy_interface import dataset_adapter as DA
import numpy as np
import pandas as pd

def flatten(input, output):
    # Copy the cells etc.
    output.ShallowCopy(input)
    newPoints = vtk.vtkPoints()
    numPoints = input.GetNumberOfPoints()
    xx = np.zeros(numPoints)
    yy = np.zeros(numPoints)
    zz = np.zeros(numPoints)

    for i in range(0, numPoints):
        coord = input.GetPoint(i)
        x, y, z = coord[:3]
        xx[i] = x
        yy[i] = y
        zz[i] = z

    df = pd.DataFrame({'x':xx, 'y':yy, 'z':zz})
    df = df.sort_values(by=['x'])
    for i in range(0, numPoints):
        x = df.x[i]
        y = df.y[i]
        z = df.z[i]
        newPoints.InsertPoint(i, x, y, z)

    output.SetPoints(newPoints)


input = self.GetInputDataObject(0, 0);
output = self.GetOutputDataObject(0);

if input.IsA("vtkMultiBlockDataSet"):
    output.CopyStructure(input)
    iter = input.NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    while not iter.IsDoneWithTraversal():
        curInput = iter.GetCurrentDataObject()
        curOutput = curInput.NewInstance()
        curOutput.UnRegister(None)
        output.SetDataSet(iter, curOutput)
        flatten(curInput, curOutput)
        iter.GoToNextItem();
else:
  flatten(input, output)
