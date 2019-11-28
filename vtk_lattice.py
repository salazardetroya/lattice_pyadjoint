from paraview.simple import *
import numpy as np

reader = PVDReader(FileName="./test_cube.pvd")

data = servermanager.Fetch(reader)
cells = data.GetCells()
cells.InitTraversal()

list_cylinders = []

import vtk
# line source
line = vtk.vtkLineSource()

# tube source
tube = vtk.vtkTubeFilter()
tube.SetInputConnection(line.GetOutputPort())
tube.SetCapping(True)
tube.SetNumberOfSides(6)

model = vtk.vtkAppendPolyData()
i = 0
for i in range(cells.GetNumberOfCells()):
    cell = data.GetCell(i)
    [xm, xM, ym, yM, zm, zM] = cell.GetBounds()

    diag1 = np.array([[xm, ym, zm], [xM, yM, zM]])
    diag2 = np.array([[xm, ym, zM], [xM, yM, zm]])
    diag3 = np.array([[xm, yM, zm], [xM, ym, zM]])
    diag4 = np.array([[xM, ym, zm], [xm, yM, zM]])

    for diag in [diag1, diag2, diag3, diag4]:
        i += 1
        input1 = vtk.vtkPolyData()
        #draw line
        line.SetPoint1(diag[0])
        line.SetPoint2(diag[1])
        line.Update()

        #set tube radius
        tube.SetRadius(0.001)
        tube.Update()

        input1.ShallowCopy(tube.GetOutput())

        model.AddInputData(input1)

model.Update()
print("Number of cylinders {}".format(i))
output = vtk.vtkXMLPolyDataWriter()
output.SetInputData(model.GetOutput())
output.SetFileName('output.vtp')
output.Write()
