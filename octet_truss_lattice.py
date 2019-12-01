from paraview.simple import *
import numpy as np

reader = PVDReader(FileName="./test_cube.pvd")

data = servermanager.Fetch(reader)
cells = data.GetCells()
cells.InitTraversal()

# Initialize faces
data.InitializeFacesRepresentation(1)
faces = data.GetFaces()

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

    vert1 = np.array([xm, ym, zm])
    vert2 = np.array([xM, ym, zm])
    vert3 = np.array([xm, yM, zm])
    vert4 = np.array([xM, yM, zm])

    vert5 = np.array([xm, ym, zM])
    vert6 = np.array([xM, ym, zM])
    vert7 = np.array([xM, yM, zM])
    vert8 = np.array([xm, yM, zM])

    midpoint1 = (vert1 + vert2 + vert3 + vert4) / 4.0
    midpoint2 = (vert5 + vert6 + vert7 + vert8) / 4.0

    midpoint3 = (vert1 + vert3 + vert8 + vert5) / 4.0
    midpoint4 = (vert2 + vert6 + vert4 + vert7) / 4.0

    midpoint5 = (vert1 + vert2 + vert5 + vert6) / 4.0
    midpoint6 = (vert3 + vert4 + vert7 + vert8) / 4.0

    # Exterior trusses
    truss1 = np.array([vert1, vert4])
    truss2 = np.array([vert3, vert2])

    truss3 = np.array([vert8, vert6])
    truss4 = np.array([vert5, vert7])

    truss5 = np.array([vert6, vert4])
    truss6 = np.array([vert2, vert7])

    truss7 = np.array([vert1, vert8])
    truss8 = np.array([vert3, vert5])

    truss9 = np.array([vert1, vert6])
    truss10 = np.array([vert2, vert5])

    truss11 = np.array([vert8, vert4])
    truss12 = np.array([vert3, vert7])

    # Inner trusses
    itruss1 = np.array([midpoint2, midpoint5])
    itruss3 = np.array([midpoint2, midpoint6])
    itruss2 = np.array([midpoint1, midpoint5])
    itruss4 = np.array([midpoint1, midpoint6])

    itruss5 = np.array([midpoint1, midpoint3])
    itruss7 = np.array([midpoint1, midpoint4])
    itruss6 = np.array([midpoint2, midpoint3])
    itruss8 = np.array([midpoint2, midpoint4])

    itruss9 = np.array([midpoint4, midpoint6])
    itruss10 = np.array([midpoint4, midpoint5])
    itruss11 = np.array([midpoint3, midpoint5])
    itruss12 = np.array([midpoint3, midpoint6])

    for diag in [truss1, truss2, truss3, truss4, truss5, truss6, truss7, truss8, truss9, truss10, truss11, truss12, itruss1, itruss3, itruss2, itruss4, itruss5, itruss7, itruss6, itruss8, itruss9, itruss10, itruss11, itruss12]:
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
output.SetFileName('octet_truss.vtp')
output.Write()
