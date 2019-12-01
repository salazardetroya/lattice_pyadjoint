from fenics import *

mesh = UnitCubeMesh.create(10, 10, 10, CellType.Type.hexahedron)

V = FunctionSpace(mesh, 'CG', 1)

u = Function(V)
u.interpolate(Constant(1.0))
File("test_cube.pvd") << u

