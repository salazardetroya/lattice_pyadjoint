from fenics import *

mesh = UnitCubeMesh.create(20, 20, 20, CellType.Type.hexahedron)

V = FunctionSpace(mesh, 'CG', 1)

u = Function(V)
u.interpolate(Constant(1.0))
File("test_cube.pvd") << u

