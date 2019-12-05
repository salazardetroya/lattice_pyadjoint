from fenics import *
from fenics_adjoint import *
import itertools
import numpy as np

height = 10.0
width = 30.0
n_elem_side_y = 10
n_elem_side_x = 30
mesh = BoxMesh(Point(0, 0, 0), Point(width, height, height), n_elem_side_x, n_elem_side_y, n_elem_side_y)

W = VectorFunctionSpace(mesh, "CG", 1)

# Building the effective constituve tensor
vectors = list(itertools.product([1.0/sqrt(2.0), -1.0/sqrt(2.0)],repeat=2))
vectors_np = np.array(vectors)
zeros_column = np.zeros((len(vectors),1))
n_xy = np.c_[vectors_np[:,[0]], vectors_np[:,[1]], zeros_column]
n_yz = np.c_[zeros_column, vectors_np[:,[0]], vectors_np[:,[1]]]
n_xz = np.c_[vectors_np[:,[0]], zeros_column, vectors_np[:,[1]]]

normals_1 = [as_vector(n_xy[row_index,:]) for row_index in range(4)]
normals_2 = [as_vector(n_yz[row_index,:]) for row_index in range(4)]
normals_3 = [as_vector(n_xz[row_index,:]) for row_index in range(4)]
all_normals = [normals_1, normals_2, normals_3]

RHO = VectorFunctionSpace(mesh, "DG", 0)
A_list = Function(RHO, name='Control')
A_list.interpolate(Constant((3.0, 1.0, 1.0)))

length_truss = 1.0

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

from ufl import i, j, k, l
def sigma(r, A_list):
    Cijkl = [length_truss*A_section*all_normals[index][row_index][i]*all_normals[index][row_index][j]*all_normals[index][row_index][k]*all_normals[index][row_index][l] for index, A_section in enumerate(split(A_list)) for row_index in range(4) ]
    C = as_tensor(sum(Cijkl), (i,j,k,l))
    return C[i,j,k,l]*epsilon(r)[k,l]


# Define class marking Dirichlet boundary (x = 0 or x = 1)
class DirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 0.0) and on_boundary

# Define class marking Neumann boundary (y = 0 or y = 1)
class NeumanBoundary(SubDomain):
  def inside(self, x, on_boundary):
    #return x[1] > height - 0.001 and between(x[0], (0.0, width/5.0)) and on_boundary
    return near(x[0], width)

# Mark facets of the mesh
boundaries = MeshFunction('size_t', mesh, mesh.geometric_dimension() - 1)
NeumanBoundary().mark(boundaries, 2)
DirichletBoundary().mark(boundaries, 1)

# Define outer surface measure aware of Dirichlet and Neumann boundaries
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define test and trial functions
u = TestFunction(W)
v = TrialFunction(W)

a = inner(sigma(u, A_list), epsilon(v)[i,j])*dx
L = inner(Constant((0.0, -1.0, 0.0)), v)*ds(2)

bc = DirichletBC(W, Constant((0.0, 0.0, 0.0)), boundaries, 1)

u_sol = Function(W)
solve(a==L, u_sol, bcs=[bc])
File("u_sol.pvd") << u_sol
