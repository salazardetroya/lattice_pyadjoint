from fenics import *
from fenics_adjoint import *
import itertools
import numpy as np

try:
    from pyadjoint import ipopt  # noqa: F401
except ImportError:
    print("""This example depends on IPOPT and Python ipopt bindings. \
  When compiling IPOPT, make sure to link against HSL, as it \
  is a necessity for practical problems.""")
    raise

height = 10.0
width = 30.0
n_elem_side_y = 20
n_elem_side_x = 60
mesh = BoxMesh.create([Point(0, 0, 0), Point(width, height, height)], [n_elem_side_x, n_elem_side_y, n_elem_side_y], CellType.Type.hexahedron)

output_dir = "./large"

# Building the effective constituve tensor
# Four struts per plane and their respective normals (not sure if they should
# be just two given that they are aligned
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
A_list.interpolate(Constant((0.0001, 0.0001, 0.0001)))

length_truss = width / n_elem_side_x * np.sqrt(2) / 2.0

# Material properties of alumina
E = 300e5 # N / cm^2
rho = 3.5 # Mg / m^3 or g / cm^3

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

from ufl import i, j, k, l
unit_cell_volume = length_truss**3 * 4.0 / np.sqrt(2.0)
def sigma(r, A_list):
    # Factor 2 because this plane of struts is repeated twice in the unit cell
    Cijkl = [Constant(E) * 2.0 * length_truss * A_section / Constant(unit_cell_volume) *
                    all_normals[index][row_index][i]*
                    all_normals[index][row_index][j]*
                    all_normals[index][row_index][k]*
                    all_normals[index][row_index][l]
                    for index, A_section in enumerate(split(A_list)) for row_index in range(4) ]
    C = as_tensor(sum(Cijkl), (i,j,k,l))
    return C[i,j,k,l]*epsilon(r)[k,l]


# Define class marking Dirichlet boundary (x = 0 or x = 1)
class DirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 0.0) and on_boundary

# Define class marking Neumann boundary (y = 0 or y = 1)
class NeumanBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], width) and between(x[1], (height*1.0/4.0, height*3.0/4.0)) \
            and between(x[2], (height*1.0/4.0, height*3.0/4.0))

# Mark facets of the mesh
boundaries = MeshFunction('size_t', mesh, mesh.geometric_dimension() - 1)
NeumanBoundary().mark(boundaries, 2)
DirichletBoundary().mark(boundaries, 1)

# Define outer surface measure aware of Dirichlet and Neumann boundaries
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define test and trial functions
W = VectorFunctionSpace(mesh, "CG", 1)
v = TestFunction(W)
u = TrialFunction(W)

a = inner(epsilon(v)[i,j], sigma(u, A_list))*dx
traction = Constant((0.0, -10.0, 0.0))
L = inner(traction, v)*ds(2)

bc = DirichletBC(W, Constant((0.0, 0.0, 0.0)), boundaries, 1)

u_sol = Function(W)
solve(a==L, u_sol, bcs=[bc])
File(output_dir + "/displacement.pvd") << u_sol
J = assemble(inner(u_sol, traction)*ds(2))

allctrls = File(output_dir + "/allcontrols.pvd")
rho_viz = Function(RHO, name="TrussArea")
def eval_derivative_cb(j, dj, vf):
    rho_viz.assign(vf)
    allctrls << rho_viz

vf = Control(A_list)
Jhat = ReducedFunctional(J, vf, derivative_cb_post=eval_derivative_cb)
Jhat.derivative()

# Bound constraints. This is the area.
lb = 0.0001
ub = 0.005


# Annotate to keep pyadjoint from recording this operation.
V = assemble(Constant(rho)*dx(domain=mesh), annotate=False) / (n_elem_side_x*n_elem_side_y**2)
delta = 0.05  # Volume ratio for the trusses
print("Total mass unit cell: {0:.5f}. Constrained mass unit cell: {1:.5f}".format(V, V*delta))
trusses_per_plane = 4 * 2
mass_per_variable = length_truss * trusses_per_plane * rho
volume_constraint = UFLInequalityConstraint((V*delta - \
                    mass_per_variable*(A_list[0] + A_list[1] + A_list[2]))*dx, vf)

problem = MinimizationProblem(Jhat, bounds=(lb, ub), constraints=volume_constraint)
parameters = {'maximum_iterations': 400}


solver = IPOPTSolver(problem, parameters=parameters)
rho_opt = solver.solve()
