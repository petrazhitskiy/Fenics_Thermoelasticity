from fenics import *

def epsilon(u):
    return sym(grad(u))

def sigma(u, theta):
    E = 210e9  # Young's modulus, Pa
    nu = 0.3   # Poisson's ratio
    alpha = 1e-5  # Coefficient of thermal expansion, 1/K
    theta_ref = 293  # Reference temperature, K
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    dim = u.geometric_dimension()
    return lambda_ * tr(epsilon(u)) * Identity(dim) + 2 * mu * epsilon(u) - alpha * (theta - theta_ref) * Identity(dim)

# Load mesh and subdomains
mesh = Mesh('./P_Fen/mesh.xml')
subdomains = MeshFunction('size_t', mesh, './P_Fen/mesh_physical_region.xml')
boundaries = MeshFunction('size_t', mesh, './P_Fen/mesh_facet_region.xml')

# Function spaces
V = VectorFunctionSpace(mesh, 'P', 2)
T = FunctionSpace(mesh, 'P', 2)

# Functions
u = Function(V)
theta = Function(T)
v = TestFunction(V)
psi = TestFunction(T)

# Variational forms
F_u = inner(sigma(u, theta), epsilon(v))*dx
F_theta = inner(grad(theta), grad(psi))*dx

# Jacobian
J_u = derivative(F_u, u)
J_theta = derivative(F_theta, theta)

# Boundary conditions
bc_u = DirichletBC(V, Constant((0, 0)), boundaries, 1)
bc_theta = DirichletBC(T, Constant(300), boundaries, 2)

# Solve

solve(F_u == 0, u, bc_u, J=J_u)
solve(F_theta == 0, theta, bc_theta, J=J_theta)

# Save results
file_u = File('displacement.pvd')
file_theta = File('temperature.pvd')
file_u << u
file_theta << theta

