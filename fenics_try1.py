from fenics import *

def epsilon(u):
    return sym(grad(u))

def sigma(u, theta):
    E = 210e9  # Модуль Юнга, Па
    nu = 0.3  # Коэффициент Пуассона
    alpha = 1e-5  # Коэффициент термического расширения, 1/K
    theta_ref = 293  # Опорная температура, K
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    eps = epsilon(u)
    dim = u.geometric_dimension()  # Размерность
    return lambda_ * tr(eps) * Identity(dim) + 2 * mu * eps - alpha * (theta - theta_ref) * Identity(dim)

# Меш, созданный и подготовленный заранее
mesh = Mesh('./P_Fen/mesh.xml')
subdomains = MeshFunction('size_t', mesh, './P_Fen/mesh_physical_region.xml')
boundaries = MeshFunction('size_t', mesh, './P_Fen/mesh_facet_region.xml')

# Пространства функций
V = VectorFunctionSpace(mesh, 'P', 1)  # Смещения
T = FunctionSpace(mesh, 'P', 1)  # Температура

# Функции
u = Function(V)
theta = Function(T)
v = TestFunction(V)
psi = TestFunction(T)

# Вариационные формы
F_u = inner(sigma(u, theta), epsilon(v))*dx
F_theta = inner(grad(theta), grad(psi))*dx

J_u = derivative(F_u, u)
J_theta = derivative(F_theta, theta)

# Граничные условия
bc_u = DirichletBC(V, Constant((0, 0)), boundaries, 1)
bc_theta = DirichletBC(T, Constant(300), boundaries, 2)

# Решение
# solve(F_u == 0, u, bc_u)
# solve(F_theta == 0, theta, bc_theta)
solve(F_u == 0, u, bc_u, J=J_u)
solve(F_theta == 0, theta, bc_theta, J=J_theta)
# Сохранение результатов
file_u = File('displacement.pvd')
file_theta = File('temperature.pvd')
file_u << u
file_theta << theta