
from fenics import *
import matplotlib.pyplot as plt

mesh = Mesh('./P_Fen/mesh.xml')
subdomains = MeshFunction('size_t', mesh, './P_Fen/mesh_physical_region.xml')
boundaries = MeshFunction('size_t', mesh, './P_Fen/mesh_facet_region.xml')

# Печать уникальных меток субдоменов
unique_subdomain_markers = set(subdomains.array())
print("Доступные метки субдоменов для dx: ", unique_subdomain_markers)

# Печать уникальных меток границ
unique_boundary_markers = set(boundaries.array())
print("Доступные метки границ для ds: ", unique_boundary_markers)

# Зададим константы для расчета
mu1 =     Constant(0.8*(10**(11)))
mu2 =     Constant(0.5*(10**(11)))
lmbda1 =  Constant(1.25*(10**(11)))
lmbda2 =  Constant(0.65*(10**(11)))
k1 =      Constant(200.0)
k2 =      Constant(100.0)
beta1 =   Constant(6.0*(10**(-6)))
beta2 =   Constant(3.0*(10**(-6)))
alphaAir= Constant(50.0)
TAir =    Constant(20.0)
THot =    Constant(150.0)
T0    =    Constant(0.)

#Задаем точки для граничных условий, условие Дирихле, они неподвижны  
lx = 0.02
ly = 0.012

V = VectorFunctionSpace(mesh, 'P', 2)
T = FunctionSpace(mesh, 'P', 2)

u = Function(V)
theta = Function(T)

v = TestFunction(V)
psi = TestFunction(T)

def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, mu, lmbda, theta, beta, T0):
    alpha = beta * (lmbda + 2. * mu / 3.)
    return 2 * mu * epsilon(u) + lmbda * tr(epsilon(u)) * Identity(u.geometric_dimension()) - alpha * (theta - T0) * Identity(u.geometric_dimension())

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Выбор материалов в зависимости от области
F_u = inner(sigma(u, mu1, lmbda1, theta, beta1, T0), epsilon(v)) * dx(1) \
    + inner(sigma(u, mu2, lmbda2, theta, beta2, T0), epsilon(v)) * dx(2)

F_theta = inner(k1 * grad(theta), grad(psi)) * dx(1) \
    + inner(k2 * grad(theta), grad(psi)) * dx(2) \
    + alphaAir * (theta - TAir) * psi * ds(3)

J_u = derivative(F_u, u)
J_theta = derivative(F_theta, theta)

# Boundary conditions
bc_u = [DirichletBC(V, Constant((0, ly)), boundaries, 2),DirichletBC(V, Constant((lx, ly)), boundaries, 2)]
bc_theta = [DirichletBC(T, Constant(300), boundaries,  2), DirichletBC(T, Constant(200), boundaries,  4)]
# 1 - лево 2 - верх  3 - право 4 - низ 


solver_parameters = {
    'newton_solver': {
        'maximum_iterations': 20,
        'absolute_tolerance': 1e-8,
        'relative_tolerance': 1e-7,
        'relaxation_parameter': 1.0
    }
}

solve(F_u == 0, u, bc_u, J=J_u, solver_parameters=solver_parameters)
solve(F_theta == 0, theta, bc_theta, J=J_theta,solver_parameters=solver_parameters)

file_u = File('displacement.pvd')
file_theta = File('temperature.pvd')
file_u << u
file_theta << theta

plot(u)
plt.title("Напряжение")
plt.show()

plot(theta)
plt.title("Температура")
plt.show()

# Преобразование результатов в массивы координат и значений
coordinates = mesh.coordinates()
u_values = u.compute_vertex_values(mesh)

# Создание сетки для визуализации
x, y = coordinates[:, 0], coordinates[:, 1]
u_x, u_y = u_values[0::2], u_values[1::2]  # предполагается, что u_values чередуют компоненты

plt.figure(figsize=(10, 6))
plt.quiver(x, y, u_x, u_y, scale=50)
plt.title('Распределение напряжений')
plt.xlabel('X координата')
plt.ylabel('Y координата')
plt.show()
