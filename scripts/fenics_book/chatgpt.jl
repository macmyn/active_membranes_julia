using FEniCS
# using PyPlot

# Define the mesh and function space
mesh = UnitIntervalMesh(100)
V = FunctionSpace(mesh, "P", 1)

# Define the boundary condition
u_D = Constant(0.0)
# def boundary(x, on_boundary):
#     return on_boundary

bc = DirichletBC(V, u_D, "on_boundary")

# Define the variational problem
u = TrialFunction(V)
v = TestFunction(V)

L = 1.0  # Define the length scale
# f = Constant(1.0)  # Define f(x), here f is constant
f_expr = "sin(x[0]/pi)*sin(x[0]/pi)"
f = Expression(f_expr, degree=2)
a = dot(grad(u), grad(v)) * dx + (1 / L^2) * u * v * dx
L = (1 / L^2) * f * v * dx

# Compute solution
u = FeFunction(V)
lvsolve(a,L, u, bc)
get_array(u)
# # Plot solution
# plot(u)
using PyPlot
FEniCS.Plot(u)
PyPlot.savefig("chatgpt.png")
# plt.show()
