using FEniCS

# Mesh and function space
mesh = UnitSquareMesh(8,8)
V = FunctionSpace(mesh, "P", 1)

# Boundary conditions
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

function boundary(x, on_boundary)
    return on_boundary
end

bc = DirichletBC(V,u_D,"on_boundary")

# Variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

u = FeFunction(V)
lvsolve(a,L,u,bc)
get_array(L)
get_array(U)
import PyPlot
vtkfile = File("poisson/solution.pvd")
vtkfile << U.pyobject #exports the solution to a vtkfile
FEniCS.Plot(U)
FEniCS.Plot(mesh)
PyPlot.savefig("fenics2.2.1.png")