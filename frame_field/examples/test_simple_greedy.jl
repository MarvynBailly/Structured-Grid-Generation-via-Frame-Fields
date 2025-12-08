using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include("../src/greedy_solver.jl")

using LinearAlgebra
using SparseArrays

println("=== Simple Greedy Solver Test ===\n")

# Create a simple test system: minimize (1/2)x'Ax - b'x
# where x[1:2] are integers and x[3:4] are reals
#
# Let's use a simple quadratic:
# E = (x1 - 0.7)^2 + (x2 - 1.3)^2 + (x3 - 2.1)^2 + (x4 - 3.6)^2
# Optimal: x1=1, x2=1, x3=2.1, x4=3.6

# The Hessian is 2*I and gradient is 2*(x - optimal)
# So Ax = b becomes: 2*I*x = 2*[0.7, 1.3, 2.1, 3.6]
# Simplifying: x = [0.7, 1.3, 2.1, 3.6]

A = sparse(2.0 * Matrix(I, 4, 4))
b = 2.0 * [0.7, 1.3, 2.1, 3.6]

println("System:")
println("A = ")
display(Matrix(A))
println("\nb = ", b)

println("\n\nSolving with greedy solver (first 2 variables are integers)...")
x = greedy_mip_solver(A, b, 2; tau=1e-6, max_iterations=100, verbose=true)

println("\n\nResults:")
println("x = ", x)
println("Expected: [1, 1, 2.1, 3.6] (approx)")
println("Error: ", norm(x - [1.0, 1.0, 2.1, 3.6]))

println("\n=== Test Complete ===")
