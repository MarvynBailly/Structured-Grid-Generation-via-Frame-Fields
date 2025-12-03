"""
    example_usage.jl

Simple examples demonstrating the GreedyMIPSolver functionality.
Run this file to see the solver in action with different scenarios.
"""

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

using GreedyMIPSolver
using LinearAlgebra

println("="^70)
println("Greedy MIP Solver - Example Usage")
println("="^70)

# Example 1: Simple 2×2 System
println("\n📊 Example 1: Simple 2×2 System")
println("-"^70)

A = [2.0 1.0;
     1.0 3.0]
b = [5.0, 7.0]

println("System: Ax = b")
println("A = ", A)
println("b = ", b)
println("\nConstraint: x[1] must be integer\n")

result = solve_greedy_mip(A, b, 1, verbose=true)

println("\n✓ Solution found:")
println("  x = ", result.x)
println("  x[1] (integer) = ", result.x[1])
println("  x[2] (continuous) = ", result.x[2])
println("  Residual norm = ", norm(A * result.x - b))

# Example 2: Larger System
println("\n\n📊 Example 2: 5×5 System with 2 Integer Variables")
println("-"^70)

A = [5.0 1.0 0.5 0.2 0.1;
     1.0 6.0 0.5 0.3 0.2;
     0.5 0.5 7.0 0.4 0.1;
     0.2 0.3 0.4 8.0 0.3;
     0.1 0.2 0.1 0.3 9.0]
b = [10.0, 15.0, 20.0, 25.0, 30.0]

println("System size: 5×5")
println("Integer variables: 2 (x[1] and x[2])")
println("Continuous variables: 3 (x[3], x[4], x[5])\n")

result = solve_greedy_mip(A, b, 2, verbose=true)

println("\n✓ Solution found:")
println("  Integer part: ", result.x[1:2])
println("  Continuous part: ", result.x[3:end])
println("  Total iterations: ", result.iterations)
println("  Residual norm = ", norm(A * result.x - b))

# Example 3: Comparison with Direct Rounding
println("\n\n📊 Example 3: Greedy vs Direct Rounding")
println("-"^70)

A = [4.0 1.0 0.5;
     1.0 5.0 0.5;
     0.5 0.5 6.0]
b = [8.0, 12.0, 15.0]

println("Comparing two approaches for rounding x[1]:\n")

# Continuous solution
x_cont = A \ b
println("1️⃣  Continuous solution: ", round.(x_cont, digits=4))

# Direct rounding
x_direct = copy(x_cont)
x_direct[1] = round(x_direct[1])
res_direct = norm(A * x_direct - b)
println("2️⃣  Direct rounding: ", round.(x_direct, digits=4))
println("    Residual: ", round(res_direct, digits=6))

# Greedy rounding
result = solve_greedy_mip(A, b, 1, verbose=false)
res_greedy = norm(A * result.x - b)
println("3️⃣  Greedy rounding: ", round.(result.x, digits=4))
println("    Residual: ", round(res_greedy, digits=6))

improvement = (res_direct - res_greedy) / res_direct * 100
println("\n💡 Greedy method improves residual by ", round(improvement, digits=2), "%")

println("\n" * "="^70)

