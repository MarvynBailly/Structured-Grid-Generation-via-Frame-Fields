"""
    test_greedy_mip_solver.jl

Comprehensive test suite for the GreedyMIPSolver module.

Tests cover:
1. Simple 2×2 system with known solution
2. Larger system with multiple integer variables
3. Sparse system (relevant for mesh generation)
4. Edge cases and validation
"""

using Test
using LinearAlgebra
using SparseArrays

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))

using GreedyMIPSolver

@testset "GreedyMIPSolver Tests" begin
    
    @testset "Test 1: Simple 2×2 System" begin
        println("\n" * "="^60)
        println("Test 1: Simple 2×2 System")
        println("="^60)
        
        # System: 2x₁ + x₂ = 5
        #         x₁ + 3x₂ = 7
        # Continuous solution: x₁ ≈ 1.6, x₂ ≈ 1.8
        # With x₁ integer: x₁ = 2, x₂ = 1
        
        A = [2.0 1.0;
             1.0 3.0]
        b = [5.0, 7.0]
        
        println("Matrix A:")
        display(A)
        println("\nVector b: ", b)
        
        # Solve with first variable as integer
        result = solve_greedy_mip(A, b, 1, verbose=true)
        
        println("\nResult:")
        println("  Solution: ", result.x)
        println("  x[1] (integer): ", result.x[1])
        println("  x[2] (continuous): ", result.x[2])
        println("  Converged: ", result.converged)
        println("  Rounding order: ", result.rounding_order)
        
        # Verify x[1] is integer
        @test isinteger(result.x[1])
        
        # Verify solution satisfies Ax = b reasonably well
        residual = A * result.x - b
        println("  Residual norm: ", norm(residual))
        @test norm(residual) < 1e-5
        
        @test result.converged == true
    end
    
    @testset "Test 2: 4×4 System with 2 Integer Variables" begin
        println("\n" * "="^60)
        println("Test 2: 4×4 System with 2 Integer Variables")
        println("="^60)
        
        # Diagonally dominant system (ensures convergence)
        A = [4.0 1.0 0.5 0.2;
             1.0 5.0 0.3 0.4;
             0.5 0.3 6.0 0.1;
             0.2 0.4 0.1 7.0]
        b = [10.0, 15.0, 20.0, 25.0]
        
        println("Matrix A:")
        display(A)
        println("\nVector b: ", b)
        
        # First two variables must be integers
        result = solve_greedy_mip(A, b, 2, verbose=true)
        
        println("\nResult:")
        println("  Solution: ", result.x)
        println("  x[1] (integer): ", result.x[1])
        println("  x[2] (integer): ", result.x[2])
        println("  Converged: ", result.converged)
        println("  Iterations: ", result.iterations)
        
        # Verify first two variables are integers
        @test isinteger(result.x[1])
        @test isinteger(result.x[2])
        
        # Verify solution quality
        residual = A * result.x - b
        println("  Residual norm: ", norm(residual))
        @test norm(residual) < 1e-4
        
        @test result.converged == true
    end
    
    @testset "Test 3: Sparse System (Mesh-like)" begin
        println("\n" * "="^60)
        println("Test 3: Sparse System (Simulating Mesh Connectivity)")
        println("="^60)
        
        # Simulate a mesh-like sparse system
        n = 10
        A = spdiagm(
            0 => 4.0 * ones(n),      # Main diagonal
            1 => -ones(n-1),          # Super-diagonal
            -1 => -ones(n-1),         # Sub-diagonal
            2 => -0.5 * ones(n-2),    # Second super-diagonal
            -2 => -0.5 * ones(n-2)    # Second sub-diagonal
        )
        A = Matrix(A)  # Convert to dense for this test
        b = collect(1.0:n)
        
        println("System size: ", n, "×", n)
        println("Matrix density: ", count(!iszero, A) / (n*n))
        println("Vector b: ", b)
        
        # First 3 variables must be integers
        k = 3
        result = solve_greedy_mip(A, b, k, verbose=true)
        
        println("\nResult:")
        println("  Integer variables: ", result.x[1:k])
        println("  Continuous variables: ", result.x[k+1:end])
        println("  Converged: ", result.converged)
        println("  Total iterations: ", result.iterations)
        
        # Verify integer constraints
        for i in 1:k
            @test isinteger(result.x[i])
        end
        
        # Verify solution quality
        residual = A * result.x - b
        println("  Residual norm: ", norm(residual))
        @test norm(residual) < 1e-3
    end
    
    @testset "Test 4: Comparison with Direct Rounding" begin
        println("\n" * "="^60)
        println("Test 4: Greedy vs Direct Rounding Comparison")
        println("="^60)
        
        A = [3.0 1.0 0.5;
             1.0 4.0 0.5;
             0.5 0.5 5.0]
        b = [8.0, 12.0, 15.0]
        
        println("Matrix A:")
        display(A)
        println("\nVector b: ", b)
        
        # Method 1: Direct rounding (baseline)
        x_continuous = A \ b
        x_direct = copy(x_continuous)
        x_direct[1:2] = round.(x_direct[1:2])
        # Don't resolve, just use the rounded values
        
        println("\nDirect rounding:")
        println("  Continuous solution: ", x_continuous)
        println("  After rounding: ", x_direct)
        residual_direct = A * x_direct - b
        println("  Residual norm: ", norm(residual_direct))
        
        # Method 2: Greedy rounding
        result = solve_greedy_mip(A, b, 2, verbose=true)
        residual_greedy = A * result.x - b
        
        println("\nGreedy rounding:")
        println("  Solution: ", result.x)
        println("  Residual norm: ", norm(residual_greedy))
        
        # Greedy should give better or equal residual
        println("\nComparison:")
        println("  Direct rounding residual: ", norm(residual_direct))
        println("  Greedy rounding residual: ", norm(residual_greedy))
        println("  Improvement: ", (norm(residual_direct) - norm(residual_greedy)) / norm(residual_direct) * 100, "%")
        
        @test norm(residual_greedy) ≤ norm(residual_direct) + 1e-10
    end
    
    @testset "Test 5: Edge Cases and Validation" begin
        println("\n" * "="^60)
        println("Test 5: Edge Cases and Input Validation")
        println("="^60)
        
        A = [2.0 1.0; 1.0 3.0]
        b = [5.0, 7.0]
        
        # Test: k = n (all variables are integers)
        println("\nTest 5a: All variables integer (k = n)")
        result = solve_greedy_mip(A, b, 2, verbose=false)
        @test all(isinteger.(result.x))
        println("  ✓ All variables are integers")
        
        # Test: k = 1 (only one integer variable)
        println("\nTest 5b: Single integer variable (k = 1)")
        result = solve_greedy_mip(A, b, 1, verbose=false)
        @test isinteger(result.x[1])
        println("  ✓ First variable is integer")
        
        # Test: Invalid k values
        println("\nTest 5c: Invalid input handling")
        @test_throws AssertionError solve_greedy_mip(A, b, 0)  # k too small
        @test_throws AssertionError solve_greedy_mip(A, b, 3)  # k too large
        println("  ✓ Invalid k values properly rejected")
        
        # Test: Non-square matrix
        A_bad = [2.0 1.0 0.5; 1.0 3.0 0.5]
        @test_throws AssertionError solve_greedy_mip(A_bad, b, 1)
        println("  ✓ Non-square matrix rejected")
    end
    
    @testset "Test 6: Tight Convergence Threshold" begin
        println("\n" * "="^60)
        println("Test 6: Different Convergence Thresholds")
        println("="^60)
        
        A = [5.0 1.0 0.5;
             1.0 6.0 0.5;
             0.5 0.5 7.0]
        b = [10.0, 15.0, 20.0]
        
        # Loose threshold
        result_loose = solve_greedy_mip(A, b, 2, τ=1e-3, verbose=false)
        println("Loose threshold (τ=1e-3):")
        println("  Iterations: ", result_loose.iterations)
        println("  Residual: ", norm(A * result_loose.x - b))
        
        # Tight threshold
        result_tight = solve_greedy_mip(A, b, 2, τ=1e-8, verbose=false)
        println("\nTight threshold (τ=1e-8):")
        println("  Iterations: ", result_tight.iterations)
        println("  Residual: ", norm(A * result_tight.x - b))
        
        # Tighter threshold should give better solution but more iterations
        @test norm(A * result_tight.x - b) ≤ norm(A * result_loose.x - b)
        @test result_tight.iterations ≥ result_loose.iterations
        println("\n  ✓ Tighter threshold produces better solution")
    end
end

println("\n" * "="^60)
println("All tests completed!")
println("="^60)
