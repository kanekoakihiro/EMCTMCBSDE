using LinearAlgebra, DifferentialEquations
include("./equation.jl")
using .equation
include("./slv.jl")
using .slv
include("./utils_src/grid_1d.jl")
include("./utils_src/difference_operators.jl")
include("./utils_src/multidimensional_operators.jl")
include("./algorithms/method_of_lines.jl")
include("./algorithms/sparse_grids.jl")